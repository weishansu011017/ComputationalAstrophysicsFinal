include("read_dumpfile.jl")
using Makie
using GLMakie
using Observables
using GeometryBasics
using Statistics
using Glob
using Dates
using Printf

function get_files()
    if Sys.iswindows()
        if isempty(ARGS)
            error("On Windows, please specify the folder containing HDF5 files as an argument.\nExample: julia visualization_2Dparticles_interactive.jl myfolder")
        end
        folder = ARGS[1]
        return sort(glob("*.h5", folder)) 
    else
        if isempty(ARGS)
            error("Please provide one or more HDF5 files as arguments.\nExample: julia visualization_2Dparticles_interactive.jl Sim_*.h5")
        end
        return ARGS
    end
end

function main()
    # Two panel: The upper is the 2D visualization, the lower is the error-time diagram(fixed)
    ## Variable of upper panel
    xulabel = "x [pc]"
    yulabel = "y [pc]"
    SimulationTag = Observable("")
    ## Trigger of upper panel
    xytitle = Observable((Float32[], Float32[], @sprintf("t = %.3f [Myr], Error = %.1f %% (%s)", 0.0, 0.0, SimulationTag[])))
    ## Split trigger
    points = map(xytitle) do (x, y, _)
        Point{2, Float32}.(x, y)
    end
    title_obs = @lift($xytitle[3])

    # Read file
    files = get_files()
    testdata = read_dumpfile(files[1])
    SimulationTag[] = testdata.params["SimulationTag"]

    # Initalize array
    utime = testdata.params["utime"]
    time_array = zeros(Float32, length(files))
    ## Upper
    x_array = zeros(Float32, testdata.params["N"], length(files))
    y_array = zeros(Float32, testdata.params["N"], length(files))
    ## Lower
    KE_array = zeros(Float32, length(files))
    U_array = zeros(Float32, length(files))

    for (i,file) in enumerate(files)
        data = read_dumpfile(file)
        time_array[i] = data.params["t"]
        
        x = data.Table[!,"x"]
        y = data.Table[!,"y"]
        vx = data.Table[!,"vx"]
        vy = data.Table[!,"vy"]
        m = data.Table[!,"m"]
        KE_array[i] = sum(0.5 .* m .* (vx.^2 .+ vy.^2))
        U_array[i] = data.params["Utot"]

        x_array[:,i] = x
        y_array[:,i] = y

    end
    time_array .*= (utime/3.15576e13)
    ## Lower panel
    E_array = KE_array .+ U_array
    E1 = E_array[1]
    Error_array = 100 .* abs.((E_array .- E1) ./ E1)

    # Initalize figure
    should_exit = Observable(false)
    fig = Figure(size=(800,1200))
    # Upper panel
    axu = Axis(fig[1:2,1], xlabel = xulabel, ylabel = yulabel, limits = ((-0.25, 0.25), (-0.25, 0.25)),  title = title_obs, xlabelsize = 20, ylabelsize = 20, titlesize = 22, xticklabelsize = 14, yticklabelsize = 14)
    # Lower panel
    axd = Axis(fig[3,1], xlabel = "t [Myr]", ylabel = "Error [%]" , xlabelsize = 20, ylabelsize = 20, titlesize = 22, xticklabelsize = 14, yticklabelsize = 14)
    # Draw lower panel first
    lines!(axd,  time_array, Error_array, color=:black, label="Error")
    Legend(fig[3, 2], axd, labelsize = 14)

    # First upper panel
    xytitle[] = (x_array[:,1], y_array[:,1], @sprintf("t = %.3f [Myr], Error = %.1f %% (%s)", 0.0, 0.0, SimulationTag[]))
    scatter!(axu, points, markersize=5)
    
    # interactive
    frame_idx = Observable(1) 
    is_playing = Observable(false)
    slider_widget = Slider(fig[4, 1], range = 1:length(files), startvalue = 1)
    slider_label = Label(fig[4, 2], string(slider_widget.value[]))
    on(slider_widget.value) do val
        slider_label.text[] = string("Frame: ", val)
    end
    play_btn = Button(fig[4, 3], label = "Play")

    range_slider = Slider(fig[1:2, 2], range = 0.1:0.01:10.0, startvalue = 0.25, width = Relative(1.0), horizontal = false)
    Label(fig[1:2, 3], "Axis range", rotation = π/2, halign = :center, valign = :center)

    on(slider_widget.value) do val
        frame_idx[] = val             
        slider_label.text[] = string("Frame: ", val)
    end

    on(range_slider.value) do r
        axu.limits[] = ((-r, r), (-r, r))
    end

    on(frame_idx) do i
        title = @sprintf("t = %.3f [Myr], Error = %.1f %% (%s)", time_array[i], Error_array[i], SimulationTag[])
        xytitle[] = (x_array[:,i], y_array[:,i], title)
        slider_label.text[] = string("Frame: ", i)
    end
    on(play_btn.clicks) do _
        is_playing[] = !is_playing[]
        play_btn.label[] = is_playing[] ? "Pause" : "Play"
    end
    
    @async begin
        while !should_exit[]
            if is_playing[]
                i = frame_idx[]
                if i < length(files)
                    frame_idx[] = i + 1
                else
                    is_playing[] = false       
                    play_btn.label[] = "Play"  
                end
            end
            sleep(0.1)
        end
    end
    
    screen = display(fig)
    @info "The interactive screen is prepared."
    wait(screen)
    should_exit[] = true 
end
main()