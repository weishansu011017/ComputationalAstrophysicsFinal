include("read_dumpfile.jl")
using GLMakie
using Observables
using Glob
using GeometryBasics
using Statistics



function main()
    xlabel = Observable{AbstractString}("x [pc]")
    ylabel = Observable{AbstractString}("y [pc]")
    xytitle = Observable((Float32[], Float32[], "t = 0.0, , δE = 0.0"))

    points = @lift(Point{2, Float32}.($xytitle[1], $xytitle[2]))
    title_obs = @lift($xytitle[3])

    files = glob("*.h5", ".")
    testdata = read_dumpfile(files[1])

    time_array = zeros(Float32, length(files))
    x_array = zeros(Float32, testdata.params["N"], length(files))
    y_array = zeros(Float32, testdata.params["N"], length(files))
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
    E_array = KE_array .+ U_array
    E1 = E_array[1]
    δE_array = abs.((E_array .- E1) ./ E1)

    fig = Figure(size=(800,800))
    # ax = Axis(fig[1,1], xlabel = "t [code unit]", ylabel = "Energy" )
    # lines!(ax,  time_array, δE_array, color=:black, label="δE/E0")
    # Legend(fig[1, 2], ax)
    # display(fig)
    # readline()

    ax = Axis(fig[1,1], xlabel = xlabel, ylabel = ylabel, limits = ((-0.25, 0.25), (-0.25, 0.25)),  title = title_obs)
    sc = scatter!(ax, points, markersize=5)
    display(fig)


    println("Enter to start the animation")
    readline()
    for i in eachindex(time_array)
        title = "t = $(time_array[i]), δE = $(δE_array[i])"
        xytitle[] = (x_array[:,i], y_array[:,i], title)
        sleep(0.1)
    end
    println("End animation")
    readline()




end
main()