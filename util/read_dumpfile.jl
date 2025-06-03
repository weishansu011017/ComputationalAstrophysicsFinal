using HDF5
using DataFrames

struct ParticlesTable
    Table::DataFrame
    params::Dict{String, Any}
end

function read_dumpfile(filepath::AbstractString)::ParticlesTable
    h5open(filepath, "r") do f
        tbl = f["Table"]
        particle_index = read(tbl["particle_index"])
        x  = Float64.(read(tbl["x"]))
        y  = Float64.(read(tbl["y"]))
        z  = Float64.(read(tbl["z"]))
        vx = Float64.(read(tbl["vx"]))
        vy = Float64.(read(tbl["vy"]))
        vz = Float64.(read(tbl["vz"]))
        m  = Float64.(read(tbl["m"]))
        h  = Float64.(read(tbl["h"]))
        dt = Float64.(read(tbl["dt"]))

        df = DataFrame(
            particle_index = particle_index,
            x = x, y = y, z = z,
            vx = vx, vy = vy, vz = vz,
            m = m, h = h, dt = dt,
        )

        param_group = f["params"]
        params = Dict{String, Any}()
        for k in keys(param_group)
            val = read(param_group[k])
            params[string(k)] = val
        end

        return ParticlesTable(df, params)
    end
end

function transfer_cgs!(pt::ParticlesTable; year=true)
    udist = pt.params["udist"]
    umass = pt.params["umass"]
    utime = pt.params["utime"]
    uv = udist / utime

    df = pt.Table
    df.x .*= udist
    df.y .*= udist
    df.z .*= udist
    df.vx .*= uv
    df.vy .*= uv
    df.vz .*= uv
    df.m .*= umass
    df.h .*= udist
    df.dt .*= utime

    pt.params["Mtot"] *= umass
    pt.params["t"] *= utime

    if year
        sec2yr = 31556926.0
        pt.params["t"] /= sec2yr
        df.dt ./= sec2yr
    end
end