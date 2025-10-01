# Plot a simple ALCAP model
module Plot_ck

using CarboKitten
using Unitful
using GLMakie
using CarboKitten.Visualization
using CarboKitten.Export: read_slice

    function plot_section(path::String, file::String)
        header, profile = read_slice(path * file, :profile)
        fig = sediment_profile(header, profile)
        save("data/$(file).png", fig)
    end

end
