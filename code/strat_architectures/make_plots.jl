module Plot
import CarboKitten.Visualization: wheeler_diagram, wheeler_diagram!, sediment_profile!, coeval_lines!, production_curve!
using CarboKitten.Export: Header, Data, DataSlice, read_data, read_slice
using CarboKitten.Utility: in_units_of
using GLMakie
using Unitful
using CarboKitten.BoundaryTrait
using CarboKitten.Stencil: convolution
using HDF5

const na = [CartesianIndex()]

pos_label = "Distance from shore [km]"
time_label = "Elapsed model time [Myr]"
depth_label = "Depth [m]"
sampling_position = 3u"km"
miller_time_lines = -4.0u"Myr":1.0u"Myr":-1.0u"Myr"
sinusoid_time_lines = 1.0u"Myr":1.0u"Myr":4.0u"Myr"

elevation(h::Header, d::DataSlice) =
    let bl = h.initial_topography[d.slice..., na],
        sr = (h.axes.t[end] - h.axes.t[1]) * h.subsidence_rate

        bl .+ d.sediment_thickness[:, 2:end] .- sr
    end

water_depth(header::Header, data::DataSlice) =
    let h = elevation(header, data),
        wi = data.write_interval,
        s = header.subsidence_rate .* (header.axes.t[wi:wi:end] .- header.axes.t[end]),
        l = header.sea_level[wi:wi:end]

        h .- (s.+l)[na, :]
    end

const Rate = typeof(1.0u"m/Myr")

function my_dominant_facies!(ax::Axis, header::Header, data::DataSlice, time_lines, sampling_position;
    smooth_size::NTuple{2,Int}=(1,1),
    colors=Makie.wong_colors())
    n_facies = size(data.production)[1]
    colormax(d) = getindex.(argmax(d; dims=1)[1, :, :], 1)

    wi = data.write_interval

    dominant_facies = colormax(data.deposition)
    blur = convolution(Shelf, ones(Float64, smooth_size...) ./ *(smooth_size...))
    #wd = zeros(Float64, length(header.axes.x), length(header.axes.t[wi:wi:end]))
    #blur(water_depth(header, data) / u"m", wd)
    wd = water_depth(header, data) |> in_units_of(u"m")
    ax.ylabel = "time [Myr]"
    ax.xlabel = "position [km]"

    dominant_facies = Matrix{Union{Missing, Int}}(dominant_facies)
    dominant_facies[ wd .> 0] .= missing

    xkm = header.axes.x |> in_units_of(u"km")
    tmyr = header.axes.t[wi:wi:end] |> in_units_of(u"Myr")
    ft = heatmap!(ax, xkm, tmyr, dominant_facies;
        colormap=cgrad(colors[1:n_facies], n_facies, categorical=true),
        colorrange=(0.5, n_facies + 0.5),
        nan_color = :white)

    hlines!(ax, time_lines |> in_units_of(u"Myr"),
        linewidth = 2, color = :black, linestyle = :dash)
    vlines!(ax, sampling_position |> in_units_of(u"km"),
        linestyle=:solid, color=:black)
    #contourf!(ax, xkm, tmyr, wd;
    #    levels=[0.0, 10000.0], colormap=Reverse(:grays))
    #contour!(ax, xkm, tmyr, wd;
    #    levels=[0], color=:black, linewidth=2)
    return ft
end

function plot_scenario_comparison()
fig = Figure(size = (1200, 1000), backgroundcolor = :gray80)
header, data = read_slice("data/strat/miller_2020.h5", :profile)
pos_max = header.axes.x[end] |> in_units_of(u"km")
pos_min = header.axes.x[1] |> in_units_of(u"km")
ax_miller_profile = Axis(fig[1, 1])
sediment_profile!(ax_miller_profile, header, data, show_coeval_lines = false)
coeval_lines!(ax_miller_profile, header, data, [miller_time_lines...];
    linewidth = 2, color = :black, linestyle = :dash)
ax_miller_profile.xlabel = pos_label
ax_miller_profile.ylabel = depth_label
ax_miller_profile.title = "Sediment profile scenario 1 (Miller et al. (2020) sea level curve)"
vlines!(ax_miller_profile, sampling_position |> in_units_of(u"km"),
    linestyle=:solid, color=:black)
Label(fig[1, 1, TopLeft()], "A", fontsize = 30)


ax_miller_wheeler = Axis(fig[2,1])
my_dominant_facies!(ax_miller_wheeler, header, data, miller_time_lines, sampling_position)
ax_miller_wheeler.xlabel = pos_label
ax_miller_wheeler.ylabel = time_label
ax_miller_wheeler.title = "Wheeler diagram scenario 1 (Miller et al. (2020) sea level curve)"
Label(fig[2, 1, TopLeft()], "C", fontsize = 30)

ax_sinusoid_profile = Axis(fig[1, 2])
header, data = read_slice("data/strat/sinusoid.h5", :profile)
pos_max = header.axes.x[end] |> in_units_of(u"km")
pos_min = header.axes.x[1] |> in_units_of(u"km")
sediment_profile!(ax_sinusoid_profile, header, data, show_coeval_lines = false)
coeval_lines!(ax_sinusoid_profile, header, data, [sinusoid_time_lines...];
    linewidth = 2, color = :black, linestyle = :dash)
ax_sinusoid_profile.xlabel = pos_label
ax_sinusoid_profile.ylabel = depth_label
ax_sinusoid_profile.title = "Sediment profile scenario 2 (sinusoidal sea level curve)"
vlines!(ax_sinusoid_profile, sampling_position |> in_units_of(u"km"),
linestyle=:solid, color=:black)
Label(fig[1, 2, TopLeft()], "B", fontsize = 30)


ax_sinusoid_wheeler = Axis(fig[2,2])
ft = my_dominant_facies!(ax_sinusoid_wheeler, header, data, sinusoid_time_lines,
    sampling_position)
ax_sinusoid_wheeler.xlabel = pos_label
ax_sinusoid_wheeler.ylabel = time_label
ax_sinusoid_wheeler.title = "Wheeler diagram scenario 2 (sinusoidal sea level curve)"
Label(fig[2, 2, TopLeft()], "D", fontsize = 30)

linkxaxes!(ax_miller_profile, ax_miller_wheeler)
linkxaxes!(ax_sinusoid_profile, ax_sinusoid_wheeler)
Colorbar(fig[3,1:2], ft, vertical = false, ticks=(1:3, ["Euphotic", "Oligophotic", "Aphotic"]), label="Dominant carbonate factory")
save("figs/ms/scenario_comparison.png", fig)
end

function my_summary_plot(input::String, fig_name::String, time_lines, sampling_position, case)
    header, data = read_slice(input, :profile)

    fig = Figure(size = (1200, 1000), backgroundcolor = :gray80)

    ax_profile = Axis(fig[1:2, 1:2])
    sediment_profile!(ax_profile, header, data, show_coeval_lines = false)
    coeval_lines!(ax_profile, header, data, [time_lines...];
    linewidth = 2, color = :black, linestyle = :dash)
    vlines!(ax_profile, sampling_position |> in_units_of(u"km"),
    linestyle=:solid, color=:black)
    ax_profile.xlabel = pos_label
    ax_profile.ylabel = depth_label
    ax_profile.title = "Sediment profile for the $(case)"
    Label(fig[1:2, 1:2, TopLeft()], "A", fontsize = 30)

    ax_wheeler = Axis(fig[3:4, 1:2])
    df = my_dominant_facies!(ax_wheeler, header, data, time_lines, sampling_position)
    ax_wheeler.xlabel = pos_label
    ax_wheeler.ylabel = time_label
    ax_wheeler.title = "Wheeler diagram for the $(case)"
    Label(fig[3:4, 1:2, TopLeft()], "C", fontsize = 30)

    ax_ws = Axis(fig[3:4, 3],
        title = "Eustatic sea level",
        xlabel = "Sea level [m]",
        limits = (nothing, (header.axes.t[1] |> in_units_of(u"Myr"), header.axes.t[end] |> in_units_of(u"Myr"))))
    lines!(ax_ws, header.sea_level |> in_units_of(u"m"), header.axes.t |> in_units_of(u"Myr"))
    hlines!(ax_ws, time_lines |> in_units_of(u"Myr"),
        linewidth = 2, color = :black, linestyle = :dash)
    ax_ws.ylabel = time_label
    Label(fig[3:4, 3, TopLeft()], "D", fontsize = 30)

    ax_production = Axis(fig[1:2, 3])
    h5open(input, "r") do fid
        production_curve!(ax_production, fid["input"], max_depth=-50.0u"m")
    end
    ax_production.xlabel = "Production [m/Myr]"
    ax_production.ylabel = "Water depth [m]"
    ax_production.title = "Production profile"
    Label(fig[1:2, 3, TopLeft()], "B", fontsize = 30)
    

    linkxaxes!(ax_profile, ax_wheeler)
    linkyaxes!(ax_wheeler, ax_ws)

    Colorbar(fig[5,1:2], df, vertical = false, ticks=(1:3, ["Euphotic", "Oligophotic", "Aphotic"]), label="Dominant carbonate factory")



    save(fig_name, fig)
end

function miller_summary()
    my_summary_plot("data/strat/miller_2020.h5", "figs/sm/miller_summary.png",
    miller_time_lines, sampling_position, "Miller et al. (2020) sea level curve (scenario 1)")
end

function sinusoidal_summary()
    my_summary_plot("data/strat/sinusoid.h5", "figs/sm/sinusoid_summary.png",
    sinusoid_time_lines, sampling_position, "sinusoidal sea level curve (scenario 2)")
end

end

Plot.plot_scenario_comparison()
Plot.miller_summary()
Plot.sinusoidal_summary()
