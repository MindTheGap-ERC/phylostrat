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

#### Print settings ####
# Figures are designed at final print size (116 mm wide, matching the R figures in code/make_figs.R).
# Figure sizes are given in pt (1 pt = 1/72 in), so font sizes and line widths below
# correspond to their printed size. Raster output is written at `DPI`.
const MM = 72 / 25.4                 # pt per mm
const TEXT_WIDTH = 116 * MM          # printing width [pt]
const DPI = 600
const PX_PER_UNIT = DPI / 72

const LW_MAIN = 1.0                  # data lines (sea level curve)
const LW_COEVAL = 0.8                # time lines / coeval lines
const LW_SAMPLING = 0.8              # sampling position
# unconformities: CarboKitten draws them with linewidth = 1, which was designed for the
# original 1200 px wide figure. Scale to keep their original size relative to the figure.
const LW_UNCONFORMITY = TEXT_WIDTH / 1200
const LW_SEALEVEL_PROFILE = 2.5       # end sea level in sediment profiles (half of CarboKitten's 5)
const PANEL_LABEL_SIZE = 10           # as ggpubr panel letters in R

const PRINT_THEME = Theme(
    fontsize = 8,
    figure_padding = (2, 6, 2, 2),
    backgroundcolor = :white,
    Axis = (
        titlesize = 9,                # as plot.title in R
        titlefont = :bold,
        titlealign = :left,
        xlabelsize = 8,
        ylabelsize = 8,
        xticklabelsize = 7,
        yticklabelsize = 7,
        spinewidth = 0.6,
        xtickwidth = 0.6,
        ytickwidth = 0.6,
        xticksize = 3,
        yticksize = 3,
        xgridvisible = false,
        ygridvisible = false,
        xlabelpadding = 2,
        ylabelpadding = 2,
    ),
    Colorbar = (
        labelsize = 8,
        ticklabelsize = 7,
        size = 8,
        spinewidth = 0.6,
        tickwidth = 0.6,
        ticksize = 3,
    ),
    Lines = (linewidth = LW_MAIN,),
)

panel_label!(pos, txt) = Label(pos, txt; fontsize = PANEL_LABEL_SIZE, font = :bold,
    padding = (0, 6, 4, 0), halign = :right)

save_print(fname, fig) = save(fname, fig; px_per_unit = PX_PER_UNIT)

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
        linewidth = LW_COEVAL, color = :black, linestyle = :dash)
    vlines!(ax, sampling_position |> in_units_of(u"km"),
        linewidth = LW_SAMPLING, linestyle=:solid, color=:black)
    #contourf!(ax, xkm, tmyr, wd;
    #    levels=[0.0, 10000.0], colormap=Reverse(:grays))
    #contour!(ax, xkm, tmyr, wd;
    #    levels=[0], color=:black, linewidth=2)
    return ft
end

const FACIES_TICKS = (1:3, ["Euphotic", "Oligophotic", "Aphotic"])
const FACIES_LABEL = "Dominant carbonate factory"

function my_profile!(ax::Axis, header::Header, data::DataSlice, time_lines, sampling_position)
    sediment_profile!(ax, header, data, show_coeval_lines = false)
    # unconformities are the only `linesegments` plot drawn by `sediment_profile!`,
    # the end sea level is the only `hlines` plot (CarboKitten default linewidth = 5)
    for p in ax.scene.plots
        p isa LineSegments && (p.linewidth = LW_UNCONFORMITY)
        p isa HLines && (p.linewidth = LW_SEALEVEL_PROFILE)
    end
    coeval_lines!(ax, header, data, [time_lines...];
        linewidth = LW_COEVAL, color = :black, linestyle = :dash)
    vlines!(ax, sampling_position |> in_units_of(u"km"),
        linewidth = LW_SAMPLING, linestyle = :solid, color = :black)
    ax.xlabel = pos_label
    ax.ylabel = depth_label
    ax.title = "Sediment profile" # set after plotting, CarboKitten overwrites it
    return ax
end

facies_colorbar!(pos, ft; width = Relative(0.6)) = Colorbar(pos, ft; vertical = false, flipaxis = false,
    ticks = FACIES_TICKS, label = FACIES_LABEL,
    width = width, tellheight = true)

function plot_scenario_comparison()
    fig = Figure(size = (TEXT_WIDTH, 115 * MM))

    # column headers
    Label(fig[0, 1], "Scenario 1: Miller et al. (2020) sea level",
        fontsize = 8.5, font = :bold, tellwidth = false)
    Label(fig[0, 2], "Scenario 2: sinusoidal sea level",
        fontsize = 8.5, font = :bold, tellwidth = false)

    header, data = read_slice("data/strat/miller_2020.h5", :profile)
    ax_miller_profile = Axis(fig[1, 1], title = "Sediment profile")
    my_profile!(ax_miller_profile, header, data, miller_time_lines, sampling_position)
    panel_label!(fig[1, 1, TopLeft()], "A")

    ax_miller_wheeler = Axis(fig[2, 1], title = "Wheeler diagram")
    my_dominant_facies!(ax_miller_wheeler, header, data, miller_time_lines, sampling_position)
    ax_miller_wheeler.xlabel = pos_label
    ax_miller_wheeler.ylabel = time_label
    panel_label!(fig[2, 1, TopLeft()], "C")

    header, data = read_slice("data/strat/sinusoid.h5", :profile)
    ax_sinusoid_profile = Axis(fig[1, 2], title = "Sediment profile")
    my_profile!(ax_sinusoid_profile, header, data, sinusoid_time_lines, sampling_position)
    panel_label!(fig[1, 2, TopLeft()], "B")

    ax_sinusoid_wheeler = Axis(fig[2, 2], title = "Wheeler diagram")
    ft = my_dominant_facies!(ax_sinusoid_wheeler, header, data, sinusoid_time_lines,
        sampling_position)
    ax_sinusoid_wheeler.xlabel = pos_label
    ax_sinusoid_wheeler.ylabel = time_label
    panel_label!(fig[2, 2, TopLeft()], "D")

    linkxaxes!(ax_miller_profile, ax_miller_wheeler)
    linkxaxes!(ax_sinusoid_profile, ax_sinusoid_wheeler)
    facies_colorbar!(fig[3, 1:2], ft)

    rowgap!(fig.layout, 6)
    colgap!(fig.layout, 10)
    save_print("figs/ms/scenario_comparison.png", fig)
end

function my_summary_plot(input::String, fig_name::String, time_lines, sampling_position, case)
    header, data = read_slice(input, :profile)

    fig = Figure(size = (TEXT_WIDTH, 112 * MM))

    Label(fig[0, 1:3], case, fontsize = 8.5, font = :bold, tellwidth = false)

    ax_profile = Axis(fig[1, 1], title = "Sediment profile")
    my_profile!(ax_profile, header, data, time_lines, sampling_position)
    panel_label!(fig[1, 1, TopLeft()], "A")

    ax_production = Axis(fig[1, 2], title = "Production profile")
    h5open(input, "r") do fid
        production_curve!(ax_production, fid["input"], max_depth=-50.0u"m")
    end
    ax_production.xlabel = "Production [m/Myr]"
    ax_production.ylabel = "Water depth [m]"
    ax_production.title = "Production profile" # set after plotting, CarboKitten overwrites it
    panel_label!(fig[1, 2, TopLeft()], "B")

    ax_wheeler = Axis(fig[2, 1], title = "Wheeler diagram")
    df = my_dominant_facies!(ax_wheeler, header, data, time_lines, sampling_position)
    ax_wheeler.xlabel = pos_label
    ax_wheeler.ylabel = time_label
    panel_label!(fig[2, 1, TopLeft()], "C")

    ax_ws = Axis(fig[2, 2],
        title = "Eustatic sea level",
        xlabel = "Sea level [m]",
        ylabel = time_label,
        xticks = LinearTicks(2),
        limits = (nothing, (header.axes.t[1] |> in_units_of(u"Myr"), header.axes.t[end] |> in_units_of(u"Myr"))))
    lines!(ax_ws, header.sea_level |> in_units_of(u"m"), header.axes.t |> in_units_of(u"Myr"),
        linewidth = 0.6, color = :blue) # thin: high-frequency signal
    hlines!(ax_ws, time_lines |> in_units_of(u"Myr"),
        linewidth = LW_COEVAL, color = :black, linestyle = :dash)
    panel_label!(fig[2, 2, TopLeft()], "D")

    linkxaxes!(ax_profile, ax_wheeler)
    linkyaxes!(ax_wheeler, ax_ws)
    # shared time axis: drop redundant labels on the sea level panel
    hideydecorations!(ax_ws, ticks = false, grid = false)

    facies_colorbar!(fig[3, 1], df; width = Relative(0.9))

    colsize!(fig.layout, 1, Relative(0.68))
    rowgap!(fig.layout, 6)
    colgap!(fig.layout, 8)
    save_print(fig_name, fig)
end

function miller_summary()
    my_summary_plot("data/strat/miller_2020.h5", "figs/sm/miller_summary.png",
    miller_time_lines, sampling_position, "Scenario 1: Miller et al. (2020) sea level")
end

function sinusoidal_summary()
    my_summary_plot("data/strat/sinusoid.h5", "figs/sm/sinusoid_summary.png",
    sinusoid_time_lines, sampling_position, "Scenario 2: sinusoidal sea level")
end

set_theme!(PRINT_THEME)

end

Plot.plot_scenario_comparison()
Plot.miller_summary()
Plot.sinusoidal_summary()
