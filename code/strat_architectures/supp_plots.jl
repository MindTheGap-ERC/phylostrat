import CarboKitten.Visualization: wheeler_diagram, wheeler_diagram!
using CarboKitten.Export: Header, Data, DataSlice, read_data, read_slice
using CarboKitten.Utility: in_units_of
using GLMakie
using Unitful
using CarboKitten.BoundaryTrait
using CarboKitten.Stencil: convolution


const na = [CartesianIndex()]

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

function my_dominant_facies!(ax::Axis, header::Header, data::DataSlice;
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
        nan_color = :transparent)
    #contourf!(ax, xkm, tmyr, wd;
    #    levels=[0.0, 10000.0], colormap=Reverse(:grays))
    #contour!(ax, xkm, tmyr, wd;
    #    levels=[0], color=:black, linewidth=2)
    return ft
end

    fig = Figure(size = (1200, 1000), backgroundcolor = :gray80)
    header, data = read_slice("data/strat/sinusoid.h5", :profile)
    ax = Axis(fig[1, 1])

    my_dominant_facies!(ax, header, data)

    save("figs/dom_facies.png", fig)