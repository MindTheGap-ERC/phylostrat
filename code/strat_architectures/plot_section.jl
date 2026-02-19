# Plot a simple ALCAP model
module Plot_ck

using CarboKitten
using Unitful
using GLMakie
using CarboKitten.Visualization
using CarboKitten.Export: read_header, read_volume, read_slice, group_datasets

function make_summary_plots()
    fig = summary_plot("data/strat/sinusoid.h5")
    save("figs/sinusoid.png", fig)
    fig = summary_plot("data/strat/miller_2020.h5")
    save("figs/miller.png", fig)
end

function combined_strat_slice()
    fig = Figure(size = (1200, 1000), backgroundcolor = :gray80)
    header, data = read_slice("data/strat/sinusoid.h5", :profile)
    ax = Axis(fig[1, 1])
    sediment_profile!(ax, header, data)
    ax = Axis(fig[1, 2])
    header, data = read_slice("data/strat/miller_2020.h5", :profile)
    sediment_profile!(ax, header, data)
    save("figs/comb_strat.png", fig)
end


end

Plot_ck.combined_strat_slice()
Plot_ck.make_summary_plots()