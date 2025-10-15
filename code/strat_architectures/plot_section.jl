# Plot a simple ALCAP model
module Plot_ck

using CarboKitten
using Unitful
using GLMakie
using CarboKitten.Visualization
using CarboKitten.Export: read_slice

function make_summary_plots()
    fig = summary_plot("data/sinusoid.h5")
    save("figs/sinusoid.png", fig)
    fig = summary_plot("data/miller_2020.h5")
    save("figs/miller.png", fig)
end

end

Plot_ck.make_summary_plots()