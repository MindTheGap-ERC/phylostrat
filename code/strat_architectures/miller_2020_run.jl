using CarboKitten
using CarboKitten.DataSets: miller_2020
using CarboKitten.Export: read_slice, data_export, CSV, read_volume
using Unitful
using DataFrames
using Interpolations
using Random

Random.seed!(475307)
const TAG = "miller_2020"

include("constants.jl")

function sea_level_miller()
    df = miller_2020()
    sort!(df, [:time])
    return linear_interpolation(
        df.time,
        df.sealevel)
end

const INPUT_MILLER = ALCAP.Input(
    tag="$TAG",
    box=Box{Coast}(grid_size=(GRID_SIZE_X, GRID_SIZE_Y), phys_scale=PHYS_SCALE),
    time=TimeProperties(
        Δt=DELTA_T,
        steps=STEPS_MAIN,
        t0=-5.0001u"Myr"),
    ca_interval=CA_INTERVAL,
    initial_topography=INIT_TOPO,
    output=Dict(
            :profile => OutputSpec(slice = (:, div(GRID_SIZE_Y, 2)), write_interval = PROFILE_WRITE_INTERVAL),
            :topography => OutputSpec(slice = (:, :), write_interval = TOPOGRAPHY_WRITE_INTERVAL)),
    sea_level=sea_level_miller(),
    subsidence_rate=SUBSIDENCE_RATE,
    disintegration_rate=DISINTGRATION_RATE,
    insolation=INSOLATION,
    sediment_buffer_size=SEDIMENT_BUFFER_SIZE,
    depositional_resolution=DEPOSITIONAL_RESOLUTION,
    cementation_time = CEMENTATION_TIME,
    facies=FACIES)

run_model(Model{ALCAP}, INPUT_MILLER, "data/strat/$(TAG).h5")


function export_files()
    header, profile = read_slice("data/strat/$(TAG).h5", :profile)
    columns = [profile[i] for i in [21, 41, 61 , 81, 101, 121, 141, 161]]
    data_export(
        CSV(
            :stratigraphic_column => "data/strat/$(TAG)_sc.csv",
            :age_depth_model      => "data/strat/$(TAG)_adm.csv",
            :metadata => "data/strat/$(TAG)_metadata.toml",
            :water_depth => "data/strat/$(TAG)_wd.csv",
            :sediment_accumulation_curve => "data/strat/$(TAG)_sac.csv"),
        header,
        columns)
end

export_files()