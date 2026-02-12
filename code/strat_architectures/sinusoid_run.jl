using CarboKitten
using Unitful
using CarboKitten.Export: read_slice, data_export, CSV, read_volume
using DataFrames
using Random

Random.seed!(475308)
const TAG = "sinusoid"

include("constants.jl")

const INPUT_SINUSOID = ALCAP.Input(
    tag="$TAG",
    box=Box{Coast}(grid_size=(GRID_SIZE_X, GRID_SIZE_Y), phys_scale=PHYS_SCALE),
    time=TimeProperties(
        Δt=DELTA_T,
        steps=STEPS_MAIN),
    ca_interval=CA_INTERVAL,
    initial_topography=INIT_TOPO,
    output=Dict(
            :profile => OutputSpec(slice = (:, div(GRID_SIZE_Y, 2)), write_interval = PROFILE_WRITE_INTERVAL),
            :topography => OutputSpec(slice = (:, :), write_interval = TOPOGRAPHY_WRITE_INTERVAL)),
    sea_level=MAIN_SL_SINUSOID,
    subsidence_rate=SUBSIDENCE_RATE,
    disintegration_rate=DISINTGRATION_RATE,
    insolation=INSOLATION,
    sediment_buffer_size=SEDIMENT_BUFFER_SIZE,
    depositional_resolution=DEPOSITIONAL_RESOLUTION,
    cementation_time = CEMENTATION_TIME,
    facies=FACIES)

run_model(Model{ALCAP}, INPUT_SINUSOID, "data/strat/$(TAG).h5")

function export_files()
    header, profile = read_slice("data/strat/$(TAG).h5", :profile)
    columns = [profile[i] for i in [21, 41, 61, 81, 101, 121, 141, 161]]
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