# Run a simple ALCAP model
module Sinusoid

using CarboKitten
using Unitful
using CarboKitten.Export: read_slice, data_export, CSV
using DataFrames

const TAG = "sinusoid"

const FACIES = [
    ALCAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate=500u"m/Myr",
        extinction_coefficient=0.8u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=50.0u"m/yr"),
    ALCAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate=400u"m/Myr",
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=25.0u"m/yr"),
    ALCAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate=100u"m/Myr",
        extinction_coefficient=0.005u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=12.5u"m/yr")
]

# matching the parameters of Holland & Patzkowsky 1999
const PERIOD1 = 2.0u"Myr"
const PERIOD2 = 0.1u"Myr"
const AMPLITUDE1 = 50.0u"m"
const AMPLITUDE2 = 10.0u"m"

const INPUT = ALCAP.Input(
    tag="$TAG",
    box=Box{Coast}(grid_size=(150, 50), phys_scale=150.0u"m"),
    time=TimeProperties(
        Δt=0.0001u"Myr",
        steps=60000),
    ca_interval=1,
    initial_topography=(x, y) -> -x / 300.0,
    output=Dict(
            :profile => OutputSpec(slice = (:, 50), write_interval = 10)),
    sea_level=t -> AMPLITUDE1 * sin(2π * t / PERIOD1) + AMPLITUDE2 * sin(2π * t / PERIOD2),
    subsidence_rate=30.0u"m/Myr",
    disintegration_rate=50.0u"m/Myr",
    insolation=400.0u"W/m^2",
    sediment_buffer_size=50,
    depositional_resolution=0.5u"m",
    facies=FACIES)

    function main()
        run_model(Model{ALCAP}, INPUT, "data/$(TAG).h5")
    end

    function export_files()
        header, profile = read_slice("data/$(TAG).h5", :profile)
        columns = [profile[i] for i in [20, 40, 60, 80]]
        data_export(
            CSV(
                :stratigraphic_column => "data/$(TAG)_sc.csv",
                :age_depth_model      => "data/$(TAG)_adm.csv",
                :metadata => "data/$(TAG)_metadata.toml"),
            header,
            columns)
    end

end

Sinusoid.main()
Sinusoid.export_files()