const FACIES = [
    ALCAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate=500u"m/Myr",
        extinction_coefficient=0.6u"m^-1",
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

const PERIOD1 = 2.5u"Myr"
const PERIOD2 = 0.1u"Myr"
const AMPLITUDE1 = 20.0u"m"
const AMPLITUDE2 = 2.0u"m"

const MAIN_SL_SINUSOID = t -> AMPLITUDE1 * cos(2π * t / PERIOD1) + AMPLITUDE2 * sin(2π * t / PERIOD2)
const INIT_TOPO = (x, y) -> -x / 300.0

const PHYS_SCALE = 150u"m"
const GRID_SIZE_Y = 30
const GRID_SIZE_X = 230

const DELTA_T = 0.0001u"Myr"
const STEPS_MAIN = 50000

const SUBSIDENCE_RATE = 20u"m/Myr"
const DISINTGRATION_RATE = 50.0u"m/Myr"
const INSOLATION = 400.0u"W/m^2"
const SEDIMENT_BUFFER_SIZE = 100
const DEPOSITIONAL_RESOLUTION = 0.5u"m"
const CEMENTATION_TIME = 1u"kyr"
const CA_INTERVAL = 1

const PROFILE_WRITE_INTERVAL = 10
const TOPOGRAPHY_WRITE_INTERVAL = 100