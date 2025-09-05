using ArgParse
include("mean_sampling.jl")
using .simulation

function main(args)
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--power"
        help = "the number that power of 2"
        default = 5
        arg_type = Int

        "--tau"
        help = "evaporation hyper parameter"
        default = 100.0
        arg_type = Float64

        "--beta"
        help = "coefficient of external field"
        default = 2
        arg_type = Float64

        "--start_seed"
        help = "random seed at the start point"
        default = 42
        arg_type = Int

        "--magnetic"
        help = "magnetic field value"
        default = 0.0
        arg_type = Float64

        "--interaction"
        help = "coefficient of interaction"
        default = 0.1
        arg_type = Float64

    end

    parsed_args = parse_args(args, s)
    power = parsed_args["power"]
    tau = parsed_args["tau"]
    beta = parsed_args["beta"]
    start_seed = parsed_args["start_seed"]
    magnetic = parsed_args["magnetic"]
    interaction = parsed_args["interaction"]

    file = simulation.sampling(power, tau, beta, start_seed, magnetic, interaction)
    println("Saved $file")

    
end

isinteractive() || main(ARGS)
