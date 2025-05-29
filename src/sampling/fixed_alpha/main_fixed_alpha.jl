using ArgParse
include("fixed_alpha_mean.jl")

function main(args)
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--power"
        help = "the number that power of 2"
        default = 5
        arg_type = Int

        "--beta"
        help = "coefficient of external field"
        default = 0.01
        arg_type = Float64

        "--start_seed"
        help = "random seed at the start point"
        default = 42
        arg_type = Int

        "--sample"
        help = "sample size"
        default = 10000
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
    beta = parsed_args["beta"]
    start_seed = parsed_args["start_seed"]
    sample = parsed_args["sample"]
    magnetic = parsed_args["magnetic"]
    interaction = parsed_args["interaction"]

    filename_spins_mean = fixed_alpha_mean.sampling(power, beta, start_seed, sample, magnetic, interaction)
    println("Saved $filename_spins_mean")
end

isinteractive() || main(ARGS)