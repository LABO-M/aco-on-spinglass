using ArgParse
include("tesannealing.jl")
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
        default = 0.01
        arg_type = Float64

        "--start_alpha"
        help = "alpha value at the start point"
        default = 0.0
        arg_type = Float64

        "--start_seed"
        help = "random seed at the start point"
        default = 42
        arg_type = Int

        "--iter"
        help = "iteration number"
        default = 100000
        arg_type = Int

        "--sample"
        help = "sample size"
        default = 100
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
    start_alpha = parsed_args["start_alpha"]
    start_seed = parsed_args["start_seed"]
    iter = parsed_args["iter"]
    sample = parsed_args["sample"]
    magnetic = parsed_args["magnetic"]
    interaction = parsed_args["interaction"]

    filename_z_mean, filename_spins_mean = simulation.sampling(power, tau, beta, start_alpha, start_seed, iter, sample, magnetic, interaction)
    println("Saved $filename_z_mean and $filename_spins_mean")

    
end

isinteractive() || main(ARGS)
