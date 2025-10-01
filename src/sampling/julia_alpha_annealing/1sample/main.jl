using ArgParse
include("long_simulation.jl")
using .simulation

function main(args)
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--power"
        help = "the number that power of 2"
        default = 8
        arg_type = Int

        "--tau"
        help = "evaporation hyper parameter"
        default = 100.0
        arg_type = Float64

        "--beta"
        help = "coefficient of external field"
        default = 5.0
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

        "--alpha_step"
        help = "the increase of alpha"
        default = 0.001
        arg_type = Float64

        "--iters_per_alpha"
        help = "iteration per alpha"
        default = 100000
        arg_type = Int

        "--alpha_max"
        help = "max alpha value"
        default = 1.0
        arg_type = Float64

    end

    parsed_args = parse_args(args, s)
    power = parsed_args["power"]
    tau = parsed_args["tau"]
    beta = parsed_args["beta"]
    start_seed = parsed_args["start_seed"]
    magnetic = parsed_args["magnetic"]
    interaction = parsed_args["interaction"]
    alpha_step = parsed_args["alpha_step"]
    iters_per_alpha = parsed_args["iters_per_alpha"]
    alpha_max = parsed_args["alpha_max"]

    file = simulation.sampling(
        power, tau, beta, start_seed, magnetic, interaction;
        alpha_step=alpha_step,
        iters_per_alpha=iters_per_alpha,
        alpha_max=alpha_max
    )
    println("Saved $file")

    
end

isinteractive() || main(ARGS)
