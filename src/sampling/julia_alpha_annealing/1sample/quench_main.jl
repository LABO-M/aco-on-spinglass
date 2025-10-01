using ArgParse
include("quench_simulation.jl")
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

        "--alpha_fixed"
        help = "the increase of alpha"
        default = 0.3
        arg_type = Float64

        "--iters_burnin"
        help = "iteration per alpha"
        default = 100000
        arg_type = Int

        "--iters_fixed"
        help = "max alpha value"
        default = 100000
        arg_type = Int

    end

    parsed_args = parse_args(args, s)
    power = parsed_args["power"]
    tau = parsed_args["tau"]
    beta = parsed_args["beta"]
    start_seed = parsed_args["start_seed"]
    magnetic = parsed_args["magnetic"]
    interaction = parsed_args["interaction"]
    alpha_fixed = parsed_args["alpha_fixed"]
    iters_burnin = parsed_args["iters_burnin"]
    iters_fixed = parsed_args["iters_fixed"]

    file = simulation.simulate_two_phase(
        power, tau, beta, start_seed, magnetic, interaction;
        alpha_fixed=alpha_fixed,
        iters_burnin=iters_burnin,
        iters_fixed=iters_fixed
    )
    println("Saved $file")

    
end

isinteractive() || main(ARGS)
