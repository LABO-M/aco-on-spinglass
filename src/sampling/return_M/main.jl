using ArgParse
include("simulation.jl")
include("output.jl")

function main(args)
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--seed"
        help = "Random seed value"
        default = 1625
        arg_type = Int

        "--N"
        help = "Number of questions (quizzes)"
        default = 100
        arg_type = Int

        "--alpha"
        help = "Exponent parameter alpha"
        default = 0.0
        arg_type = Float64

        "--end_alpha"
        help = "Final value of alpha"
        default = 0.99
        arg_type = Float64

        "--alpha_increment"
        help = "Increment of alpha"
        default = 1.0e-6
        arg_type = Float64

        "--tau"
        help = "Time scale of the pheromone evaporation. Use '-1' for infinite tau."
        default = 1000
        arg_type = Int
        
        "--sample"
        help = "Sample size."
        default = 1000
        arg_type = Int

        "--h"
        help = "Magnetic field"
        default = 0.001
        arg_type = Float64

        "--J"
        help = "Coupling constant"
        default = 0.1
        arg_type = Float64

    end

    parsed_args = parse_args(args, s)
    seed = parsed_args["seed"]
    N = parsed_args["N"]
    alpha = parsed_args["alpha"]
    end_alpha = parsed_args["end_alpha"]
    alpha_increment = parsed_args["alpha_increment"]
    tau = parsed_args["tau"]
    sample = parsed_args["sample"]
    h = parsed_args["h"]
    J = parsed_args["J"]

    # Log the simulation parameters
    tau_str = (tau == -1) ? "inf" : int_to_SI_prefix(tau)
    println("Running simulation with the following parameters:")
    println("seed = $(seed), N = $(int_to_SI_prefix(N)), alpha = $(end_alpha), tau = $(tau_str), sample = $(int_to_SI_prefix(sample)), h = $(h), J = $(J)")

    # Run the simulation
    Z_M = Simulation.sample_ants(N, alpha, end_alpha, alpha_increment, tau, sample, h, J, seed)

    # Output Z values to CSV

    alpha_increment_str = string(alpha_increment)

    dir_M = "/home/mori-lab/shimizu/ACO/ACO-on-SpinGlass/src/sampling/return_M/data/seed$(seed)/tau$(tau_str)_h$(h)_J$(J)_α_inc$(alpha_increment_str)"
    if !isdir(dir_M)
        mkpath(dir_M)
    end
    filename_M = joinpath(dir_M, "N$(int_to_SI_prefix(N))_alpha$(end_alpha)_sample$(int_to_SI_prefix(sample)).csv")
    save_Z_to_csv(Z_M, filename_M)
    end

# Entry point of the script
isinteractive() || main(ARGS)