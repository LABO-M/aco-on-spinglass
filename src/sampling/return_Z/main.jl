using ArgParse
include("simulation.jl")
include("output.jl")

function main(args)
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--N"
        help = "Number of questions (quizzes)"
        default = 20
        arg_type = Int

        "--T"
        help = "Total number of time steps"
        default = 1_100_000
        arg_type = Int

        "--ialpha"
        help = "Exponent parameter alpha"
        default = 0.0
        arg_type = Float64

        "--malpha"
        help = "Exponent parameter alpha"
        default = 0.5
        arg_type = Float64

        "--tau"
        help = "Time scale of the pheromone evaporation. Use '-1' for infinite tau."
        default = 100
        arg_type = Int
        
        "--sample"
        help = "Sample size."
        default = 1
        arg_type = Int

        "--h"
        help = "Magnetic field"
        default = 0.01
        arg_type = Float64

        "--J"
        help = "Coupling constant"
        default = 1.0
        arg_type = Float64

        "--chunk_size"
        help = "Chunk size for parallelization"
        default = 10000
        arg_type = Int

    end

    parsed_args = parse_args(args, s)
    N = parsed_args["N"]
    T = parsed_args["T"]
    ialpha = parsed_args["ialpha"]
    malpha = parsed_args["malpha"]
    tau = parsed_args["tau"]
    sample = parsed_args["sample"]
    h = parsed_args["h"]
    J = parsed_args["J"]
    chunk_size = parsed_args["chunk_size"]

    # Log the simulation parameters
    tau_str = (tau == -1) ? "inf" : int_to_SI_prefix(tau)
    println("Running simulation with the following parameters:")
    println("N = $(int_to_SI_prefix(N)), T = $(int_to_SI_prefix(T)), alpha = $(malpha), tau = $(tau_str), sample = $(int_to_SI_prefix(sample)), h = $(h), J = $(J)")

    # Run the simulation
    Z_mean, Z_std = Simulation.sample_ants(N, T, ialpha, malpha, tau, sample, h, J, chunk_size)

    # Output Z values to CSV
    dir_Z = "data/Zt"
    if !isdir(dir_Z)
        mkpath(dir_Z)
    end
    filename_Z = joinpath(dir_Z, "N$(int_to_SI_prefix(N))_T$(int_to_SI_prefix(T))_alpha$(malpha)_tau$(tau_str)_sample$(int_to_SI_prefix(sample))_h$(h)_J$(J).csv")
    save_Z_to_csv(Z_mean, Z_std, filename_Z)
    end

# Entry point of the script
isinteractive() || main(ARGS)
