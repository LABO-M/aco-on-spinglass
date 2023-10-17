using ArgParse
include("simulation.jl")
include("output.jl")

function main(args)
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--N"
        help = "Number of questions (quizzes)"
        default = 100
        arg_type = Int
        "--T"
        help = "Total number of time steps"
        default = 1_000_000
        arg_type = Int
        "--t0"
        help = "Initial time steps for initialization"
        default = 100_000
        arg_type = Int
        "--alpha"
        help = "Exponent parameter alpha"
        default = 1.0
        arg_type = Float64
        "--tau"
        help = "Time scale of the pheromone evaporation. Set to 'inf' for infinite tau."
        arg_type = Union{Int, String}
        default = 10_000
    end

    parsed_args = parse_args(args, s)
    N = parsed_args["N"]
    T = parsed_args["T"]
    t0 = parsed_args["t0"]
    alpha = parsed_args["alpha"]
    tau = parsed_args["tau"]

    println(typeof(N))
    println(typeof(T))
    println(typeof(t0))
    println(typeof(alpha))
    println(typeof(tau))


    # Log the simulation parameters
    println("Running simulation with the following parameters:")
    println("N = $(format_num(N)), T = $(format_num(T)), t0 = $(format_num(t0)), alpha = $(alpha), tau = $(tau)")

    # Run the simulation
    if tau == "inf"
        Z = Simulation.simulate_ants(N, T, t0, alpha)
    else
        Z = Simulation.simulate_ants(N, T, t0, alpha, tau)
    end

    # Output Z values to CSV
    dir_Z = "data/Zt"
    if !isdir(dir_Z)
        mkpath(dir_Z)
    end
    tau_str = (typeof(tau) == Int) ? format_num(tau) : string(tau)
    filename_Z = joinpath(dir_Z, "N$(format_num(N))_T$(format_num(T))_t0$(format_num(t0))_alpha$(alpha)_tau$(tau_str).csv")
    save_Z_to_csv(Z, filename_Z)
end

# Entry point of the script
isinteractive() || main(ARGS)
