module Simulation

using Random, ProgressMeter, Statistics, Distributed, SharedArrays, LinearAlgebra

#const DEFAULT_EPSILON = 0.01  # Define epsilon as a constant at the top of the code

# Decision function f(z)
function decision_function(z::Float64, alpha::Float64)
    return alpha * (z - 0.5) + 0.5
end

# Discount factor D(t) for finite tau
function discount_factor(t::Vector{Int}, tau::Int)::Vector{Float64}
    return (1.0 .- exp.(-t ./ tau)) ./ (1.0 - exp(-1/tau))
end

# Discount factor D(t) for infinite tau
function discount_factor(t::Vector{Int})::Vector{Float64}
    return t
end

#Create a random spin symmetric matrix
function create_random_symmetric_matrix(N::Int)::Symmetric{Int, Matrix{Int}}
    A = zeros(Int, N, N)
    A .= rand(0:1, N, N)
    B = 2*A .- 1
    J = Symmetric(B)
    return J
end

function all_combinations(N::Int)
    result = []
    for i in 0:2^N-1
        push!(result, string(i, base=2, pad=N))
    end
    return result
end

function culculate_max_energy(N::Int, J::Symmetric{Int, Matrix{Int}}, h::Float64, C::Float64)
    B = all_combinations(N)
    max_energy = -10
    max_spin = []
    for t in B
        energy = C
        for i in 1:N
            energy += h * (2 * parse(Int, t[i]) -1)
        end
        for i in 1:N, j in 1:N
            if i != j
                energy += J[i, j] * (2 * parse(Int, t[i]) -1) * (2 * parse(Int, t[j]) -1)
            end
        end
        if energy == max_energy
            max_energy = energy
            push!(max_spin, t)
        elseif energy > max_energy
            max_energy = energy
            max_spin = [t]
        end
    end
    return max_energy, max_spin
end

# Calculate the total pheromone value
function culculate_TP(N::Int, X::Vector{Int}, h::Float64, J::Symmetric{Int, Matrix{Int}}, C::Float64)
    X_spin = 2*X .- 1
    TP = C
    TP += sum(X_spin .* h)
    for i in 1:N, j in 1:N
        if i != j
            TP += J[i, j] * X_spin[i] * X_spin[j]
        end
    end
    return TP
end

# Initialize the simulation up to the initial time t0.
function initialize_simulation(N::Int, X::Vector{Int}, h::Float64, J::Symmetric{Int, Matrix{Int}}, S::Vector{Float64}, Sm::Vector{Float64}, t0::Int, C::Float64)
    for t in 1:t0
        X .= rand(0:1, N)
        TP = culculate_TP(N, X, h, J, C)
        S[t] = (t == 1 ? TP : S[t-1] + TP)
        Sm .+= X .* TP
    end
end

function initialize_simulation(N::Int, X::Vector{Int}, h::Float64, J::Symmetric{Int, Matrix{Int}}, S::Vector{Float64}, Sm::Vector{Float64}, t0::Int, exp_val::Float64, C::Float64)
    for t in 1:t0
        X .= rand(0:1, N)
        TP = culculate_TP(N, X, h, J, C)
        S[t] = (t == 1 ? TP : S[t-1] * exp_val + TP)
        Sm .= (t == 1 ? X .* TP : Sm * exp_val .+ X .* TP)
    end
end

# Main simulation function
function simulate_ants(N::Int, T::Int, t0::Int, alpha::Float64, h::Float64, C::Float64)
    X = zeros(Int, N)
    Sm = zeros(Float64, N)
    S = zeros(Float64, T + t0)
    J = create_random_symmetric_matrix(N)
    max_energy = -10
    max_spin = []

    #Find the max energy and max spin
    Max_energy, Max_spin = culculate_max_energy(N, J, h, C)

    # Initialization
    initialize_simulation(N, X, h, J, S, Sm, t0, C)

    # Main simulation loop
    for t in (t0 + 1):(t0 + T)
        Zm = Sm ./ S[t-1]
        prob = decision_function.(Zm, alpha)
        X .= rand(Float64, N) .< prob
        TP = culculate_TP(N, X, h, J, C)
        if TP > max_energy
            max_energy = TP
            max_spin = copy(X)
        end
        S[t] = S[t-1] + TP
        Sm .+= X .* TP
    end

    return J, max_energy, max_spin, Max_energy, Max_spin
end

function simulate_ants(N::Int, T::Int, t0::Int, alpha::Float64, tau::Int, h::Float64, C::Float64)
    X = zeros(Int, N)
    Sm = zeros(Float64, N)
    S = zeros(Float64, T + t0)
    exp_val = exp(-1 / tau)
    J = create_random_symmetric_matrix(N)
    max_energy = -10
    max_spin = []

    #Find the max energy and max spin
    Max_energy, Max_spin = culculate_max_energy(N, J, h, C)

    # Initialization
    initialize_simulation(N, X, h, J, S, Sm, t0, exp_val, C)

    # Main simulation loop
    for t in (t0 + 1):(t0 + T)
        Zm = Sm ./ S[t-1]
        prob = decision_function.(Zm, alpha)
        X .= rand(Float64, N) .< prob
        TP = culculate_TP(N, X, h, J, C)
        if TP > max_energy
            max_energy = TP
            max_spin = copy(X)
        end
        S[t] = S[t-1] * exp_val + TP
        Sm .= Sm * exp_val .+ X .* TP
    end

    return J, max_energy, max_spin, Max_energy, Max_spin
end

# Function to sample Z values
function sample_ants(N::Int, T::Int, t0::Int, alpha::Float64, tau::Int, h::Float64, C::Float64)
    if tau == -1
        J, max_energy, max_spin, Max_energy, Max_spin = simulate_ants(N, T, t0, alpha, h, C)
    else
        J, max_energy, max_spin, Max_energy, Max_spin = simulate_ants(N, T, t0, alpha, tau, h, C)
    end

    println("Finished simulation")
    return J, max_energy, max_spin, Max_energy, Max_spin
end

end