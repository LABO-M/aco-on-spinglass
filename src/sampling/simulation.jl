module Simulation

using Random, ProgressMeter, Statistics, Distributed, SharedArrays

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
function create_random_symmetric_matrix(N::Int)::Symmetric{Int, Matrix{Float64}}
    A = zeros(Int, N, N)
    A .= rand(0:1, N, N)
    J = Symmetric(A)
    return J
end

# Calculate the total pheromone value
function culculate_TP(N::Int, X::Vector{Int}, h::Float64, J::Float64, C::Float64)
    X_spin = 2*X .- 1
    TP = C
    TP += sum(X_spin .* h) * (1/N)
    for i in 1:N, j in 1:N
        if i != j
            TP += J * X_spin[i] * X_spin[j] * (1/(N*(N-1)))
        end
    end
    return TP
end

# Initialize the simulation up to the initial time t0.
function initialize_simulation(N::Int, X::Vector{Int}, h::Float64, J::Float64, S::Vector{Float64}, Sm::Vector{Float64}, t0::Int, C::Float64)
    for t in 1:t0
        X .= rand(0:1, N)
        TP = culculate_TP(N, X, h, J, C)
        S[t] = (t == 1 ? TP : S[t-1] + TP)
        Sm .+= X .* TP
    end
end

function initialize_simulation(N::Int, X::Vector{Int}, h::Float64, J::Float64, S::Vector{Float64}, Sm::Vector{Float64}, t0::Int, exp_val::Float64, C::Float64)
    for t in 1:t0
        X .= rand(0:1, N)
        TP = culculate_TP(N, X, h, J, C)
        S[t] = (t == 1 ? TP : S[t-1] * exp_val + TP)
        Sm .= (t == 1 ? X .* TP : Sm * exp_val .+ X .* TP)
    end
end

# Main simulation function
function simulate_ants(N::Int, T::Int, t0::Int, alpha::Float64, h::Float64, J::Float64, C::Float64)
    X = zeros(Int, N)
    Sm = zeros(Float64, N)
    S = zeros(Float64, T + t0)

    # Initialization
    initialize_simulation(N, X, h, J, S, Sm, t0, C)

    # Main simulation loop
    for t in (t0 + 1):(t0 + T)
        Zm = Sm ./ S[t-1]
        prob = decision_function.(Zm, alpha)

        X .= rand(Float64, N) .< prob
        TP = culculate_TP(N, X, h, J, C)
        S[t] = S[t-1] + TP
        Sm .+= X .* TP
    end

    # Compute z(t) values for the entire duration
    time_range = 1:(t0 + T)
    Z = S[time_range] ./ discount_factor(collect(time_range))
    return Z
end

function simulate_ants(N::Int, T::Int, t0::Int, alpha::Float64, tau::Int, h::Float64, J::Float64, C::Float64)
    X = zeros(Int, N)
    Sm = zeros(Float64, N)
    S = zeros(Float64, T + t0)
    exp_val = exp(-1 / tau)

    # Initialization
    initialize_simulation(N, X, h, J, S, Sm, t0, exp_val, C)

    # Main simulation loop
    for t in (t0 + 1):(t0 + T)
        Zm = Sm ./ S[t-1]
        prob = decision_function.(Zm, alpha)
        X .= rand(Float64, N) .< prob
        TP = culculate_TP(N, X, h, J, C)
        S[t] = S[t-1] * exp_val + TP
        Sm .= Sm * exp_val .+ X .* TP
    end

    # Compute z(t) values for the entire duration
    time_range = 1:(t0 + T)
    Z = S[time_range] ./ discount_factor(collect(time_range), tau)
    return Z
end

# Function to sample Z values
function sample_ants(N::Int, T::Int, t0::Int, alpha::Float64, tau::Int, samples::Int, h::Float64, J::Float64, C::Float64)::Tuple{Vector{Float64}, Vector{Float64}}
    Z_samples = SharedArray{Float64}((T + t0), samples)

    progressBar = Progress(samples * T, 1, "Samples: ")
    ProgressMeter.update!(progressBar, 0)

    @sync @distributed for i in 1:samples
        if tau == -1
            Z_samples[:, i] = simulate_ants(N, T, t0, alpha, h, J, C)
        else
            Z_samples[:, i] = simulate_ants(N, T, t0, alpha, tau, h, J, C)
        end
    end

    println("Finished simulation")

    # Calculate mean and standard deviation values
    Z_mean = mean(Z_samples, dims=2)
    Z_std = std(Z_samples, dims=2)

    return vec(Z_mean), vec(Z_std)
end

end
