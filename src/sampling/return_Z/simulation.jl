module Simulation

using Random, ProgressMeter, Statistics, Distributed, SharedArrays, LinearAlgebra

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

#Calculate the maximum energy
function culculate_max_energy(J::Float64, h::Float64, C::Float64)
    TP = h + J + C
    return TP
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
function simulate_ants(N::Int, T::Int, t0::Int, alpha::Float64, h::Float64,  J::Float64, C::Float64)
    X = zeros(Int, N)
    Sm = zeros(Float64, N)
    S = zeros(Float64, T + t0)
    ACO_energy = []

    #Find the max energy and max spin
    Max_energy = culculate_max_energy(J, h, C)

    # Initialization
    initialize_simulation(N, X, h, J, S, Sm, t0, C)

    # Main simulation loop
    for t in (t0 + 1):(t0 + T)
        Zm = Sm ./ S[t-1]
        prob = decision_function.(Zm, alpha)
        X .= rand(Float64, N) .< prob
        TP = culculate_TP(N, X, h, J, C)
        push!(ACO_energy, TP/Max_energy)
        S[t] = S[t-1] + TP
        Sm .+= X .* TP
        if alpha < 1.0
            if (t-t0) % 10000 == 0
                alpha += 0.01
            end
        end

    end

    return ACO_energy
end

function simulate_ants(N::Int, T::Int, t0::Int, alpha::Float64, tau::Int, h::Float64, J::Float64, C::Float64)
    X = zeros(Int, N)
    Sm = zeros(Float64, N)
    S = zeros(Float64, T + t0)
    exp_val = exp(-1 / tau)
    ACO_energy = []

    #Find the max energy and max spin
    Max_energy = culculate_max_energy(J, h, C)

    # Initialization
    initialize_simulation(N, X, h, J, S, Sm, t0, exp_val, C)

    # Main simulation loop
    for t in (t0 + 1):(t0 + T)
        Zm = Sm ./ S[t-1]
        prob = decision_function.(Zm, alpha)
        X .= rand(Float64, N) .< prob
        TP = culculate_TP(N, X, h, J, C)
        push!(ACO_energy, TP/Max_energy)
        S[t] = S[t-1] * exp_val + TP
        Sm .= Sm * exp_val .+ X .* TP
        if alpha < 1.0
            if (t-t0) % 10000 == 0
                alpha += 0.01
            end
        end
    end

    return ACO_energy
end


# Function to sample Z values
function sample_ants(N::Int, T::Int, t0::Int, alpha::Float64, tau::Int, samples::Int, h::Float64, J::Float64, C::Float64, chunk_size::Int)::Tuple{Vector{Float64}, Vector{Float64}}
    Z_samples = SharedArray{Float64}(T, samples)

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
    Z_sample_mean = mean(Z_samples, dims=2)
    Z_mean = []
    Z_std = []
    num_chunk = Int(T/chunk_size)
    for i in 1:num_chunk
        push!(Z_mean, mean(Z_sample_mean[(i-1)*chunk_size+1:i*chunk_size]))
        push!(Z_std, std(Z_sample_mean[(i-1)*chunk_size+1:i*chunk_size]))
    end

    return vec(Z_mean), vec(Z_std)
end

end