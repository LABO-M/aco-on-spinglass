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
    for t in B
        energy = C
        for i in 1:N
            energy += h * (2 * parse(Int, t[i]) -1) * (1/N)
        J[i, j] * (2 * parse(Int, t[i]) -1) * (2 * parse(Int, t[j]) -1) * (1/(N*(N-1)))
            end
        end
        if energy > max_energy
            max_energy = energy
        end
    end
    return max_energy
end

# Calculate the total pheromone value
function culculate_TP(N::Int, X::Vector{Int}, h::Float64, J::Symmetric{Int, Matrix{Int}}, C::Float64)
    X_spin = 2*X .- 1
    TP = C
    TP += sum(X_spin .* h) * (1/N)
    for i in 1:N, j in 1:N
        if i != j
            TP += J[i, j] * X_spin[i] * X_spin[j] * (1/(N*(N-1)))
        end
    end
    return TP
end

# Main simulation function
function simulate_ants(N::Int, T::Int, alpha::Float64, h::Float64, C::Float64)
    X = zeros(Int, N)
    Sm = zeros(Float64, N)
    S = zeros(Float64, T)
    Zm = ones(Float64, N) ./ 2
    J = create_random_symmetric_matrix(N)
    ACO_energy = []

    #Find the max energy and max spin
    Max_energy = culculate_max_energy(N, J, h, C)

    # Main simulation loop
    for t in 1:T
        prob = decision_function.(Zm, alpha)
        X .= rand(Float64, N) .< prob
        TP = culculate_TP(N, X, h, J, C)
        push!(ACO_energy, TP/Max_energy)
        S[t] = (t == 1 ? TP : S[t-1] + TP)
        Sm .+= X .* TP
        Zm = Sm ./ S[t]
        if alpha < 1.0
            if t % 10000 == 0
                alpha += 0.01
            end
        end

    end

    return ACO_energy
end

function simulate_ants(N::Int, T::Int, alpha::Float64, tau::Int, h::Float64, C::Float64)
    X = zeros(Int, N)
    Sm = zeros(Float64, N)
    S = zeros(Float64, T)
    Zm = ones(Float64, N) ./ 2
    exp_val = exp(-1 / tau)
    J = create_random_symmetric_matrix(N)
    ACO_energy = []

    #Find the max energy and max spin
    Max_energy = culculate_max_energy(N, J, h, C)

    # Main simulation loop
    for t in 1:T
        prob = decision_function.(Zm, alpha)
        X .= rand(Float64, N) .< prob
        TP = culculate_TP(N, X, h, J, C)
        push!(ACO_energy, TP/Max_energy)
        S[t] = (t == 1 ? TP : S[t-1] * exp_val + TP)
        Sm .= Sm * exp_val .+ X .* TP
        Zm = Sm ./ S[t]
        if alpha < 1.0
            if t % 10000 == 0
                alpha += 0.01
            end
        end
    end

    return ACO_energy
end


# Function to sample Z values
function sample_ants(N::Int, T::Int, alpha::Float64, tau::Int, samples::Int, h::Float64, C::Float64, chunk_size::Int)::Tuple{Vector{Float64}, Vector{Float64}}
    Z_samples = SharedArray{Float64}((T), samples)

    progressBar = Progress(samples * T, 1, "Samples: ")
    ProgressMeter.update!(progressBar, 0)

    @sync @distributed for i in 1:samples
        if tau == -1
            Z_samples[:, i] = simulate_ants(N, T, alpha, h, C)
        else
            Z_samples[:, i] = simulate_ants(N, T, alpha, tau, h, C)
        end
    end

    println("Finished simulation")

    # Calculate mean and standard deviation values
    Z_sample_mean = mean(Z_samples, dims=2)
    Z_mean = []
    Z_std = []
    num_chunk = Int(T / chunk_size)
    for i in 1:num_chunk
        push!(Z_mean, mean(Z_sample_mean[(i-1)*chunk_size+1:i*chunk_size]))
        push!(Z_std, std(Z_sample_mean[(i-1)*chunk_size+1:i*chunk_size]))
    end

    return vec(Z_mean), vec(Z_std)
end
