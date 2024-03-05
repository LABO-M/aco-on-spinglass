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
function culculate_max_energy(J::Float64, h::Float64)
    energy = -h + -J
    return energy
end

# Calculate the total pheromone value
function culculate_energy(N::Int, X::Vector{Int}, h::Float64, J::Float64)
    X_spin = 2*X .- 1
    energy = 0
    energy += sum(X_spin .* -h) * (1/N)
    for i in 1:N, j in 1:N
        if i != j
            energy += -J * X_spin[i] * X_spin[j] * (1/(N*(N-1)))
        end
    end

    return energy
end

# Calculate the total pheromone value
function culculate_Pheromone(energy::Float64)
    Pheromone = exp(-energy)
    return Pheromone
end

## Initialize the simulation up to the initial time t0.
#function initialize_simulation(N::Int, X::Vector{Int}, h::Float64, J::Float64, S::Vector{Float64}, Sm::Vector{Float64}, t0::Int, C::Float64)
#    for t in 1:t0
#        X .= rand(0:1, N)
#        TP = culculate_TP(N, X, h, J, C)
#        S[t] = (t == 1 ? TP : S[t-1] + TP)
#        Sm .+= X .* TP
#    end
#end

#function initialize_simulation(N::Int, X::Vector{Int}, h::Float64, J::Float64, S::Vector{Float64}, Sm::Vector{Float64}, t0::Int, exp_val::Float64, C::Float64)
#    for t in 1:t0
#        X .= rand(0:1, N)
#        TP = culculate_TP(N, X, h, J, C)
#        S[t] = (t == 1 ? TP : S[t-1] * exp_val + TP)
#        Sm .= (t == 1 ? X .* TP : Sm * exp_val .+ X .* TP)
#    end
#end

# Main simulation function
function simulate_ants(N::Int, T::Int, ialpha::Float64, malpha::Float64, h::Float64,  J::Float64)
    X = zeros(Int, N)
    Sm = zeros(Float64, N)
    S = zeros(Float64, T)
    Zm = ones(Float64, N) * 0.5
    ACO_energy = []

    #Find the max energy and max spin
    #Max_energy = culculate_max_energy(J, h)

    # Main simulation loop
    for t in 1:T
        prob = decision_function.(Zm, ialpha)
        X .= rand(Float64, N) .< prob
        energy = culculate_energy(N, X, h, J)
        Pheromone = culculate_Pheromone(energy)
        push!(ACO_energy, energy)
        S[t] = (t == 1 ? Pheromone : S[t-1] + Pheromone)
        Sm .+= X .* Pheromone
        Zm = Sm ./ S[t]
        if t == 10000
            ialpha = malpha
        end
    end

    return ACO_energy
end

function simulate_ants(N::Int, T::Int, ialpha::Float64, malpha::Float64, tau::Int, h::Float64, J::Float64)
    X = zeros(Int, N)
    Sm = zeros(Float64, N)
    S = zeros(Float64, T)
    Zm = ones(Float64, N) * 0.5
    exp_val = exp(-1 / tau)
    ACO_energy = []

    #Find the max energy and max spin
    #Max_energy = culculate_max_energy(J, h)

    # Main simulation loop
    for t in 1:T
        prob = decision_function.(Zm, ialpha)
        X .= rand(Float64, N) .< prob
        energy = culculate_energy(N, X, h, J)
        Pheromone = culculate_Pheromone(energy)
        push!(ACO_energy, energy)
        S[t] = (t == 1 ? Pheromone : S[t-1] * exp_val + Pheromone)
        Sm .= (t == 1 ? X .* Pheromone : Sm * exp_val .+ X .* Pheromone)
        Zm = Sm ./ S[t]
        if t == 10000
            ialpha = malpha
        end
    end

    return ACO_energy
end


# Function to sample Z values
function sample_ants(N::Int, T::Int, ialpha::Float64, malpha::Float64, tau::Int, samples::Int, h::Float64, J::Float64, chunk_size::Int)::Tuple{Vector{Float64}, Vector{Float64}}
    Z_samples = SharedArray{Float64}(T, samples)

    progressBar = Progress(samples * T, 1, "Samples: ")
    ProgressMeter.update!(progressBar, 0)

    @sync @distributed for i in 1:samples
        if tau == -1
            Z_samples[:, i] = simulate_ants(N, T, ialpha, malpha, h, J)
        else
            Z_samples[:, i] = simulate_ants(N, T, ialpha, malpha, tau, h, J)
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