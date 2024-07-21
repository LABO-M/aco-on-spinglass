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

##Calculate the maximum energy
#function culculate_min_energy(J::Float64, h::Float64)
#    TP = -h + -J
#    return TP
#end

# Calculate the total pheromone value
function culculate_Pheromone(N::Int, X::Vector{Int}, h::Float64, J::Float64)
    X_spin = 2*X .- 1
    TP = 0
    TP += sum(X_spin .* -h)
    for i in 1:N, j in 1:N
        if i != j
            TP += -J * X_spin[i] * X_spin[j] * (1/(N-1))
        end
    end
    Pheromone = exp(-TP)

    return Pheromone
end

# Main simulation function
# Remain for loop to modify easily
function simulate_ants(N::Int, T::Int, alpha::Float64, end_alpha::Float64, alpha_increment::Float64, tau::Int, h::Float64, J::Float64, progressBar::ProgressMeter.Progress)
    X = zeros(Int, N)
    Sm = zeros(Float64, N)
    S = zeros(Float64, T)
    Zm = ones(Float64, N) * 0.5
    exp_val = exp(-1 / tau)
    M = zeros(Float64, N)


    # Main simulation loop
    for t in 1:T
        prob = decision_function.(Zm, alpha)
        X .= rand(Float64, N) .< prob
        TP = culculate_Pheromone(N, X, h, J)
        S[t] = (t == 1 ? TP : S[t-1] * exp_val + TP)
        Sm .= (t == 1 ? X .* TP : Sm * exp_val .+ X .* TP)
        Zm = Sm ./ S[t]
        if alpha < end_alpha
            alpha += alpha_increment
        end
        next!(progressBar)
    end
    M = 2 * alpha * (Zm .- 0.5)

    return M
end

# Function to sample Z values
function sample_ants(N::Int, alpha::Float64, end_alpha::Float64, alpha_increment::Float64, tau::Int, samples::Int, h::Float64, J::Float64)
    Z_samples = SharedArray{Float64}(N, samples)

    T = convert(Int, round((end_alpha + 0.1) / alpha_increment))

    progressBar = Progress(samples * T, 1, "Samples: ")
    ProgressMeter.update!(progressBar, 0)

    @sync @distributed for i in 1:samples
        Z_samples[:, i] = simulate_ants(N, T, alpha, end_alpha, alpha_increment, tau, h, J, progressBar)
        next!(progressBar)
    end

    println("Finished simulation")

    #Create M vector
    Z_M = vcat(Z_samples...)

    return Z_M

end

end
