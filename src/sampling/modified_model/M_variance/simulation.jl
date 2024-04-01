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
function culculate_min_energy(J::Float64, h::Float64)
    TP = -h + -J
    return TP
end

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
function simulate_ants(N::Int, s::Int, alpha::Float64, falpha::Float64, tau::Int, h::Float64, J::Float64, progressBar::ProgressMeter.Progress)
    T = s * 101
    X = zeros(Int, N)
    Sm = zeros(Float64, N)
    S = zeros(Float64, T)
    Zm = ones(Float64, N) * 0.5
    exp_val = exp(-1 / tau)
    M = zeros(Float64, N)
    M_temp_matrix = zeros(Float64, N, s)
    M_variance_series = zeros(Float64, 101)


    # Main simulation loop
    for t in 1:T
        prob = decision_function.(Zm, alpha)
        X .= rand(Float64, N) .< prob
        TP = culculate_Pheromone(N, X, h, J)
        S[t] = (t == 1 ? TP : S[t-1] * exp_val + TP)
        Sm .= (t == 1 ? X .* TP : Sm * exp_val .+ X .* TP)
        Zm = Sm ./ S[t]
        M = 2 * alpha * (Zm .- 0.5)
        M_temp_matrix[:, mod1(t, s)] = M
        if alpha < falpha + 0.001
            if t % s == 0
                M_time_mean = mean(M_temp_matrix, dims=2)
                M_variance_series[convert(Int, t/s)] = var(M_time_mean)
                alpha += 0.01
            end
        end
        next!(progressBar)
    end

    return M_variance_series
end

#function simulate_ants(N::Int, T::Int, alpha::Float64, calpha::Float64, falpha::Float64, tau::Int, h::Float64,  J::Float64, progressBar::ProgressMeter.Progress)
#    T = s * 101
#    X = zeros(Int, N)
#    Sm = zeros(Float64, N)
#    S = zeros(Float64, T)
#    Zm = ones(Float64, N) * 0.5
#    exp_val = exp(-1 / tau)
#    M = zeros(Float64, N)
#    M_all_series = zeros(Float64, N, T)
#    M_series_mean = zeros(Float64, N * 101)
#
#
#    # Main simulation loop
#    for t in 1:T
#        prob = decision_function.(Zm, alpha)
#        X .= rand(Float64, N) .< prob
#        TP = culculate_Pheromone(N, X, h, J)
#        S[t] = (t == 1 ? TP : S[t-1] * exp_val + TP)
#        Sm .= (t == 1 ? X .* TP : Sm * exp_val .+ X .* TP)
#        Zm = Sm ./ S[t]
#        if alpha < falpha + 0.001
#            if t % s == 0
#                alpha += 0.01
#            end
#        end
#        M = 2 * alpha * (Zm .- 0.5)
#        M_all_series[:, t] = M
#        next!(progressBar)
#    end
#
#    for i in 1:101
#        for j in 1:s
#
#
#    return M_series
#end


# Function to sample Z values
function sample_ants(N::Int, s::Int, alpha::Float64, falpha::Float64, tau::Int, h::Float64, J::Float64)

    T = s * 101
    M_samples = zeros(Float64, 101)
    progressBar = Progress(T, 1, "Samples: ")
    ProgressMeter.update!(progressBar, 0)

    M_samples = simulate_ants(N, s, alpha, falpha, tau, h, J, progressBar)

    println("Finished simulation")

    return M_samples

end

end
