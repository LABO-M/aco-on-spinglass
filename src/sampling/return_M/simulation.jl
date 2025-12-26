module Simulation

using Random
using Distributed
using SharedArrays
using LinearAlgebra

export decision_function, discount_factor, calculate_pheromone,
       simulate_until_alpha, sample_M

# -----------------------------
# Decision function f(z)  (paper: f(z) = (1-α)/2 + α z)
# -----------------------------
@inline function decision_function(z::Float64, alpha::Float64)::Float64
    return alpha * (z - 0.5) + 0.5
end

# -----------------------------
# Discount factor D(t)
# (not required for the recursion form, but kept if you want it)
# -----------------------------
function discount_factor(t::Vector{Int}, tau::Int)::Vector{Float64}
    return (1.0 .- exp.(-t ./ tau)) ./ (1.0 - exp(-1 / tau))
end
discount_factor(t::Vector{Int})::Vector{Float64} = Float64.(t)

# -----------------------------
# Total pheromone (exp(-Energy))
# Energy style:
#   E = -h * sum_i s_i  - (J/(N-1)) * sum_{i<j} s_i s_j
# where s_i ∈ {-1,+1}
# -----------------------------
function calculate_pheromone(N::Int, X::Vector{Int}, h::Float64, J::Float64)::Float64
    # X ∈ {0,1}^N  ->  s = 2X - 1 ∈ {-1,+1}^N
    s = 2 .* X .- 1

    # field term
    E = -h * sum(s)

    # pair term (use i<j to avoid double counting)
    coef = J / (N - 1)
    for i in 1:(N-1)
        si = s[i]
        @inbounds for j in (i+1):N
            E += -coef * si * s[j]
        end
    end

    return exp(-E)
end

# -----------------------------
# Core simulation:
# Paper-matching sampling rule for Fig.3:
#   - Run annealing with α(t) increasing by alpha_increment
#   - At the FIRST time α reaches target_alpha, sample ONCE:
#       M(i,t) = 2 α(t) (Z(i,t) - 1/2)
#   - Return that M vector immediately (no extra steps)
#
# Notes:
#   - Z is updated via Sm/S recursion with evaporation exp(-1/τ)
#   - We compute M using the SAME α that was used to generate X at that step.
# -----------------------------
function simulate_until_alpha(
    N::Int,
    alpha0::Float64,
    target_alpha::Float64,
    alpha_increment::Float64,
    tau::Int,
    h::Float64,
    J::Float64;
    rng::AbstractRNG = Random.default_rng(),
    max_steps::Int = 10^9
)::Vector{Float64}

    X  = zeros(Int, N)
    Sm = zeros(Float64, N)
    Z  = fill(0.5, N)

    exp_val = exp(-1 / tau)

    # Recursion scalars (store only current S)
    S_prev = 0.0

    alpha = alpha0

    for t in 1:max_steps
        # decision + sampling X (same step α)
        prob = decision_function.(Z, alpha)
        @inbounds for i in 1:N
            X[i] = rand(rng) < prob[i] ? 1 : 0
        end

        # pheromone
        TP = calculate_pheromone(N, X, h, J)

        # update S and Sm (evaporation)
        S_curr = (t == 1) ? TP : (S_prev * exp_val + TP)

        if t == 1
            @inbounds for i in 1:N
                Sm[i] = X[i] * TP
            end
        else
            @inbounds for i in 1:N
                Sm[i] = Sm[i] * exp_val + X[i] * TP
            end
        end

        # update Z
        @inbounds for i in 1:N
            Z[i] = Sm[i] / S_curr
        end

        # compute M using current α (paper definition)
        M = 2.0 * alpha .* (Z .- 0.5)

        # --- IMPORTANT: stop exactly when α reaches target ---
        if alpha >= target_alpha
            return M
        end

        # increment α for next step
        alpha = min(alpha + alpha_increment, target_alpha)

        # shift S
        S_prev = S_curr
    end

    error("Reached max_steps without hitting target_alpha. Check alpha_increment / max_steps.")
end

# -----------------------------
# Sampling wrapper (Fig.3 style):
#   Repeat "simulate_until_alpha" for 'samples' trials
#   and concatenate M values into a single vector (length N*samples)
#
# For Distributed:
#   - Make sure this module is available on all workers (e.g. @everywhere include("..."))
# -----------------------------
function sample_ants(
    N::Int,
    alpha0::Float64,
    target_alpha::Float64,
    alpha_increment::Float64,
    tau::Int,
    samples::Int,
    h::Float64,
    J::Float64,
    seed::Int;
    use_distributed::Bool = true
)::Vector{Float64}

    M_samples = SharedArray{Float64}(N, samples)

    if use_distributed && nworkers() > 1
        @sync @distributed for i in 1:samples
            rng = MersenneTwister(seed + i - 1)
            M_samples[:, i] = simulate_until_alpha(
                N, alpha0, target_alpha, alpha_increment, tau, h, J;
                rng=rng
            )
        end
    else
        for i in 1:samples
            rng = MersenneTwister(seed + i - 1)
            M_samples[:, i] = simulate_until_alpha(
                N, alpha0, target_alpha, alpha_increment, tau, h, J;
                rng=rng
            )
        end
    end

    # concatenate to 1D vector like the paper histogram input
    return vcat(M_samples...)
end

end # module