module Simulation

using Random
using Distributed
using SharedArrays

export decision_function, calculate_pheromone,
       simulate_to_alpha_and_sample_once, sample_ants

# ------------------------------------------------------------
# Linear decision function (paper: f(z) = (1-α)/2 + α z)
# z ∈ [0,1], α ∈ [0,1]
# ------------------------------------------------------------
@inline function decision_function(z::Float64, alpha::Float64)::Float64
    return (1.0 - alpha) * 0.5 + alpha * z
end

# ------------------------------------------------------------
# Total pheromone: exp(-E)
# Paper energy (uniform h and J):
#   E = -h * Σ s_i  - (J/(N-1)) * Σ_{i≠j} s_i s_j
# where s_i = 2X_i - 1 ∈ {-1,+1}
#
# For uniform coupling:
#   Σ_{i≠j} s_i s_j = (Σ s_i)^2 - Σ s_i^2 = (Σ s_i)^2 - N
#
# So:
#   E = -h * m  - (J/(N-1)) * (m^2 - N)
# with m = Σ s_i = 2*sum(X) - N
# ------------------------------------------------------------
@inline function calculate_pheromone(N::Int, X::Vector{Int}, h::Float64, J::Float64)::Float64
    m = 2 * sum(X) - N
    E = -h * m - (J / (N - 1)) * (m*m - N)
    return exp(-E)
end

# ------------------------------------------------------------
# Core simulation for Fig.3-style sampling:
#   - Slow annealing: alpha starts from alpha0 and increases by alpha_increment
#   - At the first step where the CURRENT alpha == end_alpha (used for sampling),
#     return a single sample of M:
#        M(i) = 2α ( Z(i) - 1/2 )
#
# Notes:
#  - We do NOT add extra burn-in after reaching end_alpha (paper samples at the target α).
#  - tau = -1 means infinite tau => no evaporation (evap = 1).
# ------------------------------------------------------------
function simulate_to_alpha_and_sample_once(
    N::Int,
    alpha0::Float64,
    end_alpha::Float64,
    alpha_increment::Float64,
    tau::Int,
    h::Float64,
    J::Float64;
    rng::AbstractRNG = Random.default_rng(),
    max_steps::Int = 10^9
)::Vector{Float64}

    X  = zeros(Int, N)
    S1 = zeros(Float64, N)      # S_1(i,t)
    Z  = fill(0.5, N)           # Z(i,t) = S1/S

    # evaporation factor
    evap = (tau == -1) ? 1.0 : exp(-1 / tau)

    S_prev = 0.0

    # clamp alpha range (safety)
    alpha = clamp(alpha0, 0.0, 1.0)
    end_alpha = clamp(end_alpha, 0.0, 1.0)

    for t in 1:max_steps
        # 1) sample X using current alpha
        @inbounds for i in 1:N
            p = decision_function(Z[i], alpha)
            X[i] = (rand(rng) < p) ? 1 : 0
        end

        # 2) pheromone deposit (Boltzmann weight)
        TP = calculate_pheromone(N, X, h, J)

        # 3) update S and S1 with evaporation
        S_curr = (t == 1) ? TP : (S_prev * evap + TP)

        if t == 1
            @inbounds for i in 1:N
                S1[i] = X[i] * TP
            end
        else
            @inbounds for i in 1:N
                S1[i] = S1[i] * evap + X[i] * TP
            end
        end

        # 4) update Z
        invS = 1.0 / S_curr
        @inbounds for i in 1:N
            Z[i] = S1[i] * invS
        end

        # 5) If we've reached the target alpha (the alpha used THIS step), sample and return
        if alpha >= end_alpha - 1e-15
            # M(i) = 2α(Z(i)-1/2)
            return 2.0 * end_alpha .* (Z .- 0.5)
        end

        # 6) anneal alpha for next step
        alpha = min(alpha + alpha_increment, end_alpha)

        S_prev = S_curr
    end

    error("Reached max_steps without sampling. Check alpha_increment/end_alpha/max_steps.")
end

# ------------------------------------------------------------
# Sampling wrapper:
#   Repeat simulate_to_alpha_and_sample_once 'samples' times,
#   and return flattened vector length N*samples.
#
# Signature MUST match main.jl:
#   sample_ants(N, alpha0, end_alpha, alpha_increment, tau, samples, h, J, seed)
# ------------------------------------------------------------
function sample_ants(
    N::Int,
    alpha0::Float64,
    end_alpha::Float64,
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
        @sync @distributed for s in 1:samples
            rng = MersenneTwister(seed + s - 1)
            M_samples[:, s] = simulate_to_alpha_and_sample_once(
                N, alpha0, end_alpha, alpha_increment, tau, h, J;
                rng=rng
            )
        end
    else
        for s in 1:samples
            rng = MersenneTwister(seed + s - 1)
            M_samples[:, s] = simulate_to_alpha_and_sample_once(
                N, alpha0, end_alpha, alpha_increment, tau, h, J;
                rng=rng
            )
        end
    end

    # flatten to Vector{Float64} of length N*samples
    return vec(Array(M_samples))
end

end # module
