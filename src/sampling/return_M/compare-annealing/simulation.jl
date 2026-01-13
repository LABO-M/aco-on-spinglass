module Simulation

using Random
using Distributed
using SharedArrays

export simulate, sample_M

# -----------------------------
# decision functions
# -----------------------------
@inline function p_linear(z::Float64, alpha::Float64)::Float64
    zz = clamp(z, 0.0, 1.0)
    return (1.0 - alpha) * 0.5 + alpha * zz
end

@inline function p_sigmoid(z::Float64, alpha::Float64, beta::Float64)::Float64
    # α=0 なら 0.5
    zz = clamp(z, 0.0, 1.0)
    u = 2.0 * zz - 1.0
    x = beta * alpha * u
    return 1.0 / (1.0 + exp(-x))
end

@inline function decision_prob(z::Float64, alpha::Float64,
                              decision::Symbol, beta::Float64)::Float64
    if decision === :linear
        return p_linear(z, alpha)
    elseif decision === :sigmoid
        return p_sigmoid(z, alpha, beta)
    else
        error("Unknown decision=$decision (use :linear or :sigmoid)")
    end
end

# -----------------------------
# pheromone: exp(-E)
# E = -h sum s_i - (J/(N-1)) sum_{i<j} s_i s_j,  s_i=2X_i-1
# -----------------------------
function pheromone_weight(N::Int, X::Vector{Int}, h::Float64, J::Float64)::Float64
    s = 2 .* X .- 1
    E = -h * sum(s)
    coef = J / (N - 1)
    @inbounds for i in 1:(N-1)
        si = s[i]
        for j in (i+1):N
            E += -coef * si * s[j]
        end
    end
    return exp(-E)
end

# -----------------------------
# annealing: linear or two_stage
# -----------------------------
@inline function next_alpha(alpha::Float64, end_alpha::Float64,
                           schedule::Symbol, inc1::Float64,
                           alpha_mid::Float64, inc2::Float64)::Float64
    if schedule === :linear
        return min(alpha + inc1, end_alpha)
    elseif schedule === :two_stage
        inc = (alpha < alpha_mid) ? inc1 : inc2
        return min(alpha + inc, end_alpha)
    else
        error("Unknown schedule=$schedule (use :linear or :two_stage)")
    end
end

# -----------------------------
# simulate once -> return M (length N)
# M(i) = 2α(Z(i)-1/2)
# -----------------------------
function simulate(
    N::Int,
    alpha0::Float64,
    end_alpha::Float64,
    alpha_inc1::Float64,
    tau::Int,
    h::Float64,
    J::Float64;
    burn_in::Int = 1000,
    max_steps::Int = 10^9,
    rng::AbstractRNG = Random.default_rng(),
    decision::Symbol = :linear,
    beta::Float64 = 5.0,
    schedule::Symbol = :linear,
    alpha_mid::Float64 = 0.75,
    alpha_inc2::Float64 = 1e-7
)::Vector{Float64}

    X  = zeros(Int, N)
    S1 = zeros(Float64, N)
    Z  = fill(0.5, N)

    evap = (tau == -1) ? 1.0 : exp(-1 / tau)

    S_prev = 0.0
    alpha = alpha0
    reached = false
    burn_count = 0

    for t in 1:max_steps
        @inbounds for i in 1:N
            p = decision_prob(Z[i], alpha, decision, beta)
            X[i] = (rand(rng) < p) ? 1 : 0
        end

        TP = pheromone_weight(N, X, h, J)

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

        @inbounds for i in 1:N
            Z[i] = S1[i] / S_curr
        end

        if !reached
            if alpha >= end_alpha
                reached = true
                alpha = end_alpha
                burn_count = 0
            else
                alpha = next_alpha(alpha, end_alpha, schedule,
                                   alpha_inc1, alpha_mid, alpha_inc2)
            end
        else
            burn_count += 1
            if burn_count >= burn_in
                return 2.0 * alpha .* (Z .- 0.5)
            end
        end

        S_prev = S_curr
    end

    error("max_steps reached without sampling")
end

# -----------------------------
# sample many -> concatenate (length N*samples)
# -----------------------------
function sample_M(
    N::Int,
    alpha0::Float64,
    end_alpha::Float64,
    alpha_inc1::Float64,
    tau::Int,
    samples::Int,
    h::Float64,
    J::Float64,
    seed::Int;
    burn_in::Int = 1000,
    decision::Symbol = :linear,
    beta::Float64 = 5.0,
    schedule::Symbol = :linear,
    alpha_mid::Float64 = 0.75,
    alpha_inc2::Float64 = 1e-7
)::Vector{Float64}

    M_samples = SharedArray{Float64}(N, samples)

    if nworkers() > 1
        @sync @distributed for s in 1:samples
            rng = MersenneTwister(seed + s - 1)
            M_samples[:, s] = simulate(N, alpha0, end_alpha, alpha_inc1, tau, h, J;
                burn_in=burn_in, rng=rng,
                decision=decision, beta=beta,
                schedule=schedule, alpha_mid=alpha_mid, alpha_inc2=alpha_inc2)
        end
    else
        for s in 1:samples
            rng = MersenneTwister(seed + s - 1)
            M_samples[:, s] = simulate(N, alpha0, end_alpha, alpha_inc1, tau, h, J;
                burn_in=burn_in, rng=rng,
                decision=decision, beta=beta,
                schedule=schedule, alpha_mid=alpha_mid, alpha_inc2=alpha_inc2)
        end
    end

    return vcat(M_samples...)
end

end # module
