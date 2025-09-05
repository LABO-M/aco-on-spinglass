module simulation

using Random
using CSV
using DataFrames
using LinearAlgebra
using Statistics
using Printf
using FilePathsBase

# --- 共通関数群 ---

# 外積は transpose を使って安全に
function calculate_energy(spins::AbstractVector, h, J::AbstractMatrix)
    n = length(spins)
    hvec = isa(h, Number) ? fill(h, n) : h
    S = spins * transpose(spins)  # 外積（n×n）
    return -sum(hvec .* spins) - sum(J .* S) / max(n - 1, 1)
end

function calculate_external_effect(determined_spins::AbstractVector, J::AbstractMatrix, new_spin_index::Int)
    return -mean(determined_spins .* J[1:length(determined_spins), new_spin_index])
end

function decision_probabilities(st1::AbstractVector, st0::AbstractVector, index::Int, alpha::Float64, beta::Float64, external_effect::Float64)
    z = st1[index] / (st1[index] + st0[index])
    return (z^alpha * exp(-beta * external_effect)) / (z^alpha * exp(-beta * external_effect) + (1 - z)^alpha * exp(beta * external_effect))
end

# =========================================================
#  - フェーズA：α=0 で 1000 回（フェロモン無し初期条件）
#  - フェーズB：α を 0.001 ずつ増加、各 α で 1000 回、α=1.0 で終了
#  - 磁場あり/なし両対応（magnetic=0.0 なら対称）
# =========================================================
function simulate(power::Int, tau::Float64, beta::Float64, seed::Int;
                  magnetic::Float64=0.0, interaction::Float64=0.0,
                  warmup_iters::Int=10000, alpha_step::Float64=0.001,
                  iters_per_alpha::Int=10000)

    rng = MersenneTwister(seed)

    n = 2 ^ power
    J = ones(n, n) * interaction
    J[diagind(J)] .= 0
    h = magnetic == 0.0 ? fill(0.0, n) : fill(magnetic, n)

    st1 = fill(1.0, n)
    st0 = fill(1.0, n)

    z_mean_series   = Float64[]
    z_var_series    = Float64[]
    m_mean_series   = Float64[]
    m_var_series    = Float64[]
    alpha_series    = Float64[]

    # --- フェーズA (α=0, early stopなし) ---
    alpha = 0.0
    zs = Float64[]
    ms = Float64[]

    for _ in 1:warmup_iters
        determined_spins = [2 * (rand(rng) < 0.5) - 1]
        for new_spins in 2:n
            external_effect = -magnetic + calculate_external_effect(determined_spins, J, new_spins)
            prob = decision_probabilities(st1, st0, new_spins, alpha, beta, external_effect)
            new_spin = 2 * (rand(rng) < prob) - 1
            push!(determined_spins, new_spin)
        end

        energy = calculate_energy(determined_spins, h, J)
        st1 .= st1 .* exp(-1 / tau) .+ exp(-energy) .* ((determined_spins .+ 1) ./ 2)
        st0 .= st0 .* exp(-1 / tau) .+ exp(-energy) .* ((-determined_spins .+ 1) ./ 2)

        z = mean(st1 ./ (st1 .+ st0))
        m = mean(determined_spins)

        push!(zs, z)
        push!(ms, m)
    end

    push!(z_mean_series, mean(zs))
    push!(z_var_series,  var(zs))
    push!(m_mean_series, mean(ms))
    push!(m_var_series,  var(ms))
    push!(alpha_series,  0.0)

    # --- フェーズB (各αで early stop 判定) ---
    alpha = 0.0
    while alpha < 1.0 + 1e-12
        alpha = round(alpha + alpha_step, digits=6)
        if alpha > 1.0
            alpha = 1.0
        end

        zs = Float64[]
        ms = Float64[]

        for _ in 1:iters_per_alpha
            determined_spins = [2 * (rand(rng) < 0.5) - 1]
            for new_spins in 2:n
                external_effect = -magnetic + calculate_external_effect(determined_spins, J, new_spins)
                prob = decision_probabilities(st1, st0, new_spins, alpha, beta, external_effect)
                new_spin = 2 * (rand(rng) < prob) - 1
                push!(determined_spins, new_spin)
            end

            energy = calculate_energy(determined_spins, h, J)
            st1 .= st1 .* exp(-1 / tau) .+ exp(-energy) .* ((determined_spins .+ 1) ./ 2)
            st0 .= st0 .* exp(-1 / tau) .+ exp(-energy) .* ((-determined_spins .+ 1) ./ 2)

            z = mean(st1 ./ (st1 .+ st0))
            m = mean(determined_spins)

            push!(zs, z)
            push!(ms, m)
        end

        m_mean = mean(ms)

        push!(z_mean_series, mean(zs))
        push!(z_var_series,  var(zs))
        push!(m_mean_series, m_mean)
        push!(m_var_series,  var(ms))
        push!(alpha_series,  alpha)

        if abs(m_mean) > 0.9 || alpha >= 1.0 - 1e-12
            break
        end
    end

    return z_mean_series, z_var_series, m_mean_series, m_var_series, alpha_series
end



function sampling(power::Int, tau::Float64, beta::Float64, seed::Int, magnetic::Float64, interaction::Float64)
    # --- αも含めて受け取る ---
    z_mean_series, z_var_series, spins_mean_series, spins_var_series, alpha_series = simulate(power, tau, beta, seed, magnetic=magnetic, interaction=interaction)

    iter = length(z_mean_series)
    df = DataFrame(step = 1:iter,
                   alpha = alpha_series,
                   z_mean = z_mean_series,
                   z_var = z_var_series,
                   spins_mean = spins_mean_series,
                   spins_var = spins_var_series)

    dir_switch = magnetic == 0.0 ? "symmetric" : "asymmetric"
    dir_path = "/home/mori-lab/shimizu/aco/data/ising/annealing/onesample/$(dir_switch)"
    mkpath(dir_path)

    filename_all = @sprintf("beta%.1e_n2^%d_tau%.1e_results.csv", beta, power, tau)
    full_path_all = joinpath(dir_path, filename_all)
    CSV.write(full_path_all, df)

    return filename_all
end

end