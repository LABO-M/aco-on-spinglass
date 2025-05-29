module annealing_mean

using Random
using CSV
using DataFrames
using LinearAlgebra
using Statistics
using Base.Threads
using Printf
using FilePathsBase

# --- 共通関数群 ---
function calculate_energy(spins, h, J)
    return -sum(h .* spins) - sum(J .* (spins * spins')) / max(length(spins) - 1, 1)
end

function calculate_external_effect(determined_spins, J, new_spin_index)
    return -mean(determined_spins .* J[1:length(determined_spins), new_spin_index])
end

function decision_probabilities(st1, st0, index, alpha, beta, external_effect)
    z = st1[index] / (st1[index] + st0[index])
    return (z^alpha * exp(-beta * external_effect)) / (z^alpha * exp(-beta * external_effect) + (1 - z)^alpha * exp(beta * external_effect))
end

function sampling(power, tau, beta, start_alpha, start_seed, iter, sample, magnetic, interaction, dir_switch)
    # --- mの時系列格納用 ---
    m_mean_series_all = Vector{Vector{Float64}}(undef, sample)
    spins_mean_series_all = Vector{Vector{Float64}}(undef, sample)

    @threads for sample_id in 1:sample
        rng = MersenneTwister(start_seed + sample_id)

        # 初期化
        n = 2 ^ power
        spins = rand(rng, Bool, n)
        spins = 2 .* spins .- 1
        J = ones(n, n) * interaction
        J[diagind(J)] .= 0
        h = ones(n) * magnetic
        st1 = ones(n) * 100.0
        st0 = ones(n) * 100.0
        energy = calculate_energy(spins, h, J)
        st1 += exp(-energy) .* ((spins .+ 1) ./ 2)
        st0 += exp(-energy) .* ((-spins .+ 1) ./ 2)

        alpha = 0.0
        alpha_inc = 1.0 / iter
        m_mean_series = Float64[]
        spins_mean_series = Float64[]

        for i in 1:iter
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

            z = st1 ./ (st1 .+ st0)
            m = 2 .* alpha .* (z .- 0.5)
            push!(m_mean_series, mean(m))
            push!(spins_mean_series, mean(determined_spins))

            alpha += alpha_inc
        end

        m_mean_series_all[sample_id] = m_mean_series
        spins_mean_series_all[sample_id] = spins_mean_series

    end

    # --- DataFrame 化 ---
    df_m_mean = DataFrame(step = 1:iter)
    for i in 1:sample
        df_m_mean[!, Symbol("sample$i")] = m_mean_series_all[i]
    end
    df_spin_mean = DataFrame(step = 1:iter)
    for i in 1:sample
        df_spin_mean[!, Symbol("sample$i")] = spins_mean_series_all[i]
    end

    # --- 出力先ディレクトリ作成 ---
    iter_str = @sprintf("%.1e", iter)
    dir_path = "/home/mori-lab/shimizu/aco/data/ising/annealing/seed$(start_seed)/iter$(iter_str)/$(dir_switch)"
    mkpath(dir_path)

    # --- ファイル名と保存 ---
    filename_m_mean = @sprintf("beta%.1e_sample%.1e_n2^%d_tau%.1e_modified_m_mean.csv", beta, sample, power, tau)
    full_path_m_mean = joinpath(dir_path, filename_m_mean)
    CSV.write(full_path_m_mean, df_m_mean)
    filename_spins_mean = @sprintf("beta%.1e_sample%.1e_n2^%d_tau%.1e_modified_spins_mean.csv", beta, sample, power, tau)
    full_path_spins_mean = joinpath(dir_path, filename_spins_mean)
    CSV.write(full_path_spins_mean, df_spin_mean)

    return filename_m_mean, filename_spins_mean
end

end