module fixed_alpha_mean

using Random
using CSV
using DataFrames
using LinearAlgebra
using Statistics
using Base.Threads
using Printf
using FilePathsBase

# --- 共通関数群 ---
function calculate_external_effect(determined_spins, J, new_spin_index)
    return -mean(determined_spins .* J[1:length(determined_spins), new_spin_index])
end

function decision_probabilities(beta, external_effect)
    return exp(-beta * external_effect) / (exp(-beta * external_effect) + exp(beta * external_effect))
end

function sampling(power, beta, start_seed, sample, magnetic, interaction, dir_switch)
    spins_mean_series = []

    @threads for sample_id in 1:sample
        Random.seed!(start_seed + sample_id)

        # 初期化
        n = 2 ^ power
        spins = rand(Bool, n)
        spins = 2 .* spins .- 1
        J = ones(n, n) * interaction
        J[diagind(J)] .= 0
        h = ones(n) * magnetic

        determined_spins = [2 * (rand() < 0.5) - 1]
        for new_spins in 2:n
            external_effect = -magnetic + calculate_external_effect(determined_spins, J, new_spins)
            prob = decision_probabilities(beta, external_effect)
            new_spin = 2 * (rand() < prob) - 1
            push!(determined_spins, new_spin)
        end
        push!(spins_mean_series, mean(determined_spins))

    end

    # --- 出力先ディレクトリ作成 ---
    dir_path = "/home/mori-lab/shimizu/aco/data/ising/zero_alpha/seed$(start_seed)/$(dir_switch)"
    mkpath(dir_path)

    # --- ファイル名と保存 ---
    filename_spins_mean = @sprintf("beta%.1e_sample%.1e_n2^%d_spins_mean.csv", beta, sample, power)
    full_path_spins_mean = joinpath(dir_path, filename_spins_mean)
    CSV.write(full_path_spins_mean, spins_mean_series)

    return filename_spins_mean
end

end