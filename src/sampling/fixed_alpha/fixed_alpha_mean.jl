module fixed_alpha_mean

using Random
using CSV
using DataFrames
using LinearAlgebra
using Statistics
using Base.Threads
using Printf
using FilePathsBase
using Distributed
using SharedArrays

# --- 共通関数群 ---
function calculate_external_effect(determined_spins, interaction)
    return -mean(interaction .* determined_spins)
end

function decision_probabilities(beta, external_effect)
    return exp(-beta * external_effect) / (exp(-beta * external_effect) + exp(beta * external_effect))
end

function simulate(power, beta, seed, sample, magnetic, interaction)
    #seed
    rng = MersenneTwister(seed)

    # 初期化
    n = 2 ^ power
    J = ones(n, n) * interaction
    J[diagind(J)] .= 0
    h = ones(n) * magnetic
    determined_spins = [2 * (rand(rng) < 0.5) - 1]
    for new_spins in 2:n
        external_effect = -magnetic + calculate_external_effect(determined_spins, interaction)
        prob = decision_probabilities(beta, external_effect)
        new_spin = 2 * (rand(rng) < prob) - 1
        push!(determined_spins, new_spin)
    end
    return mean(determined_spins)
end

function simulate(power, beta, seed, sample, interaction)
    #seed
    rng = MersenneTwister(seed)

    # 初期化
    n = 2 ^ power
    J = ones(n, n) * interaction
    J[diagind(J)] .= 0
    determined_spins = [2 * (rand(rng) < 0.5) - 1]
    for new_spins in 2:n
        external_effect = calculate_external_effect(determined_spins, interaction)
        prob = decision_probabilities(beta, external_effect)
        new_spin = 2 * (rand(rng) < prob) - 1
        push!(determined_spins, new_spin)
    end
    return mean(determined_spins)
end


function sampling(power, beta, start_seed, sample, magnetic, interaction)
    spins_mean_array = SharedArray{Float64}(sample)

    @sync @distributed for i in 1:sample
        seed = start_seed + (i-1) * 100
        if magnetic == 0.0
            spins_mean_array[i] = simulate(power, beta, seed, sample, interaction)
        else
            spins_mean_array[i] = simulate(power, beta, seed, sample, magnetic, interaction)
        end
    end


    # --- 出力先ディレクトリ作成 ---
    dir_switch = magnetic == 0.0 ? "symmetric" : "asymmetric"
    dir_path = "/home/mori-lab/shimizu/aco/data/ising/zero_alpha/seed$(start_seed)/$(dir_switch)"
    mkpath(dir_path)

    # --- ファイル名と保存 ---
    filename_spins_mean = @sprintf("beta%.1e_sample%.1e_n2^%d_spins_mean.csv", beta, sample, power)
    full_path_spins_mean = joinpath(dir_path, filename_spins_mean)
    vec = Vector(spins_mean_array)
    df = DataFrame(sample = vec)
    CSV.write(full_path_spins_mean, df)

    return filename_spins_mean
end

end