module simulation

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

function simulate(power, tau, beta, start_alpha, seed, iter, sample, magnetic, interaction)
    #seed
    rng = MersenneTwister(seed)

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
    z_var_series = Float64[]
    spins_var_series = Float64[]

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
        push!(z_var_series, var(z))
        push!(spins_var_series, var(determined_spins))
        alpha += alpha_inc

    end

    return z_var_series, spins_var_series
end

function simulate(power, tau, beta, start_alpha, seed, iter, sample, interaction)
    #seed
    rng = MersenneTwister(seed)

    # 初期化
    n = 2 ^ power
    spins = rand(rng, Bool, n)
    spins = 2 .* spins .- 1
    J = ones(n, n) * interaction
    J[diagind(J)] .= 0
    st1 = ones(n) * 100.0
    st0 = ones(n) * 100.0
    energy = calculate_energy(spins, 0, J)
    st1 += exp(-energy) .* ((spins .+ 1) ./ 2)
    st0 += exp(-energy) .* ((-spins .+ 1) ./ 2)
    alpha = 0.0
    alpha_inc = 1.0 / iter
    z_var_series = Float64[]
    spins_var_series = Float64[]

    for i in 1:iter
        determined_spins = [2 * (rand(rng) < 0.5) - 1]

        for new_spins in 2:n
            external_effect = calculate_external_effect(determined_spins, J, new_spins)
            prob = decision_probabilities(st1, st0, new_spins, alpha, beta, external_effect)
            new_spin = 2 * (rand(rng) < prob) - 1
            push!(determined_spins, new_spin)
        end

        energy = calculate_energy(determined_spins, 0, J)
        st1 .= st1 .* exp(-1 / tau) .+ exp(-energy) .* ((determined_spins .+ 1) ./ 2)
        st0 .= st0 .* exp(-1 / tau) .+ exp(-energy) .* ((-determined_spins .+ 1) ./ 2)
        z = st1 ./ (st1 .+ st0)
        push!(z_var_series, var(z))
        push!(spins_var_series, var(determined_spins))
        alpha += alpha_inc

    end

    return z_var_series, spins_var_series
end

function sampling(power::Int, tau::Float64, beta::Float64, start_alpha::Float64, start_seed::Int, iter::Int, sample::Int, magnetic::Float64, interaction::Float64)
    z_var_array = SharedArray{Float64}(iter, sample)
    spins_var_array = SharedArray{Float64}(iter, sample)

    @sync @distributed for i in 1:sample
        seed = start_seed + (i-1) * 1000
        if magnetic == 0.0
            z_var_array[:, i], spins_var_array[:, i] = simulate(power, tau, beta, start_alpha, seed, iter, sample, interaction)
        else
            z_var_array[:, i], spins_var_array[:, i] = simulate(power, tau, beta, start_alpha, seed, iter, sample, magnetic, interaction)
        end
    end

    # --- DataFrame 化 ---
    df_z_var = DataFrame(step = 1:iter)
    for i in 1:sample
        df_z_var[!, Symbol("sample$i")] = z_var_array[:, i]
    end
    df_spin_var = DataFrame(step = 1:iter)
    for i in 1:sample
        df_spin_var[!, Symbol("sample$i")] = spins_var_array[:, i]
    end

    # --- 出力先ディレクトリ作成 ---
    dir_switch = magnetic == 0.0 ? "symmetric" : "asymmetric"
    iter_str = @sprintf("%.1e", iter)
    dir_path = "/home/mori-lab/shimizu/aco/data/ising/annealing/seed$(start_seed)/iter$(iter_str)/$(dir_switch)"
    mkpath(dir_path)

    # --- ファイル名と保存 ---
    filename_z_var = @sprintf("beta%.1e_sample%.1e_n2^%d_tau%.1e_z_var.csv", beta, sample, power, tau)
    full_path_z_var = joinpath(dir_path, filename_z_var)
    CSV.write(full_path_z_var, df_z_var)
    filename_spins_var = @sprintf("beta%.1e_sample%.1e_n2^%d_tau%.1e_spins_var.csv", beta, sample, power, tau)
    full_path_spins_var = joinpath(dir_path, filename_spins_var)
    CSV.write(full_path_spins_var, df_spin_var)

    return filename_z_var, filename_spins_var

end

end
