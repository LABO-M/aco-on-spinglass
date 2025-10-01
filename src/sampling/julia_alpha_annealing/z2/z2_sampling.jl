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

function simulate(power, tau, beta, seed, iter, magnetic, interaction)
    rng = MersenneTwister(seed)

    n = 2 ^ power
    if n < 100
        error("n must be ≥ 100 to record z[100]; got n=$(n)")
    end

    J = ones(n, n) * interaction
    J[diagind(J)] .= 0
    h = ones(n) * magnetic
    st1 = ones(n) .* exp(n * (magnetic + interaction)) .* (tau / 2)
    st0 = ones(n) .* exp(n * (magnetic + interaction)) .* (tau / 2)
    alpha = 0.0
    alpha_inc = 1.0 / iter

    z2_series = Float64[]
    z3_series = Float64[]
    z100_series = Float64[]

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

        push!(z2_series,   z[2])
        push!(z3_series,   z[3])
        push!(z100_series, z[100])

        alpha += alpha_inc
    end

    return z2_series, z3_series, z100_series
end


function simulate(power, tau, beta, seed, iter, interaction)
    rng = MersenneTwister(seed)

    n = 2 ^ power
    if n < 100
        error("n must be ≥ 100 to record z[100]; got n=$(n)")
    end

    J = ones(n, n) * interaction
    J[diagind(J)] .= 0
    st1 = ones(n) .* exp(n * interaction) .* (tau / 2)
    st0 = ones(n) .* exp(n * interaction) .* (tau / 2)
    alpha = 0.0
    alpha_inc = 1.0 / iter

    z2_series = Float64[]
    z3_series = Float64[]
    z100_series = Float64[]

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

        push!(z2_series,   z[2])
        push!(z3_series,   z[3])
        push!(z100_series, z[100])

        alpha += alpha_inc
    end

    return z2_series, z3_series, z100_series
end




function sampling(power::Int, tau::Float64, beta::Float64, start_seed::Int, iter::Int, sample::Int, magnetic::Float64, interaction::Float64)
    # z2, z3, z100 の3つを SharedArray で持つ
    z2_array   = SharedArray{Float64}(iter, sample)
    z3_array   = SharedArray{Float64}(iter, sample)
    z100_array = SharedArray{Float64}(iter, sample)

    @sync @distributed for i in 1:sample
        seed = start_seed + (i-1) * 1000
        local z2::Vector{Float64}
        local z3::Vector{Float64}
        local z100::Vector{Float64}
        if magnetic == 0.0
            z2, z3, z100 = simulate(power, tau, beta, seed, iter, interaction)
        else
            z2, z3, z100 = simulate(power, tau, beta, seed, iter, magnetic, interaction)
        end
        z2_array[:, i]   = z2
        z3_array[:, i]   = z3
        z100_array[:, i] = z100
    end

    # --- DataFrame 化（z2, z3, z100） ---
    df_z2   = DataFrame(step = 1:iter)
    df_z3   = DataFrame(step = 1:iter)
    df_z100 = DataFrame(step = 1:iter)
    for i in 1:sample
        df_z2[!,   Symbol("sample$i")] = z2_array[:, i]
        df_z3[!,   Symbol("sample$i")] = z3_array[:, i]
        df_z100[!, Symbol("sample$i")] = z100_array[:, i]
    end

    # --- 出力先ディレクトリ ---
    dir_switch = magnetic == 0.0 ? "symmetric" : "asymmetric"
    iter_str = @sprintf("%.1e", iter)
    dir_path = "/home/mori-lab/shimizu/aco/data/ising/annealing/seed$(start_seed)/iter$(iter_str)/$(dir_switch)"
    mkpath(dir_path)

    # --- ファイル名と保存 ---
    filename_z2   = @sprintf("beta%.1e_sample%.1e_n2^%d_tau%.1e_z2.csv",   beta, sample, power, tau)
    filename_z3   = @sprintf("beta%.1e_sample%.1e_n2^%d_tau%.1e_z3.csv",   beta, sample, power, tau)
    filename_z100 = @sprintf("beta%.1e_sample%.1e_n2^%d_tau%.1e_z100.csv", beta, sample, power, tau)

    CSV.write(joinpath(dir_path, filename_z2),   df_z2)
    CSV.write(joinpath(dir_path, filename_z3),   df_z3)
    CSV.write(joinpath(dir_path, filename_z100), df_z100)

    return filename_z2, filename_z3, filename_z100
end




end
