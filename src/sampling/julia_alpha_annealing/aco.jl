using Random
using CSV
using DataFrames
using LinearAlgebra
using Base.Threads
using Printf
using FilePathsBase

# --- パラメータの設定 ---
n = 100
tau = 100.0
evaporation_rate = exp(-1 / tau)
beta = 1
start_alpha = 0.0
start_seed = 42
iter = 100000
num_samples = 100  # サンプル数

# --- 共通関数群 ---
function initialize_spinglass(n, seed, iter, start_alpha)
    Random.seed!(seed)
    spins = rand(Bool,n)
    spins = 2 .* spins .- 1  # {+1, -1}

    J = randn(n, n)
    J = tril(J, -1) + tril(J, -1)'  # 対称行列
    h = ones(n) * 0.001
    st1 = ones(n) * 100.0
    st0 = ones(n) * 100.0
    energy = calculate_energy(spins, h, J)
    st1 += exp(-energy) .* ((spins .+ 1) ./ 2)
    st0 += exp(-energy) .* ((-spins .+ 1) ./ 2)

    alpha = start_alpha
    alpha_inc = (1 - start_alpha) / iter

    return h, J, st1, st0, alpha, alpha_inc
end

function calculate_energy(spins, h, J)
    return -sum(h .* spins) - sum(J .* (spins * spins')) / max(length(spins) - 1, 1)
end

function calculate_external_effect(determined_spins, J, new_spin_index)
    return -sum(determined_spins .* J[1:length(determined_spins), new_spin_index])
end

function decision_probabilities(st1, st0, index, alpha, beta, external_effect)
    z = st1[index] / (st1[index] + st0[index])
    return (z^alpha * exp(-beta * external_effect)) / (z^alpha * exp(-beta * external_effect) + (1 - z)^alpha * exp(beta * external_effect))
end

# --- 各サンプル分のエネルギー格納 ---
energy_all_aco = Vector{Vector{Float64}}(undef, num_samples)

@threads for sample_id in 1:num_samples
    Random.seed!(start_seed + sample_id)
    spins = rand(Bool, n)
    spins = 2 .* spins .- 1
    J = ones(n, n) * 0.1
    J[diagind(J)] .= 0
    h = ones(n) * 0.001
    st1 = ones(n) * 100.0
    st0 = ones(n) * 100.0
    energy = calculate_energy(spins, h, J)
    st1 += exp(-energy) .* ((spins .+ 1) ./ 2)
    st0 += exp(-energy) .* ((-spins .+ 1) ./ 2)

    alpha = 0.0
    alpha_inc = 1.0 / iter

    energy_series = Float64[]

    for i in 1:iter
        determined_spins = [2 * (rand() < 0.5) - 1]
        for new_spins in 2:n
            external_effect = calculate_external_effect(determined_spins, J, new_spins)
            prob = decision_probabilities(st1, st0, new_spins, alpha, beta, external_effect)
            new_spin = 2 * (rand() < prob) - 1
            push!(determined_spins, new_spin)
        end

        energy = calculate_energy(determined_spins, h, J)
        st1 .= st1 .* evaporation_rate .+ exp(-energy) .* ((determined_spins .+ 1) ./ 2)
        st0 .= st0 .* evaporation_rate .+ exp(-energy) .* ((-determined_spins .+ 1) ./ 2)
        alpha += alpha_inc

        push!(energy_series, energy)
    end
    energy_all_aco[sample_id] = energy_series
end

# --- データフレーム化 ---
df_aco = DataFrame(step = 1:iter)
for i in 1:num_samples
    df_aco[!, Symbol("sample$i")] = energy_all_aco[i]
end

# --- 出力先ディレクトリ生成 ---
iter_str = @sprintf("%.1e", iter)
dir_path = "/home/mori-lab/shimizu/aco/data/ising/annealing/seed$(start_seed)/iter$(iter_str)"
mkpath(dir_path)  # ディレクトリ作成

# --- ファイル名（科学記法でフォーマット） ---
filename = @sprintf("beta%.1e_sample%.1e_tau%.1e_energy.csv", beta, num_samples, tau)
full_path = joinpath(dir_path, filename)

# --- CSV出力 ---
CSV.write(full_path, df_aco)
println("ACO energy series written to aco_all_energy.csv")
