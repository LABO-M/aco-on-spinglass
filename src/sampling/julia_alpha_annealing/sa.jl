using Random
using CSV
using DataFrames
using LinearAlgebra
using Base.Threads
using Printf
using FilePathsBase

# --- パラメータ ---
n = 100
start_seed = 42
start_temp = 100.0
end_temp = 0.01
cooling_rate = (end_temp / start_temp) ^ (1 / 100_000)
iter = 100_000
num_samples = 2

# --- エネルギー関数 ---
function calculate_energy(spins::Vector{Int8}, h, J)
    return -sum(h .* spins) - sum(J .* (spins * spins')) / max(length(spins) - 1, 1)
end

# --- 単純なスピン反転（候補生成） ---
function generate_neighbor(spins::Vector{Int8})
    new_spins = copy(spins)
    idx = rand(1:length(spins))
    new_spins[idx] *= -1
    return new_spins
end

# --- SA 各サンプル実行関数（エネルギー系列を返す） ---
function run_single_sa(seed::Int)
    Random.seed!(seed)
    
    spins = rand(Bool, n)
    spins = Int8.(2 .* spins .- 1)
    h = ones(n) .* 0.001
    J = ones(n, n) .* 0.1
    J[diagind(J)] .= 0

    energy_series = Vector{Float64}(undef, iter)
    temp = start_temp
    current_energy = calculate_energy(spins, h, J)

    for i in 1:iter
        new_spins = generate_neighbor(spins)
        new_energy = calculate_energy(new_spins, h, J)
        ΔE = new_energy - current_energy

        if ΔE < 0 || rand() < exp(-ΔE / temp)
            spins = new_spins
            current_energy = new_energy
        end

        energy_series[i] = current_energy
        temp = max(temp * cooling_rate, end_temp)
    end

    return energy_series
end

# --- 並列実行＆データ収集 ---
all_energy = Vector{Vector{Float64}}(undef, num_samples)

@threads for sample_id in 1:num_samples
    seed = start_seed + sample_id
    energies = run_single_sa(seed)
    all_energy[sample_id] = energies
end

# --- データフレーム化 ---
df = DataFrame(step = 1:iter)
for i in 1:num_samples
    df[!, Symbol("sample$i")] = all_energy[i]
end

# --- 出力先ディレクトリ生成 ---
iter_str = @sprintf("%.1e", iter)
dir_path = "/home/mori-lab/shimizu/aco/data/ising/sa/seed$(start_seed)/iter$(iter_str)"
mkpath(dir_path)  # ディレクトリ作成

# --- ファイル名（科学記法でフォーマット） ---
filename = "energy.csv"
full_path = joinpath(dir_path, filename)

# --- CSV出力 ---
CSV.write(full_path, df)
println("SA energy series written to $full_path")
