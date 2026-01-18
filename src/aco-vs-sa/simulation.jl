using Random
using Statistics
using Printf
using Dates
using LinearAlgebra

# --- 実験設定 (定数) ---
const N_SPINS = 100         # スピン数
const J_VAL = 0.1           # 相互作用係数
const H_VAL = 0.001         # 外部磁場
const MAX_EVALS = 200000    # 統一計算予算 (関数評価回数)
const N_TRIALS = 20         # 試行回数

# ACO設定
const N_ANTS = 20
const TAU = 1000.0          # 蒸発の時定数

# SA設定
const BETA_START = 0.0
const BETA_END = 5.0

# --- データ構造 ---
struct IsingModel
    n::Int
    J::Float64
    h::Float64
    J_mat::Matrix{Float64}
end

function IsingModel(n, J, h)
    mat = fill(J, n, n)
    for i in 1:n; mat[i,i] = 0.0; end
    return IsingModel(n, J, h, mat)
end

# エネルギー計算 (全計算)
function calc_energy(model::IsingModel, spins::Vector{Int})
    # spins ∈ {-1, 1}
    term_h = -model.h * sum(spins)
    term_J = -0.5 * dot(spins, model.J_mat * spins)
    return term_h + term_J
end

# 差分エネルギー計算 (SA用高速化)
function calc_delta_energy(model::IsingModel, spins::Vector{Int}, flip_idx::Int)
    s_i = spins[flip_idx]
    # interaction = sum(J_ij * s_j)
    interaction = dot(view(model.J_mat, :, flip_idx), spins)
    return 2.0 * s_i * (model.h + interaction)
end

# --- Simulated Annealing (SA) ---
function run_sa(model::IsingModel, max_evals::Int, seed::Int)
    Random.seed!(seed)
    t_start = time()
    
    current_spins = rand((-1, 1), model.n)
    current_energy = calc_energy(model, current_spins)
    best_energy = current_energy
    
    # 理論上の基底状態 (全て+1)
    ground_truth_energy = calc_energy(model, ones(Int, model.n))
    
    success = false
    hitting_time = NaN
    history = Float64[]
    
    # 履歴を間引いて記録するためのインターバル
    record_interval = max(1, max_evals ÷ 1000)

    for t in 1:max_evals
        beta = BETA_START + (BETA_END - BETA_START) * (t / max_evals)
        
        flip_idx = rand(1:model.n)
        dE = calc_delta_energy(model, current_spins, flip_idx)
        
        if dE < 0 || rand() < exp(-beta * dE)
            current_spins[flip_idx] *= -1
            current_energy += dE
            if current_energy < best_energy
                best_energy = current_energy
            end
        end
        
        # 成功判定
        if !success && abs(current_energy - ground_truth_energy) < 1e-5
            success = true
            hitting_time = t
        end
        
        if t % record_interval == 0
            push!(history, best_energy)
        end
    end
    
    t_end = time()
    return (success, hitting_time, best_energy, t_end - t_start, history)
end

# --- ACO (Alpha-Annealing) ---
function run_aco(model::IsingModel, max_evals::Int, seed::Int)
    Random.seed!(seed)
    t_start = time()
    
    ground_truth_energy = calc_energy(model, ones(Int, model.n))
    
    # フェロモン初期化 (微小値)
    s1 = fill(1e-10, model.n)
    s0 = fill(1e-10, model.n)
    
    best_energy = Inf
    success = false
    hitting_time = NaN
    history = Float64[]
    
    generations = max_evals ÷ N_ANTS
    evals_count = 0
    record_interval = max(1, max_evals ÷ 1000)
    
    for t in 1:generations
        # アニーリングスケジュール
        max_alpha = 1.0 - (1.0 / TAU)
        alpha = (t / generations) * max_alpha
        
        # フェロモン率 Z
        total_s = s1 .+ s0
        z = s1 ./ total_s
        
        # 意思決定確率 P(X=1)
        prob_one = @. (1.0 - alpha) * 0.5 + alpha * z
        clamp!(prob_one, 0.0, 1.0)
        
        # バッチごとのエネルギー保存用
        batch_energies = Vector{Float64}(undef, N_ANTS)
        batch_spins = Matrix{Int}(undef, N_ANTS, model.n)
        
        # アリの生成と評価
        for i in 1:N_ANTS
            # 各スピンを確率に基づいて決定
            for j in 1:model.n
                batch_spins[i, j] = rand() < prob_one[j] ? 1 : -1
            end
            
            e = calc_energy(model, batch_spins[i, :])
            batch_energies[i] = e
            
            if e < best_energy
                best_energy = e
                if !success && abs(best_energy - ground_truth_energy) < 1e-5
                    success = true
                    hitting_time = evals_count + i
                end
            end
        end
        
        evals_count += N_ANTS
        if evals_count % record_interval == 0
            push!(history, best_energy)
        end

        # フェロモン更新
        s1 .*= exp(-1.0 / TAU)
        s0 .*= exp(-1.0 / TAU)
        
        min_e = minimum(batch_energies)
        weights = exp.(-(batch_energies .- min_e))
        
        for i in 1:N_ANTS
            w = weights[i]
            for j in 1:model.n
                if batch_spins[i, j] == 1
                    s1[j] += w
                else
                    s0[j] += w
                end
            end
        end
    end
    
    t_end = time()
    return (success, hitting_time, best_energy, t_end - t_start, history)
end

# --- メイン実行とCSV出力 ---
function main()
    model = IsingModel(N_SPINS, J_VAL, H_VAL)
    
    println("Simulation Start: N=$(N_SPINS), MaxFEs=$(MAX_EVALS), Trials=$(N_TRIALS)")
    
    # 結果保存用ファイル
    open("results_summary.csv", "w") do io
        println(io, "Method,Trial,Success,HittingTime,FinalEnergy,WallClockTime")
        
        # SA実行
        println("Running SA...")
        for i in 1:N_TRIALS
            res = run_sa(model, MAX_EVALS, i)
            # CSV書き込み (NaNは空文字等にせずそのまま書く)
            println(io, "SA,$i,$(res[1]),$(res[2]),$(res[3]),$(res[4])")
            
            # 最後の試行だけ履歴を別途保存
            if i == N_TRIALS
                open("history_sa.csv", "w") do h_io
                    println(h_io, "Step,BestEnergy")
                    for (step, e) in enumerate(res[5])
                        println(h_io, "$step,$e")
                    end
                end
            end
        end
        
        # ACO実行
        println("Running ACO...")
        for i in 1:N_TRIALS
            res = run_aco(model, MAX_EVALS, i + 1000)
            println(io, "ACO,$i,$(res[1]),$(res[2]),$(res[3]),$(res[4])")
            
            if i == N_TRIALS
                open("history_aco.csv", "w") do h_io
                    println(h_io, "Step,BestEnergy")
                    for (step, e) in enumerate(res[5])
                        println(h_io, "$step,$e")
                    end
                end
            end
        end
    end
    println("Done. CSV files generated.")
end

main()