using Random
using DataFrames
using CSV
using ProgressMeter

# --- パラメータ設定 ---
const N = 100           # システムサイズ
const J = 0.1           # 相互作用
const h = 0.001         # 外場
const T_max = 1000000     # アニーリングステップ数
const tau = 100.0       # 蒸発時定数
const n_trials = 10     # 取得するサンプル数（軌道数）

# --- 関数定義 ---
# エネルギー計算 (全結合Ising)
function calc_energy(X, J, h)
    mag = sum(2 .* X .- 1) # magnetization sum sigma_i
    # E = -h * sum(sigma) - J/(N-1) * sum(sigma_i * sigma_j)
    # sum(sigma_i * sigma_j) = ( (sum sigma)^2 - sum(sigma^2) ) / 2
    # sum(sigma^2) is always N since sigma^2 = 1
    interaction = ((mag^2 - length(X)) / 2)
    return -h * mag - (J / (length(X) - 1)) * interaction
end

function run_aco_trajectory(trial_id)
    # 初期化
    X = rand(0:1, N)
    S1 = zeros(Float64, N) .+ 1e-10 # ゼロ除算防止の微小値
    S0 = zeros(Float64, N) .+ 1e-10
    
    decay = exp(-1.0 / tau)
    
    # 記録用配列
    alphas = Float64[]
    magnetizations = Float64[]
    ids = Int[]
    
    # アニーリングループ
    for t in 1:T_max
        alpha = t / T_max
        
        # 1. 解の構成 (確率的選択)
        for k in 1:N
            denom = S0[k] + S1[k]
            Z = S1[k] / denom
            
            # 線形意思決定関数: P = (1-a)/2 + a*Z
            prob = (1.0 - alpha) * 0.5 + alpha * Z
            
            X[k] = rand() < prob ? 1 : 0
        end
        
        # 2. エネルギー計算
        E = calc_energy(X, J, h)
        w = exp(-E)
        
        # 3. フェロモン更新 (再帰式)
        # S(t) = S(t-1)*decay + exp(-E)
        for k in 1:N
            S1[k] = S1[k] * decay + (X[k] == 1 ? w : 0.0)
            S0[k] = S0[k] * decay + (X[k] == 0 ? w : 0.0)
        end
        
        # 4. データ記録 (データ量削減のため間引く)
        if t % 50 == 0 || t == 1
            m = sum(2 .* X .- 1) / N
            push!(alphas, alpha)
            push!(magnetizations, m)
            push!(ids, trial_id)
        end
    end
    
    return DataFrame(trial_id=ids, alpha=alphas, m=magnetizations)
end

# --- メイン処理 ---
println("軌道計算を開始します...")
df_all = DataFrame()
@showprogress for i in 1:n_trials
    df = run_aco_trajectory(i)
    append!(df_all, df)
end

CSV.write("trajectory_data.csv", df_all)
println("trajectory_data.csv を保存しました。")