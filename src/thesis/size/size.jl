using Random
using DataFrames
using CSV
using ProgressMeter
using Base.Threads

# --- パラメータ ---
const J = 0.1
const h = 0.001
const tau = 100.0
const N_list = [10, 100, 1000]
# Tを対数刻みで設定
const T_steps = exp10.(range(3, 6, length=20)) # 10^3 ~ 10^7 (元のコードに合わせて調整してください)
const n_samples = 100 # 各点あたりの試行回数

function check_success(N, T_total, J, h, tau)
    # 初期化の高速化: fillを使う
    X = rand(0:1, N)
    S1 = fill(1e-10, N)
    S0 = fill(1e-10, N)
    decay = exp(-1.0 / tau)
    
    # ループ回数を整数化
    limit_t = trunc(Int, T_total)
    
    # 定数計算をループ外へ
    J_norm = J / (N - 1)

    for t in 1:limit_t
        alpha = t / T_total
        
        # Xの更新 (ここは元のままですが、ここも最適化可能です)
        @inbounds for k in 1:N
            denom = S0[k] + S1[k]
            # ゼロ除算回避（念のため）
            Z = denom > 1e-20 ? S1[k] / denom : 0.5
            
            prob = (1.0 - alpha) * 0.5 + alpha * Z
            X[k] = rand() < prob ? 1 : 0
        end
        
        # --- 【ここが最大の修正点】 ---
        # 配列演算 (2 .* X .- 1) を避ける
        # sum(σ) = sum(2X - 1) = 2*sum(X) - N
        
        sum_X = sum(X) # これは高速でメモリ確保しません
        sum_sigma = 2 * sum_X - N # スカラー計算
        
        # エネルギー計算 (一時配列を作らない)
        # E = -h * Σσ - (J/N-1) * ((Σσ)^2 - N)/2
        E = -h * sum_sigma - J_norm * ((sum_sigma^2 - N) / 2.0)
        
        w = exp(-E)
        
        @inbounds for k in 1:N
            if X[k] == 1
                S1[k] = S1[k] * decay + w
                S0[k] = S0[k] * decay
            else
                S1[k] = S1[k] * decay
                S0[k] = S0[k] * decay + w
            end
        end
    end
    
    # 成功判定: 全て+1 (m=1.0) かどうか
    return sum(X) == N ? 1 : 0
end

println("サイズ依存性の計算を開始します...")
results = DataFrame(N=Int[], T=Float64[], success_rate=Float64[])

# マルチスレッドで計算
@showprogress for N in N_list
    for T_val in T_steps
        success_count = Atomic{Int}(0)
        @threads for _ in 1:n_samples
            s = check_success(N, T_val, J, h, tau)
            atomic_add!(success_count, s)
        end
        rate = success_count[] / n_samples
        push!(results, (N, T_val, rate))
    end
end

CSV.write("finite_size_data.csv", results)
println("finite_size_data.csv を保存しました。")