using Random
using DataFrames
using CSV
using ProgressMeter
using Base.Threads

# --- パラメータ ---
const J = 0.1
const h = 0.001
const tau = 100.0
const N_list = [50, 100, 1000]
# Tを対数刻みで設定
const T_steps = exp10.(range(3, 6, length=20)) # 10^3 ~ 10^5
const n_samples = 100 # 各点あたりの試行回数

function check_success(N, T_total, J, h, tau)
    X = rand(0:1, N)
    S1 = zeros(Float64, N) .+ 1e-10
    S0 = zeros(Float64, N) .+ 1e-10
    decay = exp(-1.0 / tau)
    
    for t in 1:trunc(Int, T_total)
        alpha = t / T_total
        for k in 1:N
            Z = S1[k] / (S0[k] + S1[k])
            prob = (1.0 - alpha) * 0.5 + alpha * Z
            X[k] = rand() < prob ? 1 : 0
        end
        
        E = -h * sum(2 .* X .- 1) - (J / (N - 1)) * ((sum(2 .* X .- 1)^2 - N) / 2)
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