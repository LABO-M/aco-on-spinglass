using Random
using DataFrames
using CSV
using ProgressMeter
using Base.Threads

# --- パラメータ ---
const N = 100
const J = 0.1
const h = 0.001

# グリッド設定
const T_list = exp10.(range(3, 7, length=20)) # 10^3 ~ 10^5.5
const tau_list = exp10.(range(1, 6, length=20)) # 10^1 ~ 10^4
const n_samples = 50 # 計算時間を考慮して少し減らす

# 関数は実験2と同じものを使用 (再定義は省略可だが、単体動作のため書いておく)
function check_success_hm(N, T_total, J, h, tau)
    X = rand(0:1, N)
    S1 = zeros(Float64, N) .+ 1e-10
    S0 = zeros(Float64, N) .+ 1e-10
    decay = exp(-1.0 / tau)
    
    T_int = trunc(Int, T_total)
    for t in 1:T_int
        alpha = t / T_total
        for k in 1:N
            Z = S1[k] / (S0[k] + S1[k])
            prob = (1.0 - alpha) * 0.5 + alpha * Z
            X[k] = rand() < prob ? 1 : 0
        end
        # 高速化のためmag計算を展開
        mag = sum(2 .* X .- 1)
        E = -h * mag - (J / (N - 1)) * ((mag^2 - N) / 2)
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
    return sum(X) == N ? 1 : 0
end

println("ヒートマップデータの計算を開始します...")
results = DataFrame(tau=Float64[], T=Float64[], success_rate=Float64[])

# 総組み合わせ数
total_comb = length(T_list) * length(tau_list)
p = Progress(total_comb)

for tau_val in tau_list
    for T_val in T_list
        success_count = Atomic{Int}(0)
        @threads for _ in 1:n_samples
            s = check_success_hm(N, T_val, J, h, tau_val)
            atomic_add!(success_count, s)
        end
        push!(results, (tau_val, T_val, success_count[] / n_samples))
        next!(p)
    end
end

CSV.write("heatmap_data.csv", results)
println("heatmap_data.csv を保存しました。")