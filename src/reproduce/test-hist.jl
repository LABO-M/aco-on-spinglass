using Random
using Statistics
using Printf

# --- モデル定義 (論文のEq.1に基づき係数0.5なし) ---
function calculate_energy(X::Vector{Int}, J::Float64, h::Float64, N::Int)
    sigma_sum = sum(2 .* X .- 1)
    # 論文 Eq(1) は 1/2 の係数がない定義
    term1 = -h * sigma_sum
    term2 = - (J / (N - 1)) * (sigma_sum^2 - N) 
    return term1 + term2
end

mutable struct ACOState
    N::Int
    S_total::Float64
    S1::Vector{Float64}
    tau::Float64
end

function initialize_aco(N::Int, tau::Float64)
    S_init = 1.0
    return ACOState(N, S_init, fill(0.5 * S_init, N), tau)
end

function step_ant!(state::ACOState, alpha::Float64, J::Float64, h::Float64)
    N = state.N
    X = zeros(Int, N)
    
    for k in 1:N
        Z = state.S1[k] / state.S_total
        prob_1 = (1.0 - alpha) * 0.5 + alpha * Z
        X[k] = rand() < prob_1 ? 1 : 0
    end
    
    E = calculate_energy(X, J, h, N)
    
    decay = exp(-1.0 / state.tau)
    weight = exp(-E)
    
    state.S_total = state.S_total * decay + weight
    for k in 1:N
        state.S1[k] = state.S1[k] * decay + (X[k] == 1 ? weight : 0.0)
    end
    
    # 各サイトの磁化 M(i,t) をベクトルで返す
    Zs = state.S1 ./ state.S_total
    Ms = @. 2 * alpha * (Zs - 0.5)
    
    return Ms
end

# --- データ生成 ---
const N_SPINS = 100
const J_VAL = 0.1
const H_VAL = 1.0e-3

function generate_fig3_custom(filename::String)
    println("Generating Figure 3 data with separate alpha_c...")

    # ▼▼▼ ここで tau ごとの臨界点を指定してください ▼▼▼
    # 論文では tau=100 で 0.835 付近
    ALPHA_C_TAU100  = 0.835
    
    # tau=1000 用の値を指定 (例: 少しずらす、あるいは同じにする)
    ALPHA_C_TAU1000 = 0.84
    # ▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲

    taus = [100.0, 1000.0]
    
    # 論文再現のためにはサンプル数を多くする必要があります
    n_samples = 2000     
    T_schedule = 1000000 
    
    open(filename, "w") do io
        write(io, "tau,alpha_target,magnetization\n")
        
        for tau in taus
            # tau に応じてターゲットリストを切り替え
            if tau == 100.0
                # 左:0.8, 中:tau100用Critical, 右:0.9
                target_alphas = [0.8, ALPHA_C_TAU100, 0.9]
                println("Processing tau=100. Targets: $target_alphas")
            else
                # 左:0.8, 中:tau1000用Critical, 右:0.9
                target_alphas = [0.8, ALPHA_C_TAU1000, 0.9]
                println("Processing tau=1000. Targets: $target_alphas")
            end
            
            sorted_targets = sort(target_alphas)

            for s in 1:n_samples
                if s % 500 == 0; print("."); end
                state = initialize_aco(N_SPINS, tau)
                current_target_idx = 1
                
                for t in 1:T_schedule
                    alpha = t / T_schedule
                    
                    if current_target_idx <= length(sorted_targets) && alpha >= sorted_targets[current_target_idx]
                        Ms = step_ant!(state, alpha, J_VAL, H_VAL)
                        
                        target_val = sorted_targets[current_target_idx]
                        
                        # ベクトル全要素を書き込み
                        for m_i in Ms
                            @printf(io, "%.1f,%.3f,%.6f\n", tau, target_val, m_i)
                        end
                        
                        current_target_idx += 1
                    else
                        step_ant!(state, alpha, J_VAL, H_VAL)
                    end
                end
            end
            println(" Done tau=$tau")
        end
    end
    println("Saved $filename")
end

# 実行
generate_fig3_custom("figure3_data_custom.csv")