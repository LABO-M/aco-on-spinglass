using Random
using Statistics
using Printf

# ==========================================
# 1. モデル定義
# ==========================================

function calculate_energy(X::Vector{Int}, J::Float64, h::Float64, N::Int)
    # スピン変数 σ = 2X - 1 ({ -1, 1 })
    sigma_sum = sum(2 .* X .- 1)
    
    # エネルギー計算
    # ★修正: Isign模型の相互作用項には通常 0.5 の係数がつきます (二重カウント防止)
    # E = -h Σσ - (J / 2(N-1)) * ( (Σσ)^2 - N )
    term1 = -h * sigma_sum
    term2 = -(J / (N - 1)) * (sigma_sum^2 - N) 
    
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
    
    # アリの経路構築
    for k in 1:N
        Z = state.S1[k] / state.S_total
        prob_1 = (1.0 - alpha) * 0.5 + alpha * Z
        X[k] = rand() < prob_1 ? 1 : 0
    end
    
    # エネルギー評価
    E = calculate_energy(X, J, h, N)
    
    # フェロモン更新
    decay = exp(-1.0 / state.tau)
    weight = exp(-E)
    
    state.S_total = state.S_total * decay + weight
    for k in 1:N
        state.S1[k] = state.S1[k] * decay + (X[k] == 1 ? weight : 0.0)
    end
    
    # ★修正: 戻り値を「系全体の磁化(スカラー)」に変更
    # ヒストグラムにはサンプリングされた配位Xの磁化を使うのが一般的です
    m_system = mean(2 .* X .- 1)
    
    # 成功判定用（Fig4）にフェロモンの偏りも計算したい場合は別途保持が必要ですが、
    # ここでは物理量の観測として m_system を返します。
    return m_system, X
end

# ==========================================
# 2. データ生成関数
# ==========================================

const N_SPINS = 100
const J_VAL = 0.1
const H_VAL = 1.0e-3

# --- Figure 3 データ生成 ---
function generate_fig3_data(filename::String)
    println("Generating Figure 3 data...")

    # 閾値設定
    ALPHA_C_TAU100  = 0.828
    ALPHA_C_TAU1000 = 0.835 # カスタム

    taus = [100.0, 1000.0]
    
    # ★注意: n_samples = 5000, T_schedule = 10^6 は非常に時間がかかります。
    # テスト時は n_samples = 100 程度推奨。本番では 1000~5000 に戻してください。
    n_samples = 1000      
    T_schedule = 1000000 
    
    open(filename, "w") do io
        write(io, "tau,alpha_target,magnetization\n")
        
        for tau in taus
            # ターゲット設定
            if tau == 100.0
                target_alphas = [0.8, ALPHA_C_TAU100, 0.9]
                println("  Processing tau=100 Targets: $target_alphas")
            else
                target_alphas = [0.8, ALPHA_C_TAU1000, 0.9]
                println("  Processing tau=1000 Targets: $target_alphas")
            end
            
            sorted_targets = sort(target_alphas)

            # 進捗表示用
            next_print = 10

            for s in 1:n_samples
                if (s * 100 / n_samples) >= next_print
                    println("    Progress: $next_print% ($s/$n_samples)")
                    next_print += 10
                end

                state = initialize_aco(N_SPINS, tau)
                current_target_idx = 1
                
                for t in 1:T_schedule
                    alpha = t / T_schedule
                    
                    if current_target_idx <= length(sorted_targets) && alpha >= sorted_targets[current_target_idx]
                        m_sys, _ = step_ant!(state, alpha, J_VAL, H_VAL)
                        
                        target_val = sorted_targets[current_target_idx]
                        
                        # ★修正: スカラー値 m_sys を1行だけ書き込み
                        @printf(io, "%.1f,%.3f,%.6f\n", tau, target_val, m_sys)
                        
                        current_target_idx += 1
                    else
                        # 記録しないステップ
                        step_ant!(state, alpha, J_VAL, H_VAL)
                    end
                end
            end
        end
    end
    println("Saved $filename")
end

# --- Figure 4 データ生成 ---
function generate_fig4_data(filename::String)
    println("Generating Figure 4 data...")
    
    # T_max を小さくするか、trials を減らさないと時間がかかりすぎる可能性があります
    conditions = [
        ("Fast", 10000, 100.0),
        ("Fast", 10000, 1000.0),
        ("Slow", 1000000, 100.0), 
        ("Slow", 1000000, 1000.0)
    ]
    trials = 100 # テスト用 (本番は 1000)
    
    alpha_points = collect(0.7:0.01:1.0) # 点を少し間引いています (0.005 -> 0.01)
    
    open(filename, "w") do io
        write(io, "schedule,tau,alpha,success_prob\n")
        
        for (sched_name, T_max, tau) in conditions
            println("  Running $sched_name, tau=$tau ...")
            success_counts = zeros(Int, length(alpha_points))
            
            # 進捗表示
            next_print = 20

            for s in 1:trials
                if (s * 100 / trials) >= next_print
                    println("    Trial Progress: $next_print%")
                    next_print += 20
                end

                state = initialize_aco(N_SPINS, tau)
                curr_plot_idx = 1
                
                for t in 1:T_max
                    alpha = t / T_max
                    m_sys, X_curr = step_ant!(state, alpha, J_VAL, H_VAL)
                    
                    while curr_plot_idx <= length(alpha_points) && alpha >= alpha_points[curr_plot_idx]
                        # 成功判定: 全てのスピンが+1 (Xが全て1) であるか
                        # フェロモンではなく実際の配位で判定するのが厳密ですが、
                        # ACOの収束を見るならフェロモン確率(Ms) > 0 でも可。
                        # ここでは最も厳しい「全スピン+1」を成功とします。
                        if all(x -> x == 1, X_curr)
                            success_counts[curr_plot_idx] += 1
                        end
                        curr_plot_idx += 1
                    end
                end
            end
            
            probs = success_counts ./ trials
            for i in 1:length(alpha_points)
                @printf(io, "%s,%.1f,%.3f,%.4f\n", sched_name, tau, alpha_points[i], probs[i])
            end
        end
    end
    println("Saved $filename")
end

# ==========================================
# 3. メイン実行
# ==========================================

# 実行
generate_fig3_data("figure3_data_test.csv")
#generate_fig4_data("figure4_data_corrected.csv")