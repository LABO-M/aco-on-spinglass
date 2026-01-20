using Random
using Statistics
using Printf

# ==========================================
# 1. 実験設定 (論文準拠)
# ==========================================
const N_SPINS = 100         # N [cite: 391]
const J_VAL   = 0.1         # J [cite: 391]
const H_VAL   = 0.001       # h [cite: 391]
const T_SLOW  = 1_000_000   # T (Slow Annealing) [cite: 388]
# テスト用に少なくしていますが、論文通りの滑らかなグラフには 1000 推奨 [cite: 389]
const N_TRIALS = 100        

# エネルギー計算: 論文 Eq.(1) 
# 通常のIsing模型の定義(係数0.5)と異なり、重複和をとる定義になっています。
# 高速化のため、総和 S = Σσ を用いて計算します。
function calc_energy_inf_range(sum_sigma::Int, N::Int, J::Float64, h::Float64)
    # 第1項: -h * Σσ
    term1 = -h * sum_sigma
    # 第2項: -(J / (N-1)) * Σ_{i≠j} σ_i σ_j
    # Σ_{i≠j} σ_i σ_j = (Σσ)^2 - Σ(σ^2) = S^2 - N
    term2 = -(J / (N - 1)) * (sum_sigma^2 - N)
    return term1 + term2
end

# 基底状態のエネルギー (h>0 なので 全スピン+1 が解) [cite: 76]
const E_GROUND = calc_energy_inf_range(N_SPINS, N_SPINS, J_VAL, H_VAL)

# ==========================================
# 2. ACO (α-Annealing) 実装
# ==========================================
# [cite: 101, 136, 385]
function run_aco_paper(tau::Float64, seed::Int)
    Random.seed!(seed)
    
    # 初期化: フェロモンはS1=S0からスタート (M=0)
    S_total = 1.0
    S1 = fill(0.5, N_SPINS) 
    decay = exp(-1.0 / tau) # 蒸発項 [cite: 90]
    
    success = false
    hit_step = -1
    
    for t in 1:T_SLOW
        # α-annealing schedule [cite: 385]
        alpha = t / T_SLOW
        
        # --- 1. アリによる解の構築 ---
        sum_sigma = 0
        current_X = Vector{Int}(undef, N_SPINS)
        all_magnetization_positive = true
        inv_S_total = 1.0 / S_total
        
        for i in 1:N_SPINS
            # 決定関数 Eq.(4, 5) [cite: 109]
            Z = S1[i] * inv_S_total
            
            # ACOの成功判定条件: 全ての M(i,t) > 0 
            # M = 2α(Z - 0.5) なので、Z > 0.5 と同義
            if Z <= 0.5; all_magnetization_positive = false; end
            
            prob = (1.0 - alpha) * 0.5 + alpha * Z
            prob = clamp(prob, 0.0, 1.0) # 数値誤差対策
            
            if rand() < prob
                current_X[i] = 1; sum_sigma += 1
            else
                current_X[i] = -1; sum_sigma -= 1
            end
        end
        
        # 成功判定
        if !success && all_magnetization_positive
            success = true
            hit_step = t
            # 論文のFig4のデータとしては、成功した時点でその試行は完了とみなせます
            break 
        end

        # --- 2. フェロモン更新 ---
        # エネルギー計算
        E = calc_energy_inf_range(sum_sigma, N_SPINS, J_VAL, H_VAL)
        # 論文では exp(-E) を加算 [cite: 86]
        weight = exp(-E)
        
        # Eq.(6) 近辺の更新則
        S_total = S_total * decay + weight
        for i in 1:N_SPINS
            incr = (current_X[i] == 1) ? weight : 0.0
            S1[i] = S1[i] * decay + incr
        end
    end
    return success, hit_step
end

# ==========================================
# 3. SA (Simulated Annealing) 実装
# ==========================================
# [cite: 394-403]
function run_sa_paper(tau::Float64, seed::Int)
    Random.seed!(seed)
    
    # 初期化 (ランダム)
    spins = rand((-1, 1), N_SPINS)
    sum_sigma = sum(spins)
    current_E = calc_energy_inf_range(sum_sigma, N_SPINS, J_VAL, H_VAL)
    
    success = false
    hit_step = -1
    
    # 最終逆温度 β = 2τ 
    final_beta = 2.0 * tau
    
    for t in 1:T_SLOW
        # 線形アニーリング [cite: 398]
        beta = (final_beta / T_SLOW) * t
        
        # 1 MCS (Monte Carlo Step) = N回の更新試行 [cite: 403]
        # (論文では N trial per MCS とあるためループします)
        for _ in 1:N_SPINS
            idx = rand(1:N_SPINS)
            s_i = spins[idx]
            
            # エネルギー差分計算 (無限レンジ模型の高速計算)
            # スピン s_i をフリップすると、総和 S は S - 2*s_i になる
            S_new = sum_sigma - 2 * s_i
            E_new = calc_energy_inf_range(S_new, N_SPINS, J_VAL, H_VAL)
            dE = E_new - current_E
            
            # Metropolis判定 [cite: 394]
            if dE < 0 || rand() < exp(-beta * dE)
                spins[idx] *= -1
                sum_sigma = S_new
                current_E = E_new
            end
        end
        
        # 成功判定: 基底状態を見つけたか 
        if !success && abs(current_E - E_GROUND) < 1e-5
            success = true
            hit_step = t
            break
        end
    end
    return success, hit_step
end

# ==========================================
# 4. メイン実行 & CSV出力
# ==========================================
function main()
    filename = "aco_sa_comparison_corrected.csv"
    println("Generating $filename ...")
    println("Settings: N=$N_SPINS, J=$J_VAL, h=$H_VAL, T=$T_SLOW, Trials=$N_TRIALS")

    open(filename, "w") do io
        # method: 手法名
        # hit_alpha: 成功した時点での進行度 (SAの場合は t/T)
        write(io, "method,tau,trial,success,hit_step,hit_alpha\n")
        
        tau_list = [100.0, 1000.0]
        
        for tau in tau_list
            println("  Processing τ = $tau")
            
            # --- ACO ---
            print("    Running ACO... ")
            for i in 1:N_TRIALS
                s, h_step = run_aco_paper(tau, i + 1000)
                # 成功したらその時点のα、失敗したら1.0として記録（プロット時に除外）
                h_alpha = s ? (h_step / T_SLOW) : 1.0
                write(io, "ACO,$tau,$i,$s,$h_step,$h_alpha\n")
            end
            println("Done.")

            # --- SA ---
            print("    Running SA...  ")
            for i in 1:N_TRIALS
                s, h_step = run_sa_paper(tau, i + 2000)
                h_alpha = s ? (h_step / T_SLOW) : 1.0
                write(io, "SA,$tau,$i,$s,$h_step,$h_alpha\n")
            end
            println("Done.")
        end
    end
    println("Finished! Saved to $filename")
end

main()