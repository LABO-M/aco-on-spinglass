using Distributed
using CSV
using DataFrames
using Dates

# プロセスが足りない場合は追加（コマンドライン引数で指定しなかった場合用）
if nprocs() == 1
    addprocs(4) 
end

@everywhere begin
    using Random
    using SharedArrays
    using Statistics
    
    # --- 定数と物理モデル ---
    # 論文式(1)に基づく全結合イジング模型のエネルギー計算 (O(N)高速化)
    function calculate_boltzmann_weight_fast(N::Int, current_spin_sum::Int, h::Float64, J::Float64)
        # 相互作用項: -1/(N-1) * Σ_{i≠j} J Si Sj
        # Σ_{i≠j} Si Sj = (Σ Si)^2 - Σ(Si^2) = (Σ Si)^2 - N
        interaction = (current_spin_sum^2 - N) / (N - 1)
        E = -h * current_spin_sum - J * interaction
        return exp(-E)
    end

    # 決定関数 (論文式(5))
    function decision_function(z::Float64, alpha::Float64)
        return alpha * (z - 0.5) + 0.5
    end

    # --- Figure 3用: 定常状態分布のサンプリング ---
    function sample_stationary_distribution(N::Int, alpha::Float64, tau::Int, h::Float64, J::Float64, steps::Int, samples_needed::Int)
        X = rand(0:1, N)
        S_scalar = 0.0
        Sm = zeros(Float64, N)
        decay = exp(-1.0 / tau)
        
        # 結果格納用
        collected_m = Float64[]
        
        # 緩和時間（適当にtauの5倍）
        burn_in = 5 * tau
        
        for t in 1:(burn_in + steps)
            # 現在の分布に基づく確率計算
            Zm = S_scalar > 1e-100 ? Sm ./ S_scalar : fill(0.5, N)
            probs = decision_function.(Zm, alpha)
            
            # アリの選択更新
            X .= rand(N) .< probs
            
            # 重み計算
            spin_sum = sum(2 .* X .- 1)
            weight = calculate_boltzmann_weight_fast(N, spin_sum, h, J)
            
            # フェロモン更新
            if t == 1
                S_scalar = weight
                Sm .= X .* weight
            else
                S_scalar = S_scalar * decay + weight
                Sm .= Sm .* decay .+ X .* weight
            end
            
            # サンプリング（緩和後、かつ相関を切るため間引き）
            if t > burn_in && (t % 10 == 0)
                # 磁化 m = 2α(z - 0.5)
                # 論文のヒストグラムは「個々のスピンの磁化」または「システム全体の磁化」ですが
                # 全結合で均質なので、全サイトのmを集めて分布を作ります。
                current_ms = 2 * alpha .* (Zm .- 0.5)
                append!(collected_m, current_ms)
            end
            
            if length(collected_m) >= samples_needed * N
                break
            end
        end
        
        return collected_m
    end

    # --- Figure 4用: アニーリング成功確率 ---
    function run_annealing(N::Int, T_max::Int, tau::Int, h::Float64, J::Float64)
        X = rand(0:1, N)
        S_scalar = 0.0
        Sm = zeros(Float64, N)
        decay = exp(-1.0 / tau)
        
        # 線形スケジュール: alpha = t / T
        alpha_step = 1.0 / T_max
        alpha = 0.0
        
        # 履歴記録用（1:成功, 0:失敗）
        # メモリ節約のため、Int8またはBitVectorでも良いが、集計用にInt
        success_history = zeros(Int, T_max)
        
        for t in 1:T_max
            # Zm計算
            Zm = S_scalar > 1e-100 ? Sm ./ S_scalar : fill(0.5, N)
            
            # 決定
            probs = decision_function.(Zm, alpha)
            X .= rand(N) .< probs
            
            # 更新
            spin_sum = sum(2 .* X .- 1)
            weight = calculate_boltzmann_weight_fast(N, spin_sum, h, J)
            
            if t == 1
                S_scalar = weight
                Sm .= X .* weight
            else
                S_scalar = S_scalar * decay + weight
                Sm .= Sm .* decay .+ X .* weight
            end
            
            # 成功判定: 全てのiで M(i) > 0  <=>  全てのiで Zm[i] > 0.5
            if all(Zm .> 0.5)
                success_history[t] = 1
            end
            
            alpha += alpha_step
        end
        return success_history
    end
end

function main()
    println("Simulation started at $(now())")
    
    # パラメータ (論文に準拠) [cite: 391, 406]
    N = 100
    J = 0.1
    h = 0.001 
    
    # ==========================================
    # 1. Figure 3 データ生成 (定常状態分布)
    # ==========================================
    println("Generating data for Figure 3...")
    
    # 論文 [cite: 303, 304] より J=0.1でのalpha_cは約0.83...
    # tau=100での厳密なalpha_cは0.83515
    
    
    conditions_fig3 = [
        (tau=100, alpha=0.85), 
        (tau=100, alpha=0.9087), 
        (tau=100, alpha=0.95),
        (tau=1000, alpha=0.85), 
        (tau=1000, alpha=0.915), 
        (tau=1000, alpha=0.95)
    ]
    
    df_fig3 = DataFrame(tau=Int[], alpha=Float64[], m_value=Float64[])
    
    # データ収集ループ
    for cond in conditions_fig3
        println("  Sampling: tau=$(cond.tau), alpha=$(cond.alpha)")
        # サンプル数: デモ用に少なめに設定。論文並みに綺麗にするなら samples_needed=100000 程度推奨
        ms = sample_stationary_distribution(N, cond.alpha, cond.tau, h, J, 50000, 5000)
        
        # DataFrameに追加 (メモリ節約のため間引いて保存しても良い)
        # ここでは最初の5000個だけ保存
        save_limit = min(length(ms), 5000)
        temp_df = DataFrame(tau=cond.tau, alpha=cond.alpha, m_value=ms[1:save_limit])
        append!(df_fig3, temp_df)
    end
    
    CSV.write("aco_fig3_data.csv", df_fig3)
    println("Saved aco_fig3_data.csv")

    # ==========================================
    # 2. Figure 4 データ生成 (成功確率)
    # ==========================================
    println("Generating data for Figure 4...")
    
    # アニーリング速度 T (論文: Fast=10^4, Slow=10^6)
    # デモ用に T_slow = 5*10^4 にしています。
    # 論文の結果(成功率1.0)を再現するには T_slow = 1000000 に戻してください。
    T_fast = 10000
    T_slow = 1000000  # ★重要: 時間がかかるので減らしています
    
    trials = 50     # ★重要: 論文は1000回ですが、デモ用に50回
    
    scenarios = [
        (tau=100, T=T_fast, label="tau=100, Fast"),
        (tau=1000, T=T_fast, label="tau=1000, Fast"),
        (tau=100, T=T_slow, label="tau=100, Slow"),
        (tau=1000, T=T_slow, label="tau=1000, Slow") # 時間があれば追加
    ]
    
    df_fig4 = DataFrame(alpha=Float64[], success_prob=Float64[], label=String[])
    
    for scen in scenarios
        println("  Annealing: $(scen.label) (T=$(scen.T))")
        
        # 並列計算 (pmapでも良いが、和を取るだけなので簡易的に実装)
        results = @distributed (+) for _ in 1:trials
            run_annealing(N, scen.T, scen.tau, h, J)
        end
        
        prob_curve = results ./ trials
        
        # alpha軸
        alpha_axis = range(0.0, 1.0, length=scen.T)
        
        # CSVサイズ削減のため、プロットに必要な点数(例えば500点)に間引く
        step_size = max(1, div(scen.T, 500))
        idx = 1:step_size:scen.T
        
        temp_df = DataFrame(
            alpha = collect(alpha_axis[idx]),
            success_prob = collect(prob_curve[idx]),
            label = scen.label
        )
        append!(df_fig4, temp_df)
    end
    
    CSV.write("aco_fig4_data.csv", df_fig4)
    println("Saved aco_fig4_data.csv")
    
    println("All simulations finished at $(now())")
end

main()