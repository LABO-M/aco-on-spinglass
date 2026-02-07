using Distributed
using CSV
using DataFrames
using Dates
using Plots
using Statistics
using StatsPlots

# プロセス追加 (未追加の場合のみ)
if nprocs() == 1
    addprocs(4) 
end

@everywhere begin
    using Random
    using SharedArrays

    # --- 物理モデル再定義 ---
    function calculate_boltzmann_weight_fast(N::Int, current_spin_sum::Int, h::Float64, J::Float64)
        interaction = (current_spin_sum^2 - N) / (N - 1)
        E = -h * current_spin_sum - J * interaction
        return exp(-E)
    end

    function decision_function(z::Float64, alpha::Float64)
        return alpha * (z - 0.5) + 0.5
    end

    # 比較用：単一のアニーリング/固定α試行
    function run_single_simulation(N::Int, T_max::Int, tau::Int, h::Float64, J::Float64; mode=:annealing, fixed_alpha=1.0)
        X = rand(0:1, N)
        S_scalar = 0.0
        Sm = zeros(Float64, N)
        decay = exp(-1.0 / tau)
        
        success_history = zeros(Int, T_max)
        # m_historyは今回成功率比較のみに使うため、必要最小限に留める
        
        for t in 1:T_max
            current_alpha = (mode == :annealing) ? (t / T_max) : fixed_alpha
            
            Zm = S_scalar > 1e-100 ? Sm ./ S_scalar : fill(0.5, N)
            probs = decision_function.(Zm, current_alpha)
            X .= rand(N) .< probs
            
            spin_sum = sum(2 .* X .- 1)
            weight = calculate_boltzmann_weight_fast(N, spin_sum, h, J)
            
            if t == 1
                S_scalar = weight
                Sm .= Float64.(X) .* weight
            else
                S_scalar = S_scalar * decay + weight
                Sm .= Sm .* decay .+ Float64.(X) .* weight
            end
            
            # 全てのスピンが正の方向(Zm > 0.5)を向いているかを判定
            if all(z -> z > 0.5, Zm) 
                success_history[t] = 1 
            end
        end
        return success_history
    end
end

# --- メイン処理 ---
function main_visualize()
    N = 100
    J = 0.1
    h = 0.001
    tau = 100
    T_max = 100000  # 100万は多すぎるため、動作確認用に10万に調整（必要に応じて戻してください）
    every = 1000    # 1000ステップごとに描画

    # ==========================================
    # 1. アニメーション作成 (1000ステップ毎)
    # ==========================================
    println("Creating animation (sampling every $every steps)...")
    X = rand(0:1, N)
    S_scalar = 0.0
    Sm = zeros(Float64, N)
    decay = exp(-1.0 / tau)
    
    # フレームを格納するアニメーションオブジェクトを先に定義
    anim = Animation()

    for t in 1:T_max
        alpha = t / T_max
        Zm = S_scalar > 1e-100 ? Sm ./ S_scalar : fill(0.5, N)
        probs = decision_function.(Zm, alpha)
        X .= rand(N) .< probs
        spin_sum = sum(2 .* X .- 1)
        weight = calculate_boltzmann_weight_fast(N, spin_sum, h, J)
        
        if t == 1
            S_scalar = weight; Sm .= Float64.(X) .* weight
        else
            S_scalar = S_scalar * decay + weight; Sm .= Sm .* decay .+ Float64.(X) .* weight
        end

        # 1000ステップごとにプロットしてフレームを追加
        if t % every == 0
            current_ms = Zm
            p_anim = histogram(current_ms, bins=-0.1:0.05:1.1, xlims=(-0.1, 1.1), ylims=(0, N),
                title="t=$t, α=$(round(alpha, digits=2))",
                xlabel="Z", ylabel="Count", legend=false, color=:skyblue, fillalpha=0.7)
            frame(anim, p_anim) # 明示的にフレームを追加
        end
    end
    gif(anim, "magnetization_evolution.gif", fps = 10)
    println("Saved magnetization_evolution.gif")

    # ==========================================
    # 2. 成功率の比較 (型安全な分散処理)
    # ==========================================
    println("Comparing Success Probabilities...")
    trials = 100
    T_comp = 1000000
    
    # 型を明示的に指定して分散和をとる
    results_anneal = @distributed (+) for _ in 1:trials
        run_single_simulation(N, T_comp, tau, h, J, mode=:annealing)
    end
    
    results_fixed = @distributed (+) for _ in 1:trials
        run_single_simulation(N, T_comp, tau, h, J, mode=:fixed, fixed_alpha=0.9)
    end

    # プロット
    p = plot(title="Success Probability Comparison", xlabel="Steps", ylabel="Success Rate", ylims=(0, 1.1))
    plot!(p, 1:T_comp, results_anneal ./ trials, label="Annealing (α: 0→1)", lw=2)
    plot!(p, 1:T_comp, results_fixed ./ trials, label="Fixed α (=1.0)", lw=2, linestyle=:dash)
    
    savefig(p, "comparison_success_prob.png")
    println("Saved comparison_success_prob.png")
end

main_visualize()