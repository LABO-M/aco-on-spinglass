using Distributed
using CSV
using DataFrames
using Dates
using Plots
using Statistics
using StatsPlots
using Plots.Measures # 余白調整のために追加

if nprocs() == 1
    addprocs(4) 
end

@everywhere begin
    using Random
    using SharedArrays

    function calculate_boltzmann_weight_fast(N::Int, current_spin_sum::Int, h::Float64, J::Float64)
        interaction = (current_spin_sum^2 - N) / (N - 1)
        E = -h * current_spin_sum - J * interaction
        return exp(-E)
    end

    function decision_function(z::Float64, alpha::Float64)
        return alpha * (z - 0.5) + 0.5
    end

    function run_single_simulation(N::Int, T_max::Int, tau::Int, h::Float64, J::Float64; mode=:annealing, fixed_alpha=1.0)
        X = rand(0:1, N)
        S_scalar = 0.0
        Sm = zeros(Float64, N)
        decay = exp(-1.0 / tau)
        success_history = zeros(Int, T_max)
        
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
            if all(z -> z > 0.5, Zm) success_history[t] = 1 end
        end
        return success_history
    end
end

function main_visualize()
    N = 100
    J = 0.1
    h = 0.001
    tau = 100
    T_max = 300000 # 動作確認のため一旦短縮していますが、必要に応じて戻してください
    every = 1000    
    fixed_alphas = [0.8, 0.9]

    # ==========================================
    # 1. アニメーション作成 (並列比較)
    # ==========================================
    println("Creating MP4 animation...")
    
    states = [ (mode=m, val=v, Sm=fill(500.0, N), S=1000.0) for (m, v) in [(:anneal, 0.0), (:fixed, 0.8), (:fixed, 0.95)] ]
    anim = Animation()

    for t in 1:T_max
        plots_step = []
        for i in 1:length(states)
            mode, alpha_val, Sm, S_scalar = states[i].mode, states[i].val, states[i].Sm, states[i].S
            curr_alpha = (mode == :anneal) ? (t / T_max) : alpha_val
            
            Zm = S_scalar > 1e-100 ? Sm ./ S_scalar : fill(0.5, N)
            X = rand(N) .< decision_function.(Zm, curr_alpha)
            weight = calculate_boltzmann_weight_fast(N, sum(2 .* X .- 1), h, J)
            
            # 更新
            decay = exp(-1.0 / tau)
            if t == 1
                new_S = weight; new_Sm = Float64.(X) .* weight
            else
                new_S = S_scalar * decay + weight
                new_Sm = Sm .* decay .+ Float64.(X) .* weight
            end
            states[i] = (mode=mode, val=alpha_val, Sm=new_Sm, S=new_S)

            if t % every == 0
                title_str = mode == :anneal ? "Annealing (α→$(round(curr_alpha, digits=2)))" : "α=$alpha_val"
                p = histogram(Zm, bins=-0.05:0.05:1.05, xlims=(-0.05, 1.05), ylims=(0, N),
                              title=title_str, xlabel="Z", ylabel="Count", legend=false, 
                              color=:skyblue, 
                              bottom_margin=8Plots.mm, # 下側の余白を確保
                              left_margin=5Plots.mm)   # 左側の余白を確保
                push!(plots_step, p)
            end
        end
        
        if t % every == 0
            # 全体のレイアウトでも余白を調整
            combined_p = plot(plots_step..., layout=(1, 3), size=(1200, 400), margin=3Plots.mm)
            frame(anim, combined_p)
        end
    end
    mp4(anim, "30_evolution_comparison.mp4", fps = 15)
    println("Saved evolution_comparison.mp4")

    # ==========================================
    # 2. 成功率の比較
    # ==========================================
    println("Comparing Success Probabilities...")
    trials = 100 
    T_comp = 1000000
    
    p_comp = plot(title="Success Rate Comparison", xlabel="Steps", ylabel="Rate", 
                  ylims=(0, 1.1), legend=:bottomright, 
                  bottom_margin=5Plots.mm, left_margin=5Plots.mm)

    # アニーリング
    println("Processing Annealing...")
    res_anneal = @distributed (+) for _ in 1:trials
        run_single_simulation(N, T_comp, tau, h, J, mode=:annealing)
    end
    plot!(p_comp, 1:every:T_comp, (res_anneal ./ trials)[1:every:T_comp], label="Annealing", lw=2)

    # 固定αのリスト
    for a in fixed_alphas
        println("Processing Fixed α = $a...")
        res_fixed = @distributed (+) for _ in 1:trials
            run_single_simulation(N, T_comp, tau, h, J, mode=:fixed, fixed_alpha=a)
        end
        plot!(p_comp, 1:every:T_comp, (res_fixed ./ trials)[1:every:T_comp], label="Fixed α=$a", lw=1.5, linestyle=:dash)
    end

    savefig(p_comp, "success_prob_detailed.png")
    println("Saved success_prob_detailed.png")
end

main_visualize()