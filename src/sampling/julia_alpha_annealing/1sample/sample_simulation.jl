module simulation

using Random
using CSV
using DataFrames
using LinearAlgebra
using Statistics
using Printf
using FilePathsBase
using Base.Threads  # ← 追加：スレッド並列

# --- 共通関数群 ---

# エネルギーは i<j のみ集計（重複回避）
# 呼び出し側から triu(J,1) を渡す
function calculate_energy(spins::AbstractVector, h, J_upper::AbstractMatrix)
    n = length(spins)
    hvec = isa(h, Number) ? fill(h, n) : h
    S = spins * transpose(spins)         # 外積（n×n）
    pair_sum = sum(J_upper .* S)         # i<j のみ
    return -sum(hvec .* spins) - pair_sum / max(n - 1, 1)
end

# 既定：これまでに決まったスピンとの平均相互作用（列 new_spin_index）
function calculate_external_effect(determined_spins::AbstractVector, J::AbstractMatrix, new_spin_index::Int)
    return -mean(determined_spins .* J[1:length(determined_spins), new_spin_index])
end

function decision_probabilities(st1::AbstractVector, st0::AbstractVector, index::Int, alpha::Float64, beta::Float64, external_effect::Float64)
    z = st1[index] / (st1[index] + st0[index])
    num = (z^alpha) * exp(-beta * external_effect)
    den = num + ((1 - z)^alpha) * exp(beta * external_effect)
    return num / den
end

# 安全に成分 z[i] を取り出す（範囲外は missing）
@inline function get_component(zvec::AbstractVector{<:Real}, i::Int)
    return (1 <= i <= length(zvec)) ? zvec[i] : missing
end

# =========================================================
#  - フェーズAは廃止
#  - α を 0.0→alpha_max まで alpha_step=0.001 刻みで走らせ、
#    各 α で iters_per_alpha 回、毎回:
#      * スピン1も含め 1..n を順に決定（最初のスピンはフェロモン参照）
#      * st1, st0 を更新
#      * zvec の成分 z[2], z[3], z[100], z[end], z_mean を蓄積
#    最後の1000イテレーションを「1行=1イテレーション」で保存用に返す
# =========================================================
function simulate(power::Int, tau::Float64, beta::Float64, seed::Int;
                  magnetic::Float64=0.0, interaction::Float64=0.0,
                  alpha_step::Float64=0.001, iters_per_alpha::Int=10000, alpha_max::Float64=1.0)

    rng = MersenneTwister(seed)

    n = 2 ^ power
    J = ones(n, n) * interaction
    J[diagind(J)] .= 0
    J_upper = triu(J, 1)                 # i<j のみ
    h = magnetic == 0.0 ? fill(0.0, n) : fill(magnetic, n)

    # フェロモン（記憶）
    st1 = fill(1.0, n)
    st0 = fill(1.0, n)

    # 返却用（αごとに最後の1000を窓として保持）
    z2_windows    = Vector{Vector{Union{Missing,Float64}}}()
    z3_windows    = Vector{Vector{Union{Missing,Float64}}}()
    z100_windows  = Vector{Vector{Union{Missing,Float64}}}()
    zend_windows  = Vector{Vector{Union{Missing,Float64}}}()
    zmean_windows = Vector{Vector{Float64}}()
    alpha_series  = Float64[]

    alpha = 0.0
    while alpha <= alpha_max + 1e-12
        # 各αで蓄積（全 iters_per_alpha 分）
        z2_list    = Union{Missing,Float64}[]
        z3_list    = Union{Missing,Float64}[]
        z100_list  = Union{Missing,Float64}[]
        zend_list  = Union{Missing,Float64}[]
        zmean_list = Float64[]

        for _ in 1:iters_per_alpha
            determined_spins = Int[]
            # 1..n の順に決定（最初のスピンもフェロモン参照）
            for new_index in 1:n
                if new_index == 1
                    # 最初のスピン：相互作用相手がいないので外場のみ
                    external_effect = -magnetic
                else
                    external_effect = -magnetic + calculate_external_effect(determined_spins, J, new_index)
                end
                prob = decision_probabilities(st1, st0, new_index, alpha, beta, external_effect)
                new_spin = 2 * (rand(rng) < prob) - 1
                push!(determined_spins, new_spin)
            end

            # フェロモン更新
            energy = calculate_energy(determined_spins, h, J_upper)
            decay = exp(-1 / tau)
            weight = exp(-energy)
            st1 .= st1 .* decay .+ weight .* ((determined_spins .+ 1) ./ 2)
            st0 .= st0 .* decay .+ weight .* ((-determined_spins .+ 1) ./ 2)

            # 観測
            zvec = st1 ./ (st1 .+ st0)
            push!(z2_list,    get_component(zvec, 2))
            push!(z3_list,    get_component(zvec, 3))
            push!(z100_list,  get_component(zvec, 100))
            push!(zend_list,  get_component(zvec, n))
            push!(zmean_list, mean(zvec))
        end

        # 最後の1000イテレーションだけ切り出して保存
        w = min(length(zmean_list), 10000)
        range_idx = (length(zmean_list) - w + 1):length(zmean_list)

        push!(z2_windows,    collect(z2_list[range_idx]))
        push!(z3_windows,    collect(z3_list[range_idx]))
        push!(z100_windows,  collect(z100_list[range_idx]))
        push!(zend_windows,  collect(zend_list[range_idx]))
        push!(zmean_windows, collect(zmean_list[range_idx]))
        push!(alpha_series,  alpha)

        # αを進める
        alpha = round(alpha + alpha_step, digits=6)
        if alpha > alpha_max
            break
        end
    end

    return z2_windows, z3_windows, z100_windows, zend_windows, zmean_windows, alpha_series
end


function sampling(power::Int, tau::Float64, beta::Float64, seed::Int,
                  magnetic::Float64, interaction::Float64;
                  alpha_step::Float64=0.001, iters_per_alpha::Int=1000, alpha_max::Float64=1.0)

    z2_windows, z3_windows, z100_windows, zend_windows, zmean_windows, alpha_series =
        simulate(power, tau, beta, seed;
                 magnetic=magnetic, interaction=interaction,
                 alpha_step=alpha_step, iters_per_alpha=iters_per_alpha, alpha_max=alpha_max)

    # αごとに最後の1000行（または iters_per_alpha が小さければその数）を展開
    rows = Vector{NamedTuple{(:alpha,:iter_in_alpha,:z2,:z3,:z100,:zend,:z_mean),
                             Tuple{Float64,Int,Union{Missing,Float64},Union{Missing,Float64},Union{Missing,Float64},Union{Missing,Float64},Float64}}}()

    for (k, α) in enumerate(alpha_series)
        w = length(zmean_windows[k])
        for t in 1:w
            push!(rows, (
                alpha        = α,
                iter_in_alpha= t,
                z2           = z2_windows[k][t],
                z3           = z3_windows[k][t],
                z100         = z100_windows[k][t],
                zend         = zend_windows[k][t],
                z_mean       = zmean_windows[k][t],
            ))
        end
    end

    df = DataFrame(rows)

    dir_switch = magnetic == 0.0 ? "symmetric" : "asymmetric"
    dir_path = "/home/mori-lab/shimizu/aco/data/ising/annealing/onesample/$(dir_switch)"
    mkpath(dir_path)

    filename_all = @sprintf("beta%.1e_n2^%d_tau%.1e_results.csv", beta, power, tau)
    full_path_all = joinpath(dir_path, filename_all)
    CSV.write(full_path_all, df)

    return filename_all
end

# =========================================================
# 追加関数：複数サンプルをスレッド並列で実行し、1つのCSVに統合
# =========================================================
function sampling_parallel(power::Int, tau::Float64, beta::Float64, start_seed::Int,
                           magnetic::Float64, interaction::Float64;
                           alpha_step::Float64=0.001, iters_per_alpha::Int=1000, alpha_max::Float64=1.0,
                           samples::Int=10, seed_stride::Int=1000)

    # 各サンプルの行（NamedTuple配列）を保持
    rows_per_sample = Vector{Vector{NamedTuple}}(undef, samples)

    @threads for s in 1:samples
        seed = start_seed + (s-1) * seed_stride

        z2_windows, z3_windows, z100_windows, zend_windows, zmean_windows, alpha_series =
            simulate(power, tau, beta, seed;
                     magnetic=magnetic, interaction=interaction,
                     alpha_step=alpha_step, iters_per_alpha=iters_per_alpha, alpha_max=alpha_max)

        # αごとに最後のウィンドウを展開し、sample列を付与
        local_rows = Vector{NamedTuple{(:sample,:alpha,:iter_in_alpha,:z2,:z3,:z100,:zend,:z_mean),
                          Tuple{Int,Float64,Int,Union{Missing,Float64},Union{Missing,Float64},Union{Missing,Float64},Union{Missing,Float64},Float64}}}()

        for (k, α) in enumerate(alpha_series)
            w = length(zmean_windows[k])
            for t in 1:w
                push!(local_rows, (
                    sample       = s,
                    alpha        = α,
                    iter_in_alpha= t,
                    z2           = z2_windows[k][t],
                    z3           = z3_windows[k][t],
                    z100         = z100_windows[k][t],
                    zend         = zend_windows[k][t],
                    z_mean       = zmean_windows[k][t],
                ))
            end
        end

        rows_per_sample[s] = local_rows
    end

    # 結合して DataFrame 化
    all_rows = reduce(vcat, rows_per_sample)
    df = DataFrame(all_rows)

    # 出力
    dir_switch = magnetic == 0.0 ? "symmetric" : "asymmetric"
    dir_path = "/home/mori-lab/shimizu/aco/data/ising/annealing/onesample/$(dir_switch)"
    mkpath(dir_path)

    filename_all = @sprintf("beta%.1e_sample%d_n2^%d_tau%.1e_results.csv", beta, samples, power, tau)
    full_path_all = joinpath(dir_path, filename_all)
    CSV.write(full_path_all, df)

    return filename_all
end

end
