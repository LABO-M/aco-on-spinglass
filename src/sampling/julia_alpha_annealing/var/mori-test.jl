using Random
using Printf
using CSV
using DataFrames
# ユーティリティ関数
clamp_val(x, lo, hi) = x < lo ? lo : (x > hi ? hi : x)
# P(X=1 | h_eff, z, alpha, beta) = e^{βh} z^α / ( e^{βh} z^α + e^{-βh} (1-z)^α )
function prob_one(beta, h_eff, z, alpha)
    z = clamp_val(z, 1e-12, 1.0 - 1e-12)
    if alpha <= 0.0
        # alpha = 0 の極限: z^0 = (1-z)^0 = 1 → e^{βh} / (e^{βh} + e^{-βh})
        eh = exp(beta * h_eff)
        emh = 1.0 / eh
        return eh / (eh + emh)
    end
    if beta == 0.0
        za = z^alpha
        oma = (1.0 - z)^alpha
        return za / (za + oma)
    else
        eh = exp(beta * h_eff)
        emh = 1.0 / eh
        za = z^alpha
        oma = (1.0 - z)^alpha
        num = eh * za
        den = num + emh * oma
        return num / den
    end
end
# パラメータ構造体
mutable struct Params
    N::Int                  # スピン数
    J::Float64             # 交換相互作用（同質）
    h::Float64             # 外場（同質）
    alpha::Float64         # Final alpha
    beta::Float64          # 固定 β
    tau::Float64           # 蒸発タイムスケール
    steps::Int             # 各αのMCS数（既定 2000）
    blocks::Int            # αのブロック数（既定 100 = 0.00..0.99）
    output_skip::Int       # 何回に1回出力するか
    out_csv::String        # 出力CSV
end
function parse_args(args::Vector{String})
    # 既定値
    P = Params(
        1024,              # N
        0.01,              # J
        0.00,              # h
        1.0,               # alpha
        0.0,               # beta
        100.0,             # tau
        2000,              # steps
        10000,             # blocks
        10,                # output_skip
        "traj.csv"         # out_csv
    )
    i = 1
    while i <= length(args)
        if args[i] == "--N" && i + 1 <= length(args)
            P.N = parse(Int, args[i+1])
            i += 2
        elseif args[i] == "--J" && i + 1 <= length(args)
            P.J = parse(Float64, args[i+1])
            i += 2
        elseif args[i] == "--h" && i + 1 <= length(args)
            P.h = parse(Float64, args[i+1])
            i += 2
        elseif args[i] == "--alpha" && i + 1 <= length(args)
            P.alpha = parse(Float64, args[i+1])
            i += 2
        elseif args[i] == "--beta" && i + 1 <= length(args)
            P.beta = parse(Float64, args[i+1])
            i += 2
        elseif args[i] == "--tau" && i + 1 <= length(args)
            P.tau = parse(Float64, args[i+1])
            i += 2
        elseif args[i] == "--steps" && i + 1 <= length(args)
            P.steps = parse(Int, args[i+1])
            i += 2
        elseif args[i] == "--blocks" && i + 1 <= length(args)
            P.blocks = parse(Int, args[i+1])
            i += 2
        elseif args[i] == "--skip" && i + 1 <= length(args)
            P.output_skip = parse(Int, args[i+1])
            i += 2
        elseif args[i] == "--out" && i + 1 <= length(args)
            P.out_csv = args[i+1]
            i += 2
        elseif args[i] == "-h" || args[i] == "--help"
            println("""
            Usage: julia aca_beta_fixed.jl [--N 1024] [--J 0.01] [--h 0.001] [--alpha 1.0] [--beta 0.0]
                      [--tau 100] [--steps 2000] [--blocks 100] [--out traj.csv]
            Note: α runs 0.0000..0.9999 by default with $(P.blocks) blocks.
            """)
            exit(0)
        else
            i += 1
        end
    end
    return P
end
function main()
    # コマンドライン引数のパース
    P = parse_args(ARGS)
    N = P.N
    J = P.J
    h = P.h
    alpha_final = P.alpha
    beta = P.beta
    tau = P.tau
    STEP = P.steps
    OSkip = P.output_skip
    B = P.blocks
    Ttot = STEP * B
    # 状態配列の初期化
    X = zeros(Bool, N)      # スピン状態 (0/1 → false/true)
    Z = fill(0.5, N)        # フェロモン濃度
    S1 = zeros(Float64, N)  # フェロモン累積
    S = 0.0  # 全体の重み累積
    # 物理量
    decay = exp(-1.0 / tau)
    E_min = -N * abs(h) - 0.5 * N * J
    # 結果出力
    println("N=$N,J=$J,h=$h,alpha=$alpha_final,beta=$beta,tau=$tau,STEP=$STEP,blocks=$B")
    # CSVファイル準備
    open(P.out_csv, "w") do fp
        println(fp, "# N=$N,J=$J,h=$h,alpha=$alpha_final,beta=$beta,tau=$tau,SPA=$STEP,blocks=$B")
        println(fp, "t,alpha,beta,block,half_flag,m,E,Z_mean,EXN")
        t_global = 0
        da = alpha_final / B
        for ib in 0:(B-1)
            current_alpha = da * ib
            for t in 0:(STEP-1)
                t_global += 1
                half_flag = (t < (STEP ÷ 2)) ? 1 : 0  # 前半SPA/2→1, 後半SPA/2→0
                sum_sigma = 0
                # 逐次 i = 1..N (Juliaは1始まり)
                for i in 1:N
                    if i == 1
                        X[1] = (half_flag == 1)  # X(1,t) を強制
                        sum_sigma += X[1] ? 1 : -1
                        continue
                    end
                    # h_eff(i,t) = h + (J / (i-1)) * sum_{j<i} sigma_j
                    h_eff = h + J * sum_sigma / (i - 1)
                    p1 = prob_one(beta, h_eff, Z[i], current_alpha)
                    u = rand()  # [0,1) の一様乱数
                    X[i] = (u < p1)
                    sum_sigma += X[i] ? 1 : -1
                end
                # m(t), E(t)
                m = sum_sigma / N
                E = -N * h * m - (J / (2.0 * (N - 1.0))) * (N * N * m * m - N)
                # m(N-1,t) : 最後のサイトを除いた平均
                m_N1 = (sum_sigma - (X[N] ? 1 : -1)) / (N - 1)
                h_eff_final = h + J * m_N1
                EXN = prob_one(beta, h_eff_final, Z[N], current_alpha)
                # フェロモン更新
                w = exp(-(E - E_min))
                S = S * decay + w
                for i in 1:N
                    S1[i] = S1[i] * decay + w * (X[i] ? 1.0 : 0.0)
                    Z[i] = S1[i] / S
                end
                # Z_mean計算
                Z_mean = sum(Z) / N
                # 出力
                if ib % OSkip == 0
                    @printf(fp, "%d,%.4f,%.4f,%d,%d,%.10g,%.10g,%.10g,%.10g\n",
                            t_global, current_alpha, beta, ib, half_flag, m, E, Z_mean, EXN)
                end
            end
        end
    end
    println("Done. wrote $Ttot rows to $(P.out_csv)")
    return nothing
end
# スクリプトとして実行された場合にmain()を呼び出す
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end