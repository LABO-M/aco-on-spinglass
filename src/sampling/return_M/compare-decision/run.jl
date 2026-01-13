using Distributed
# 例: julia -p8 run.jl でワーカー作成済みなら addprocs は不要
# addprocs(8)

include(joinpath(@__DIR__, "simulation.jl"))
@everywhere include($(joinpath(@__DIR__, "simulation.jl")))

using .Simulation
using DelimitedFiles
using Printf

# ---- 設定（ここだけ変えればOK） ----
seed   = 42
N      = 100
tau    = 100
h      = 0.001
J      = 0.1
samples = 5000
burn_in = 1000

alpha_inc = 1e-6
alphas = [0.8, 0.835118240479, 0.9]   # αcはあなたの理論値に差し替え可

# 比較：linear と sigmoid
betas = [2.0, 5.0, 10.0]

outdir = joinpath(@__DIR__, "data")
mkpath(outdir)

function save_vec(path, v)
    open(path, "w") do io
        for x in v
            println(io, x)
        end
    end
end

# linear
for a in alphas
    M = Simulation.sample_M(N, 0.0, a, alpha_inc, tau, samples, h, J, seed;
        burn_in=burn_in, decision=:linear)
    fn = @sprintf("M_linear_tau%d_N%d_h%.4g_J%.4g_alpha%.12g_seed%d.csv", tau, N, h, J, a, seed)
    save_vec(joinpath(outdir, fn), M)
    @info "saved" fn
end

# sigmoid
for beta in betas
    for a in alphas
        M = Simulation.sample_M(N, 0.0, a, alpha_inc, tau, samples, h, J, seed;
            burn_in=burn_in, decision=:sigmoid, beta=beta)
        fn = @sprintf("M_sigmoid_beta%.3g_tau%d_N%d_h%.4g_J%.4g_alpha%.12g_seed%d.csv", beta, tau, N, h, J, a, seed)
        save_vec(joinpath(outdir, fn), M)
        @info "saved" fn
    end
end
