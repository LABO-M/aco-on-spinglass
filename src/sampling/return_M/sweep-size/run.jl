using Distributed
include(joinpath(@__DIR__, "simulation.jl"))
@everywhere include($(joinpath(@__DIR__, "simulation.jl")))

using .Simulation
using Printf

seed   = 42
tau    = 100
h      = 0.001
J      = 0.1
samples = 5000
burn_in = 1000

alpha_inc = 1e-6
alpha = 0.835118240479     # αc近傍など
Ns = [50, 100, 200, 500, 1000]

outdir = joinpath(@__DIR__, "data")
mkpath(outdir)

function save_vec(path, v)
    open(path, "w") do io
        for x in v
            println(io, x)
        end
    end
end

for N in Ns
    M = Simulation.sample_M(N, 0.0, alpha, alpha_inc, tau, samples, h, J, seed;
        burn_in=burn_in, decision=:linear)
    fn = @sprintf("M_size_tau%d_N%d_h%.4g_J%.4g_alpha%.12g_seed%d.csv", tau, N, h, J, alpha, seed)
    save_vec(joinpath(outdir, fn), M)
    @info "saved" fn
end
