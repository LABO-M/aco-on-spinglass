using Distributed
include(joinpath(@__DIR__, "simulation.jl"))
@everywhere include($(joinpath(@__DIR__, "simulation.jl")))

using .Simulation
using Printf

seed   = 42
N      = 100
tau    = 1000
h      = 0.001
J      = 0.1
samples = 5000
burn_in = 1000

alpha = 0.835118240479

# linear
alpha_inc1 = 1e-6

# two-stage
alpha_mid = 0.75
alpha_inc2 = 1e-7

outdir = joinpath(@__DIR__, "data")
mkpath(outdir)

function save_vec(path, v)
    open(path, "w") do io
        for x in v
            println(io, x)
        end
    end
end

# schedule=linear
M = Simulation.sample_M(N, 0.0, alpha, alpha_inc1, tau, samples, h, J, seed;
    burn_in=burn_in, decision=:linear, schedule=:linear)
fn = @sprintf("M_anneal_linear_tau%d_N%d_h%.4g_J%.4g_alpha%.12g_seed%d.csv", tau, N, h, J, alpha, seed)
save_vec(joinpath(outdir, fn), M)
@info "saved" fn

# schedule=two_stage
M = Simulation.sample_M(N, 0.0, alpha, alpha_inc1, tau, samples, h, J, seed;
    burn_in=burn_in, decision=:linear,
    schedule=:two_stage, alpha_mid=alpha_mid, alpha_inc2=alpha_inc2)
fn = @sprintf("M_anneal_two_stage_mid%.3g_inc2%.3g_tau%d_N%d_h%.4g_J%.4g_alpha%.12g_seed%d.csv",
              alpha_mid, alpha_inc2, tau, N, h, J, alpha, seed)
save_vec(joinpath(outdir, fn), M)
@info "saved" fn
