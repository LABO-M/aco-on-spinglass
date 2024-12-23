using Distributions, Random, LinearAlgebra

function initilize(n, seed)
    Random.seed!(seed)
    spins = 2 .* (rand(n) .> 0.5) .- 1
    spins = reshape(spins, :, 1)
    h = 0.001
    J = fill(0.1, n, n)
    J = J - Diagonal(J)
    st = exp(-calculate_energy(spins, h, J))
    st1 = st .* ((spins .+ 1) ./ 2)
    return spins, h, J, st, st1
end

function calculate_energy(spins, h, J)
    return -(sum(h .* spins) + sum(J .* (spins * spins')) ./ (length(spins) - 1))
end

function decision_probability(alpha, st, st1)
    return alpha .* ((st1 ./ st) .- 0.5) .+ 0.5
end

function calculate_m(alpha, st, st1)
    return 2 .* alpha .* ((st1 ./ st) .- 0.5)
end


n = 100
seed = 42
iter = 1e6
alpha = 0.0
tau = 1000
evapoate = exp(-1/tau)
energy = []
m = []
z = []
spins, h, J, st, st1 = initilize(n, seed)

for i in 1:iter
    global spins = 2 .* (decision_probability(alpha, st, st1) .> rand(n)) .- 1
    global st = st * evapoate + exp(-calculate_energy(spins, h, J))
    global st1 = st1 .* evapoate + exp(-calculate_energy(spins, h, J)) .* ((spins .+ 1) ./ 2)
    push!(energy, calculate_energy(spins, h, J))
    push!(m, sum(calculate_m(alpha, st, st1) .> 0))
    global alpha = min(alpha + 1/iter, 0.999)
end
println(energy[end])
println(m[end])
