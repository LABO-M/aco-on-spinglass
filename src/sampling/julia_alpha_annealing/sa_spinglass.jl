using Distributions, Random, LinearAlgebra

function initilize(n, seed)
    Random.seed!(seed)
    spins = 2 .* (rand(n) .> 0.5) .- 1
    spins = reshape(spins, :, 1)
    h = 0.001
    J = randn(n, n)
    J = Symmetric(J - Diagonal(J))
    return spins, h, J
end

function calculate_energy(spins, h, J)
    return -(sum(h .* spins) + sum(J .* (spins * spins')) ./ (length(spins) - 1))
end

function delta_energy(spins, flip_index, h, J)
    flip_spin = spins[flip_index]
    local_field = h + 2 * sum(J[flip_index, :] .* spins)
    return 2 * flip_spin * local_field
end

n = 100
seed = 43
initial_temperature = 100
goal_temperature = 0.01
iter = 1e6

spins, h, J = initilize(n, seed)
T = initial_temperature
cooling_rate = (goal_temperature / initial_temperature) ^ (1 / iter)
energy_series = []

for i in 1:iter
    flip_index = rand(1:n)
    delta = delta_energy(spins, flip_index, h, J)

    if delta < 0 || rand() < exp(-delta / T)
        spins[flip_index] = -spins[flip_index]
    end

    global T *= cooling_rate
    push!(energy_series, calculate_energy(spins, h, J))
end
println(minimum(energy_series))