import torch

def initialize_random_parameters(n, seed):
    torch.manual_seed(seed)
    x = torch.tensor([1 if torch.rand(1).item() < 0.5 else -1 for _ in range(n)])
    J = torch.normal(mean=0.0, std=1.0, size=(n, n))
    J = torch.tril(J, diagonal=-1) + torch.tril(J, diagonal=-1).T  # 下三角のみ
    h = torch.zeros(n)  # 外部磁場
    return x, J, h

def calculate_energy(x, h, J):
    interaction = -torch.sum(J * torch.outer(x, x)) / (len(x) - 1)  # 相互作用エネルギー
    magnetic_field = -torch.sum(h * x)  # 磁場エネルギー
    return interaction + magnetic_field

def delta_energy(x, flipped_index, h, J):
    x_i = x[flipped_index]
    local_field = h[flipped_index] + 2 * torch.sum(J[flipped_index, :] * x)
    return 2 * x_i * local_field

def simulation(n, initial_temperature, goal_temperature, seed, iter):
    x, J, h = initialize_random_parameters(n, seed)
    T = initial_temperature
    cooling_rate = (goal_temperature / initial_temperature) ** (1 / iter)
    energy_series = []

    for _ in range(iter):
        random_index = torch.randint(0, n, (1,)).item()
        delta = delta_energy(x, random_index, h, J)
        
        if delta < 0 or torch.rand(1).item() < torch.exp(-delta / T):
            x[random_index] *= -1
        
        T *= cooling_rate
        energy_series.append(calculate_energy(x, h, J).item())

    return energy_series