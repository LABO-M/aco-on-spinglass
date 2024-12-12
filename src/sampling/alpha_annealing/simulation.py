import torch

def initialize(n, seed, device):
    torch.manual_seed(seed)
    spins = torch.tensor([1 if torch.rand(1).item() < 0.5 else -1 for _ in range(n)], device=device)
    J = torch.normal(mean=0.0, std=1.0, size=(n, n), device=device)
    J = torch.tril(J, diagonal=-1) + torch.tril(J, diagonal=-1).T  # 下三角のみ
    h = torch.zeros(n, device=device)  # 外部磁場
    st = torch.ones(n, device=device) * torch.exp(calculate_energy(spins, h, J))
    st1 = st * ((spins + 1) / 2)
    return spins, h, J, st, st1

def decision_probabilities(st, st1, alpha):
    return alpha * ((st1 / st) - 0.5) + 0.5

def calculate_energy(spins, h, J):
    return -(torch.sum(h * spins) + (torch.sum(J * torch.outer(spins, spins)) / (len(spins) - 1)))

