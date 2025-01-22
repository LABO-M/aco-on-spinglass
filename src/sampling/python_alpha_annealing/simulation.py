import torch

def initialize_spinglass(n, seed, iter, start_alpha, device):
    torch.manual_seed(seed)
    spins = torch.tensor([1 if torch.rand(1).item() < 0.5 else -1 for _ in range(n)], device=device)
    J = torch.normal(mean=0.0, std=1.0, size=(n, n), device=device)
    J = torch.tril(J, diagonal=-1) + torch.tril(J, diagonal=-1).T  # 下三角のみ
    h = torch.ones(n, device=device) * 0.001 # 外部磁場
    st1 = torch.exp(calculate_energy(spins, h, J)) * ((spins + 1) / 2)
    st0 = torch.exp(calculate_energy(spins, h, J)) * ((- spins + 1) / 2)
    alpha = torch.tensor(start_alpha, device=device)
    alpha_inc = torch.tensor((1 - start_alpha)/iter, device=device)

    return h, J, st1, st0, alpha, alpha_inc

#def initialize_homogeneous(n, seed, iter, start_alpha, device):
#    torch.manual_seed(seed)
#    spins = torch.tensor([1 if torch.rand(1).item() < 0.5 else -1 for _ in range(n)], device=device)
#    h = torch.ones(n, device=device) * 0.001 # 外部磁場
#    J = (torch.ones(size=(n, n), device=device) - torch.eye(n, device=device)) * 0.1 # 交互作用
#    st1 = torch.exp(calculate_energy(spins, h, J)) * ((spins + 1) / 2)
#    st0 = torch.exp(calculate_energy(spins, h, J)) * ((- spins + 1) / 2)
#    alpha_inc = torch.tensor((1 - start_alpha)/iter, device=device)
#    alpha = torch.tensor(start_alpha, device=device)
#
#    return spins, h, J, st1, st0, alpha, alpha_inc

def calculate_external_effect(determined_spins, J):
    result = 0
    for i in range(len(determined_spins)):
        result -= determined_spins[i] * J[i][len(determined_spins)]
    
    return result

def decision_probabilities(st1, st0, index, alpha, beta, external_effect):
    return (1 - alpha) * 0.5 + alpha * (st1[index] * torch.exp(-beta * external_effect) / (st1[index] * torch.exp(-beta * external_effect) + st0[index] * torch.exp(beta * external_effect)))

def calculate_energy(spins, h, J):
    return -(torch.sum(h * spins) + (torch.sum(J * torch.outer(spins, spins)) / (len(spins) - 1)))

