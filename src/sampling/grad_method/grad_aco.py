import torch

# 変数の設定
def initialize_random_parameters(n, seed, device):
    torch.manual_seed(seed)
    m = torch.zeros(n, device=device)
    h = torch.ones(n, device=device) * 0.001
    J = torch.normal(mean=0.0, std=1.0, size=(n, n), device=device)
    J = torch.tril(J, diagonal=-1) + torch.tril(J, diagonal=-1).T
    initial_alpha = torch.tensor(0.0, device=device)
    return m, h, J, initial_alpha

# mの正負に応じてイジングモデルのエネルギーを計算
def energy(m, h, J):
    m = 2 * (m > 0.0) - 1
    return -torch.sum(h * m) - (torch.sum(J * torch.outer(m, m)) / (len(m) - 1))

# 自由エネルギー計算関数 (ベクトル化)
def free_energy(m, h, J, alpha):
    AS = -(1 - alpha) * torch.sum(torch.log(1 - m**2))
    return AS + alpha * energy(m, h, J)

# 勾配計算関数 (ベクトル化)
def update(m, h, J, alpha):
    interaction_grad = -2 * alpha * torch.sum(J * m, dim=1) / (len(m) - 1) 
    dif = interaction_grad + (1 - alpha) * m / (1 - m**2) - alpha * h
    return dif

# 勾配降下法（α-annealing）
def gradient_descent_annealing(m, h, J, alpha, alpha_inc, iter, lr, tol):
    alpha_increment_list = []
    energy_series = []
    for iteration in range(iter):
        grad = update(m, h, J, alpha)
        m = m - lr * grad  # 勾配降下ステップ

        # 停滞のチェック（勾配が小さくなる場合）
        if torch.max(abs(grad)) < tol:
            alpha = min(alpha + alpha_inc, 0.9999)  # alpha を増加
            alpha_increment_list.append(1)
        else:
            alpha_increment_list.append(0)
        energy_series.append(energy(m, h, J).cpu().item())

    return energy_series, alpha_increment_list

# 勾配降下法（fixed α）
def gradient_descent(m, h, J, alpha, iter, lr, tol):
    energy_series = []
    for iteration in range(iter):
        grad = update(m, h, J, alpha)
        m = m - lr * grad  # 勾配降下ステップ
        energy_series.append(energy(m, h, J).cpu().item())

        # 停滞のチェック（勾配が小さくなる場合）
        if torch.max(abs(grad)) < tol:
            break
    return energy_series
