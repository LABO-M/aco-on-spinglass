import torch
import numpy as np
import output
import importlib
importlib.reload(output)

# mの正負に応じてイジングモデルのエネルギーを計算
# Not devided by N-1
def energy(m, h, J):
    m = 2 * (m > 0.0) - 1
    return -torch.sum(h * m) - (torch.sum(J * torch.outer(m, m)) - torch.sum(torch.diagonal(J * torch.outer(m, m))))

# エネルギー計算関数 (ベクトル化)
# Not devided by N-1
def f(m, h, J, tau, alpha):
    AS = -(1 - alpha) * torch.sum(torch.log(1 - m**2))
    interaction = (torch.sum(J * torch.outer(m, m)) - torch.sum(torch.diagonal(J * torch.outer(m, m)))) 
    ising = -alpha * (torch.sum(h * m) + interaction)
    return AS + ising

# 勾配計算関数 (ベクトル化)
# Not devided by N-1
def update(m, h, J, tau, alpha):
    interaction_grad = (-2 * alpha * (torch.sum(J * m, dim=1) - torch.diagonal(J) * m)) / (len(m) - 1) 
    dif = interaction_grad + (1 - alpha) * m / (1 - m**2) - alpha * h
    return dif

# 勾配降下法 (改良版)
def gradient_descent(m, h, J, tau, alpha, alpha_inc, max_iter, lr, tol):
    best_energy = 0
    alpha_increment_list = []
    energy_series = []
    for iteration in range(max_iter):
        grad = update(m, h, J, tau, alpha)
        m = m - lr * grad  # 勾配降下ステップ
        if best_energy > energy(m, h, J):
            best_energy = energy(m, h, J)

        # 停滞のチェック（勾配が小さくなる場合）
        if torch.max(abs(grad)) < tol:
            alpha = min(alpha + alpha_inc, 0.9999)  # alpha を増加
            alpha_increment_list.append(1)
        else:
            alpha_increment_list.append(0)
        ## alpha が 0.99 に達した場合、終了
        #if alpha >= 0.999:
        #    break
        energy_series.append(energy(m, h, J).cpu().item())
    #print(type(energy(m, h, J).cpu().item()))
    #print(len(energy_series))

    return energy_series, best_energy, alpha_increment_list

# シード値の設定
def initialize_random_parameters(n, seed):
    torch.manual_seed(seed)
    J = torch.normal(mean=0.0, std=1.0, size=(n, n))
    J = torch.tril(J) + torch.tril(J).T
    return J