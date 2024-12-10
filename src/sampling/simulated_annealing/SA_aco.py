#simulated annealing
import torch
import output

# set J
def initialize_random_parameters(n, seed1):
    torch.manual_seed(seed1)
    x = torch.tensor([1 if torch.rand(1).item() < 0.5 else -1 for _ in range(n)])
    J = torch.normal(mean=0.0, std=1.0, size=(n, n))
    J = torch.tril(J) + torch.tril(J).T
    return x, J

# caluculate energy
def calculate_energy(x, h, J):
    magnetic_field = torch.sum(h * x)
    interaction = torch.sum(J * torch.outer(x, x)) - torch.sum(torch.diagonal(J * torch.outer(x, x)))
    
    return -(magnetic_field + interaction)

def delta_energy(current_spin, temp_spin, h, J):
    return calculate_energy(temp_spin, h, J) - calculate_energy(current_spin, h, J)

def calculate_cooling_rate(initial_temperature, goal_temperature, iter):
    return (goal_temperature / initial_temperature) ** (1 / iter)

def accept_prob(current_spin, temp_spin, h, J, T):
    delta = delta_energy(current_spin, temp_spin, h, J)

    if delta < 0:
        accept_prob = 1
    else:
        accept_prob = torch.exp(- delta / T)

    return accept_prob

def judge(current_spin, temp_spin, h, J, T, seed2):
    torch.manual_seed(seed2)
    
    if torch.rand(1) < accept_prob(current_spin, temp_spin, h, J, T):
        return True
    else:
        return False

def simulation(n, h, initial_temperature, goal_temperature, seed1, seed2, iter):
    
    x, J = initialize_random_parameters(n, seed1)
    T = initial_temperature
    cooling_rate = calculate_cooling_rate(initial_temperature, goal_temperature, iter)
    energy_series = []
    
    for i in range(iter):
        current_spin = x.clone() # current spin

        print(calculate_energy(x, h, J))

        random_index = torch.randint(0, len(x), (1,)).item() # create random index
        x[random_index] *= -1

        if judge(current_spin, x, h, J, T, seed2) == False:
            x = current_spin
        
        T *= cooling_rate
        energy_series.append(calculate_energy(x, h, J))

    return energy_series