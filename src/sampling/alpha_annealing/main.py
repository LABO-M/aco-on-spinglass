import torch
import simulation
import output
import importlib
importlib.reload(simulation)
importlib.reload(output)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
n = 100
start_seed = 42
seed = 42
iter = 1e6
energy_series = []
z_series = []
z_series = []

def m_spinglass(n, start_seed, iter, energy_series, m_series):
    seed = start_seed
    spins, h, J, st, st1, alpha, alpha_inc = simulation.initialize_spinglass(n, seed, iter, device)

    for i in range(int(iter)):
        torch.manual_seed(seed)
        spins = 2 * (simulation.decision_probabilities(st, st1, alpha) > torch.rand(n, device=device)) - 1
        st += torch.exp(simulation.calculate_energy(spins, h, J))
        st1 += torch.exp(simulation.calculate_energy(spins, h, J)) * ((spins + 1) / 2)
        z = (st1 / st)
        alpha += alpha_inc
        m = 2 * alpha * (z - 0.5)
        energy_series.append(simulation.calculate_energy(spins, h, J).cpu().item())
        m_series.append(m.tolist())
        seed += 1

    energy_series.insert(0, 'energy')
    energy_file_name = 'tauinf_energy'
    m_file_name = 'tauinf_m'
    output.list_to_csv(energy_file_name, energy_series, start_seed, iter)
    output.matrix_to_csv(m_file_name, m_series, start_seed, iter)


def z_spinglass(n, start_seed, iter, energy_series, z_series):
    seed = start_seed
    spins, h, J, st, st1, alpha, alpha_inc = simulation.initialize(n, seed, iter, device)

    for i in range(int(iter)):
        torch.manual_seed(seed)
        spins = 2 * (simulation.decision_probabilities(st, st1, alpha) > torch.rand(n, device=device)) - 1
        st += torch.exp(simulation.calculate_energy(spins, h, J))
        st1 += torch.exp(simulation.calculate_energy(spins, h, J)) * ((spins + 1) / 2)
        z = (st1 / st)
        alpha += alpha_inc
        energy_series.append(simulation.calculate_energy(spins, h, J).cpu().item())
        z_series.append(z.tolist())
        seed += 1

    energy_series.insert(0, 'energy')
    energy_file_name = 'tauinf_energy'
    m_file_name = 'tauinf_m'
    output.list_to_csv(energy_file_name, energy_series, start_seed, iter)
    output.matrix_to_csv(m_file_name, z_series, start_seed, iter)


def m_homogenerous(n, start_seed, iter, energy_series, m_series):
    seed = start_seed
    spins, h, J, st, st1, alpha, alpha_inc = simulation.initialize(n, seed, iter, device)

m_spinglass