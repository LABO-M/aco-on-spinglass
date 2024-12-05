import csv

def list_to_csv(name, data):
    with open(f'/home/mori-lab/shimizu/ACO/ACO-on-SpinGlass/src/sampling/grad_method/data/{name}.csv', 'w', newline='') as f:
        for _ in data:
            print(_, sep=',', file=f)

def matrix_to_csv(name, data):
    with open(f'/home/mori-lab/shimizu/ACO/ACO-on-SpinGlass/src/sampling/grad_method/data/{name}.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(data)
