import os
import csv

def ensure_dir_exists(file_path):
    # ディレクトリ部分を取得して、存在しない場合は作成
    dir_name = os.path.dirname(file_path)
    if dir_name and not os.path.exists(dir_name):
        os.makedirs(dir_name)

def list_to_csv(name, data):
    file_path = f'/home/mori-lab/shimizu/ACO/ACO-on-SpinGlass/src/sampling/simulated_annealing/data/{name}.csv'
    ensure_dir_exists(file_path)  # ディレクトリの存在確認
    with open(file_path, 'w', newline='') as f:
        for _ in data:
            print(_, sep=',', file=f)

def matrix_to_csv(name, data):
    file_path = f'/home/mori-lab/shimizu/ACO/ACO-on-SpinGlass/src/sampling/simulated_annealing/data/{name}.csv'
    ensure_dir_exists(file_path)  # ディレクトリの存在確認
    with open(file_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(data)
