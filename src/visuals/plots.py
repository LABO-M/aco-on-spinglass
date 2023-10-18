import pandas as pd
import matplotlib.pyplot as plt
import os
import re

def plot_zt(csv_filepath):
    # CSVからデータを読み込む
    df = pd.read_csv(csv_filepath)

    # ファイル名からパラメータを取得
    params = extract_params_from_filename(os.path.basename(csv_filepath))

    # プロットのスタイルを設定
    linestyle_dict = {'1.0': ':', '2.0': '--'}
    linestyle = linestyle_dict[str(params['alpha'])]
    color = 'k' if params['tau'] == 'inf' else '0.5'

    # データのプロット
    plt.plot(df['t'], df['Z'], color=color, linestyle=linestyle, label=f'α={params["alpha"]}, t0={params["t0"]}, τ={params["tau"]}')

    # 軸とタイトルの設定
    plt.xlabel('t', fontsize=30)
    plt.ylabel('Z(t)', fontsize=30)
    plt.legend(fontsize=26)

    # y軸の指数表示をオフ
    plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter())
    # plt.xlim(1e3, 1e9)

    plt.grid(True, which="both", ls="--")
    plt.xscale('log')

    plt.xticks(fontsize=26)
    plt.yticks(fontsize=26, rotation=90, va='center')

    plt.show()
# 例：plot_zt('data/Zt/N100_T1M_t01K_alpha2.0_tau100.csv', 'plot/fig3/Zt_lineplot.eps')

def extract_params_from_filename(filename: str) -> dict:
    # 正規表現でファイル名からパラメータを抽出
    match = re.search(r'N(\d+)_T(\d+M?)_t0(\d+K?)_alpha([\d.]+)_tau(\d+|inf)', filename)
    if not match:
        raise ValueError(f"Invalid filename format: {filename}")

    # 抽出した各グループを適切な型に変換
    params = {
        'N': int(match.group(1)),
        'T': int(match.group(2).replace('M', '000000')),
        't0': int(match.group(3).replace('K', '000')),
        'alpha': float(match.group(4)),
        'tau': match.group(5)
    }
    return params
