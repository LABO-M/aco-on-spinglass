import pandas as pd
import matplotlib.pyplot as plt
import os
import re

def plot_zt(csv_filepaths):
    print(f"csv_filepaths: {csv_filepaths}")
    for csv_filepath in csv_filepaths:
        print(f"Processing: {csv_filepath}")
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

    plt.grid(True, which="both", ls="--")
    plt.xscale('log')

    plt.xticks(fontsize=26)
    plt.yticks(fontsize=26, rotation=90, va='center')

    plt.show()

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

def generate_filename_from_params(param):
    """
    パラメータの辞書を受け取り、対応するファイル名を生成する関数
    """
    tau_str = "inf" if param["tau"] == -1 else str(param["tau"])
    return f"data/Zt/N{param['N']}_T{param['T']}_t0{param['t0']}_alpha{param['alpha']}_tau{tau_str}.csv"

def load_and_merge_data(parameters):
    dfs = []
    for param in parameters:
        filename = generate_filename_from_params(param)
        df_temp = pd.read_csv(filename)

        # αとt0の値を使用して新しい列名を作成
        col_name = f"a{str(param['alpha']).replace('.','')}_t0{param['t0']}"
        df_temp = df_temp.rename(columns={"Z": col_name})
        dfs.append(df_temp)

    # すべてのデータフレームを1つにマージ
    df = pd.concat(dfs, axis=1)
    df = df.loc[:,~df.columns.duplicated()]
    return df
