#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path
import re

import numpy as np
import matplotlib.pyplot as plt


def run_julia(main_jl: Path, nprocs: int, julia_bin: str, args_list: list[str]) -> None:
    cmd = [julia_bin, f"-p{nprocs}", str(main_jl)] + args_list
    print("\n[RUN]", " ".join(cmd))
    subprocess.run(cmd, check=True)


def newest_csv_under(dir_path: Path) -> Path:
    """
    main.jlが作ったディレクトリ配下で、最も新しいcsvを拾う。
    命名（SI prefixなど）に依存しないため安全。
    """
    csvs = list(dir_path.rglob("*.csv"))
    if not csvs:
        raise FileNotFoundError(f"No CSV found under: {dir_path}")
    csvs.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return csvs[0]


def load_M(csv_path: Path) -> np.ndarray:
    # 1行目がヘッダの場合があるのでまず skiprows=1 で試す
    try:
        data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
    except ValueError:
        # ヘッダが無い/別形式の可能性 → skipなしで再トライ
        data = np.loadtxt(csv_path, delimiter=",")

    return np.asarray(data).reshape(-1)



def plot_hist(M: np.ndarray, out_path: Path, bins: int, density: bool, title: str):
    plt.figure()
    plt.hist(M, bins=bins, density=density)
    plt.title(title)
    plt.xlabel("m")
    plt.ylabel("p(m)" if density else "count")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=200)
    plt.close()


def main():
    parser = argparse.ArgumentParser()

    # project paths
    parser.add_argument("--project_root", type=str,
                        default="~/shimizu/aco-on-spinglass/src/sampling/return_M",
                        help="return_M directory (contains main.jl).")
    parser.add_argument("--julia_bin", type=str, default="julia")
    parser.add_argument("--nprocs", type=int, default=8)

    # common sim params
    parser.add_argument("--seed", type=int, default=1000)
    parser.add_argument("--N", type=int, default=100)
    parser.add_argument("--tau", type=int, default=100)
    parser.add_argument("--sample", type=int, default=5000)
    parser.add_argument("--h", type=float, default=0.001)
    parser.add_argument("--J", type=float, default=0.1)
    parser.add_argument("--alpha_increment", type=float, default=1.0e-6)

    # target alphas (fig3: e.g. 0.8, alpha_c(h), 0.9)
    parser.add_argument("--alphas", type=str, default="0.8,0.9",
                        help="Comma-separated end_alpha values to run. Example: '0.8,0.9'")

    # histogram settings
    parser.add_argument("--bins", type=int, default=80)
    parser.add_argument("--density", action="store_true", help="Normalize histogram to density.")
    parser.add_argument("--fig_dir", type=str, default="./figs", help="Where to save figures (relative to project_root).")

    args = parser.parse_args()

    project_root = Path(args.project_root).expanduser().resolve()
    main_jl = project_root / "main.jl"
    if not main_jl.exists():
        raise FileNotFoundError(f"main.jl not found: {main_jl}")

    # main.jl の出力先は絶対パスで固定なので、ここも固定で見る
    data_root = Path("/home/mori-lab/shimizu/aco-on-spinglass/src/sampling/return_M/data")

    alphas = [float(x.strip()) for x in args.alphas.split(",") if x.strip()]

    # 実行のまとまり単位：seed ごとにディレクトリが切られるので、ここをベースに探索する
    seed_dir = data_root / f"seed{args.seed}"

    fig_dir = (project_root / args.fig_dir).resolve()
    fig_dir.mkdir(parents=True, exist_ok=True)

    for end_alpha in alphas:
        # ---- 1) Julia run (CSV生成) ----
        julia_args = [
            "--seed", str(args.seed),
            "--N", str(args.N),
            "--alpha", "0.0",  # alpha0（main.jlはこの値を使う）
            "--end_alpha", str(end_alpha),
            "--alpha_increment", str(args.alpha_increment),
            "--tau", str(args.tau),
            "--sample", str(args.sample),
            "--h", str(args.h),
            "--J", str(args.J),
        ]
        run_julia(main_jl=main_jl, nprocs=args.nprocs, julia_bin=args.julia_bin, args_list=julia_args)

        # ---- 2) CSV locate（SI prefix等に依存せず、最新のCSVを拾う）----
        # seed配下で最新のcsvを拾うが、alphaごとに確実に対応させるため、
        # "alpha{end_alpha}" を含むファイルを優先して拾う
        all_csvs = list(seed_dir.rglob("*.csv"))
        if not all_csvs:
            raise FileNotFoundError(f"No CSVs found under: {seed_dir}")

        # alphaフィルタ（ファイル名に alpha{end_alpha} が入る想定）
        alpha_key = f"alpha{end_alpha}"
        cand = [p for p in all_csvs if alpha_key in p.name]
        if cand:
            cand.sort(key=lambda p: p.stat().st_mtime, reverse=True)
            csv_path = cand[0]
        else:
            # fallback: newest overall
            csv_path = newest_csv_under(seed_dir)

        # ---- 3) Plot histogram ----
        M = load_M(csv_path)
        out_path = fig_dir / f"fig3_hist_alpha{end_alpha}.png"
        title = f"Histogram of m (alpha={end_alpha}, N={args.N}, sample={args.sample}, tau={args.tau}, h={args.h}, J={args.J})"
        plot_hist(M, out_path, bins=args.bins, density=args.density, title=title)

        print(f"[OK] alpha={end_alpha}")
        print(f"  CSV : {csv_path}")
        print(f"  FIG : {out_path}")
        print(f"  M size: {M.size}")

    print("\nDone.")


if __name__ == "__main__":
    main()
