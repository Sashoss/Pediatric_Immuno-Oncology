import re
import argparse
import logging
from pathlib import Path
from typing import Dict

import pandas as pd
import matplotlib.pyplot as plt

_FLOAT_RE = re.compile(r"-?\d+\.\d+")

CHAIN_TO_MOLECULE: Dict[str, str] = {
    "17": "SPP1-Heavy",
    "18": "SPP1-Light",
}


def parse_args() -> argparse.Namespace:
    """Parse command‐line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse gmx_MMPBSA ΔTOTAL energies and plot comparison boxplots."
    )
    parser.add_argument(
        "-i", "--input-dir",
        type=Path,
        default=Path("./out"),
        help="Root dir containing sample/Run_MD_chain*/FINAL_RESULTS_MMPBSA_*.dat"
    )
    parser.add_argument(
        "-o", "--output-dir",
        type=Path,
        default=Path("./out/plots"),
        help="Where to write summary plot(s)"
    )
    parser.add_argument(
        "-p", "--pattern",
        type=str,
        default="FINAL_RESULTS_MMPBSA_*.dat",
        help="Glob pattern for MMPBSA result files"
    )
    return parser.parse_args()


def read_delta_total(file_path: Path, label: str = "ΔTOTAL") -> float:
    """
    Read a single .dat file and extract the ΔTOTAL binding free energy.

    Raises:
        ValueError: if the label or a float cannot be found.
    """
    with file_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if label in line:
                m = _FLOAT_RE.search(line)
                if m:
                    return float(m.group())
                break

    raise ValueError(f"Could not parse '{label}' from {file_path}")


def collect_results(root_dir: Path, pattern: str) -> pd.DataFrame:
    """
    Walk through <root_dir>/<sample>/Run_MD_chain*/<pattern>, parse ΔTOTAL,
    and assemble a DataFrame with columns:
      ['sample', 'run', 'chain_id', 'molecule', 'energy']
    """
    records = []
    for sample_dir in sorted(root_dir.iterdir()):
        if not sample_dir.is_dir():
            continue

        sample_name = sample_dir.name
        for run_dir in sorted(sample_dir.glob("Run_MD_chain*")):
            run_name = run_dir.name
            for fn in sorted(run_dir.glob(pattern)):
                chain_id = fn.stem.split("_")[-1]
                molecule = CHAIN_TO_MOLECULE.get(chain_id, f"Chain-{chain_id}")

                try:
                    energy = read_delta_total(fn)
                except ValueError as e:
                    logging.warning(e)
                    continue

                records.append({
                    "sample": sample_name,
                    "run": run_name,
                    "chain_id": chain_id,
                    "molecule": molecule,
                    "energy": energy
                })

    if not records:
        logging.error("No valid ΔTOTAL entries found!")
    return pd.DataFrame.from_records(records)


def plot_boxplot(df: pd.DataFrame, out_path: Path) -> None:
    """
    Given a DataFrame with columns ['molecule','sample','energy'], draw
    a boxplot grouped first by molecule, then by sample.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(10, 6))
    ax = plt.gca()
    df.boxplot(
        column="energy",
        by=["molecule", "sample"],
        ax=ax,
        patch_artist=True,
        grid=False
    )
    plt.suptitle("DDG box-whiskers")  
    plt.title("Binding Free Energy (ΔG_bind) by Molecule and Sample", fontsize=14)
    plt.xlabel("Molecule / Sample", fontsize=12)
    plt.ylabel("ΔG_bind (kcal/mol)", fontsize=12)
    plt.xticks(rotation=45, ha="right")
    plt.grid(axis="y", linestyle="--", alpha=0.6)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    logging.info(f"Plot saved to {out_path}")


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    args = parse_args()
    logging.info(f"Searching for files in {args.input_dir} matching '{args.pattern}'")
    df = collect_results(args.input_dir, args.pattern)

    if df.empty:
        logging.error("No data to plot; exiting.")
        return

    summary = df.groupby(["molecule", "sample"])["energy"].describe()
    logging.info("\n" + summary.to_string())

    plot_fp = args.output_dir / "binding_energy_boxplot.png"
    plot_boxplot(df, plot_fp)


if __name__ == "__main__":
    main()
