from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

from fractured_land.config import build_config
from fractured_land.data import load_map_data
from fractured_land.plotting import plot_fan, plot_map_snapshot, plot_mean_lines, write_china_origin
from fractured_land.simulation import run_simulation


def _write_matrix(matrix: np.ndarray, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(matrix).to_csv(output_path, index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Fractured-Land Hypothesis replication in Python.")
    parser.add_argument("--scenario", default="Preferred", help="Scenario name from Start.R.")
    parser.add_argument("--runs", type=int, default=3, help="Number of simulation runs.")
    parser.add_argument("--seed", type=int, default=123, help="Random seed base.")
    parser.add_argument(
        "--map-dir",
        default="ReplicationData/R/map",
        help="Path to the map data directory.",
    )
    parser.add_argument(
        "--output-dir",
        default="Output",
        help="Output directory for CSVs and plots.",
    )
    parser.add_argument("--t-max", type=int, default=None, help="Override Tmax.")
    parser.add_argument("--plots", action="store_true", help="Generate plots.")
    parser.add_argument(
        "--map-steps",
        default="",
        help="Comma-separated time steps to render maps (e.g. 50,300,500).",
    )
    parser.add_argument(
        "--progress-interval",
        type=int,
        default=10,
        help="Print progress every N steps; use 0 to disable.",
    )
    args = parser.parse_args()

    config = build_config(args.scenario)
    if args.t_max is not None:
        config = replace(config, t_max=args.t_max)
    config = replace(config, try_no=args.runs)

    map_data = load_map_data(Path(args.map_dir), config.max_sea_dist)

    map_steps: set[int] = set()
    if args.map_steps:
        map_steps = {int(step) for step in args.map_steps.split(",") if step.strip()}

    results = []
    for run_idx in range(args.runs):
        rng = np.random.default_rng(args.seed + run_idx)
        progress_interval = args.progress_interval if args.progress_interval > 0 else None

        with tqdm(
            total=config.t_max,
            desc=f"run {run_idx + 1}/{args.runs}",
            unit="step",
        ) as bar:
            last_step = 0

            def _progress(run_id: int, step: int, t_max: int) -> None:
                nonlocal last_step
                delta = step - last_step
                if delta > 0:
                    bar.update(delta)
                    last_step = step

            result = run_simulation(
                map_data,
                config,
                rng=rng,
                run_id=run_idx + 1,
                progress_interval=progress_interval,
                progress_callback=_progress if progress_interval else None,
                map_steps=map_steps,
            )
        results.append(result)

    t_max = config.t_max
    herf_eu = np.vstack([res.metrics["EU"] for res in results])
    herf_cn = np.vstack([res.metrics["CN"] for res in results])
    herf_med = np.vstack([res.metrics["Med"] for res in results])
    herf_mideast = np.vstack([res.metrics["Mideast"] for res in results])
    herf_india = np.vstack([res.metrics["India"] for res in results])
    herf_seasia = np.vstack([res.metrics["SEAsia"] for res in results])
    herf_europax = np.vstack([res.metrics["EuropaX"] for res in results])

    cn_max = np.array([res.cn_max for res in results])
    eu_max = np.array([res.eu_max for res in results])
    cn_max_pixel = [res.cn_max_orgpixel for res in results]

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    herflast = pd.DataFrame(
        {
            "round": np.arange(1, args.runs + 1),
            "HerfCN": herf_cn[:, t_max - 1],
            "HerfEU": herf_eu[:, t_max - 1],
            "CNmax": cn_max,
            "EUmax": eu_max,
        }
    )
    herflast.to_csv(output_dir / "Herflast.csv", index=False)

    _write_matrix(herf_eu.T, output_dir / "herfEU.csv")
    _write_matrix(herf_cn.T, output_dir / "herfCN.csv")
    _write_matrix(herf_med.T, output_dir / "herfMed.csv")
    _write_matrix(herf_mideast.T, output_dir / "herfMideast.csv")
    _write_matrix(herf_india.T, output_dir / "herfIndia.csv")
    _write_matrix(herf_seasia.T, output_dir / "herfSEAsia.csv")
    _write_matrix(herf_europax.T, output_dir / "herfEuropaX.csv")

    intersect_cn = np.vstack([res.intersects["CN"] for res in results])
    intersect_eu = np.vstack([res.intersects["EU"] for res in results])
    union_cn = np.vstack([res.unions["CN"] for res in results])
    union_eu = np.vstack([res.unions["EU"] for res in results])

    _write_matrix(intersect_cn.T, output_dir / "IntstCN.csv")
    _write_matrix(intersect_eu.T, output_dir / "IntstEU.csv")
    _write_matrix(union_cn.T, output_dir / "UnionCN.csv")
    _write_matrix(union_eu.T, output_dir / "UnionEU.csv")
    _write_matrix((intersect_cn / union_cn).T, output_dir / "MatchIdxCN.csv")
    _write_matrix((intersect_eu / union_eu).T, output_dir / "MatchIdxEU.csv")

    for key in results[0].regimeno:
        series = np.vstack([res.regimeno[key] for res in results])
        _write_matrix(series.T, output_dir / f"regimeno.{key}.csv")
    for key in results[0].regimesizes:
        series = np.vstack([res.regimesizes[key].reshape(10, t_max) for res in results])
        _write_matrix(series.T, output_dir / f"regimesize.{key}.csv")

    write_china_origin(cn_max_pixel, output_dir / "ChinaOrigin.csv")

    if args.plots:
        plot_fan(herf_eu, output_dir / "EUplot.png", title="Europe")
        plot_fan(herf_cn, output_dir / "CNplot.png", title="China")
        plot_mean_lines(
            {
                "Europe": herf_eu.mean(axis=0),
                "China": herf_cn.mean(axis=0),
                "Mediterranean Region": herf_med.mean(axis=0),
            },
            output_dir / "meanEach_EUCN.png",
        )
        plot_mean_lines(
            {
                "Middle East": herf_mideast.mean(axis=0),
                "India": herf_india.mean(axis=0),
                "Southeast Asia": herf_seasia.mean(axis=0),
                "EuropaX": herf_europax.mean(axis=0),
            },
            output_dir / "meanEach_Others.png",
        )

    if map_steps:
        maps_dir = output_dir / "maps"
        bounds = map_data.eurasia_bounds if config.draw_eurasia == 1 else None
        for run_idx, result in enumerate(results, start=1):
            run_dir = maps_dir / f"run_{run_idx}"
            rng = np.random.default_rng(args.seed + run_idx)
            for step in sorted(map_steps):
                snapshot = result.map_snapshots.get(step)
                if snapshot is None:
                    continue
                plot_map_snapshot(
                    map_data.map_gdf,
                    snapshot,
                    run_dir / f"map_t{step}.png",
                    bounds=bounds,
                    rng=rng,
                )


if __name__ == "__main__":
    main()
