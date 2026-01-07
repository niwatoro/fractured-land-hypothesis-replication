from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


@dataclass
class FanPlotConfig:
    quantiles: tuple[float, ...] = (0.1, 0.25, 0.5, 0.75, 0.9)


def plot_fan(
    data: np.ndarray,
    output_path: Path,
    *,
    title: str,
    xlabel: str = "Period",
    ylabel: str = "Herfindahl Index",
    config: FanPlotConfig | None = None,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if config is None:
        config = FanPlotConfig()
    qs = np.quantile(data, config.quantiles, axis=0)
    x = np.arange(1, data.shape[1] + 1)
    plt.figure(figsize=(8, 6))
    for low, high in zip(config.quantiles[: len(config.quantiles) // 2], reversed(config.quantiles)):
        if low >= high:
            continue
        low_idx = config.quantiles.index(low)
        high_idx = config.quantiles.index(high)
        plt.fill_between(x, qs[low_idx], qs[high_idx], alpha=0.2, color="steelblue")
    median_idx = config.quantiles.index(0.5)
    plt.plot(x, qs[median_idx], color="steelblue", linewidth=1.5)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_mean_lines(
    series: dict[str, np.ndarray],
    output_path: Path,
    *,
    xlabel: str = "Period",
    ylabel: str = "Herfindahl Index",
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    t_max = next(iter(series.values())).shape[0]
    x = np.arange(1, t_max + 1)
    plt.figure(figsize=(8, 4))
    for name, data in series.items():
        plt.plot(x, data, label=name, linewidth=1.5)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.ylim(0, 1)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def write_china_origin(cn_max_pixels: list[int], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame({"fid": cn_max_pixels})
    out = df.value_counts().reset_index()
    out.columns = ["fid", "Freq"]
    out.to_csv(output_path, index=False)


def plot_map_snapshot(
    map_gdf: pd.DataFrame,
    nation_ids: np.ndarray,
    output_path: Path,
    *,
    bounds: tuple[float, float, float, float] | None = None,
    rng: np.random.Generator | None = None,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if rng is None:
        rng = np.random.default_rng(0)
    sizes = np.bincount(nation_ids)
    jitter = rng.random(len(sizes)) * 1e-6
    order = np.argsort(sizes + jitter)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, len(sizes) + 1)
    color_id = (len(sizes) - ranks) % 9
    palette = [
        "#8dd3c7",
        "#ffffb3",
        "#bebada",
        "#fb8072",
        "#80b1d3",
        "#fdb462",
        "#b3de69",
        "#fccde5",
        "#d9d9d9",
    ]
    colors = [palette[color_id[nation_id]] for nation_id in nation_ids]

    fig, ax = plt.subplots(figsize=(8, 4))
    map_gdf.plot(ax=ax, color=colors, linewidth=0)
    if bounds is not None:
        xmin, ymin, xmax, ymax = bounds
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal")
    ax.set_axis_off()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
