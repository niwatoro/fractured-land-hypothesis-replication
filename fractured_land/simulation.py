from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable

import numpy as np
import pandas as pd

from fractured_land.config import ModelConfig
from fractured_land.data import MapData


@dataclass
class SimulationResult:
    metrics: dict[str, np.ndarray]
    regimesizes: dict[str, np.ndarray]
    regimeno: dict[str, np.ndarray]
    unions: dict[str, np.ndarray]
    intersects: dict[str, np.ndarray]
    map_snapshots: dict[int, np.ndarray]
    eu_max_orgpixel: int
    cn_max_orgpixel: int
    eu_max: float
    cn_max: float


def _which_equal(arr: np.ndarray, value: int) -> np.ndarray:
    return np.where(arr == value)[0]


def _quantile(arr: np.ndarray, q: float) -> float:
    return float(np.quantile(arr, q))


def _duplicated_random(values: Iterable[int], rng: np.random.Generator) -> np.ndarray:
    values = list(values)
    if len(values) <= 1:
        return np.zeros(len(values), dtype=bool)
    order = rng.permutation(len(values))
    seen = set()
    dup = [False] * len(values)
    for idx in order:
        val = values[idx]
        if val in seen:
            dup[idx] = True
        else:
            seen.add(val)
    return np.array(dup, dtype=bool)


def _length_na(seq: Iterable[int] | None) -> int:
    if seq is None:
        return 0
    if isinstance(seq, float) and np.isnan(seq):
        return 0
    return len(seq)


def _theta(
    pixel_idx: np.ndarray,
    *,
    theta_rugged: float,
    theta_cold: float,
    theta_hot: float,
    theta_steppe: float,
    pixel_ruggedness: np.ndarray,
    pixel_cold: np.ndarray,
    pixel_hot: np.ndarray,
    pixel_steppe: np.ndarray,
) -> np.ndarray:
    val = (
        theta_rugged * pixel_ruggedness[pixel_idx]
        + theta_cold * pixel_cold[pixel_idx]
        + theta_hot * pixel_hot[pixel_idx]
        + theta_steppe * pixel_steppe[pixel_idx]
    )
    return np.maximum(val, 0)


def _pixel_land_connected(
    start: int,
    territory: Iterable[int],
    neighbors_land: list[list[int]],
) -> list[int]:
    territory_set = set(territory)
    if start not in territory_set:
        return []
    connected = {start}
    frontier = [start]
    while frontier:
        new_frontier: list[int] = []
        for node in frontier:
            for neighbor in neighbors_land[node]:
                if neighbor in territory_set and neighbor not in connected:
                    connected.add(neighbor)
                    new_frontier.append(neighbor)
        frontier = new_frontier
    return list(connected)


def _split_disconnect(
    pixels: Iterable[int],
    neighbors: list[list[int]],
) -> list[list[int]]:
    pixels_list = list(pixels)
    if len(pixels_list) <= 1:
        return [pixels_list]
    remaining = set(pixels_list)
    components: list[list[int]] = []
    while remaining:
        start = next(iter(remaining))
        component = {start}
        frontier = [start]
        remaining.remove(start)
        while frontier:
            new_frontier = []
            for node in frontier:
                for neighbor in neighbors[node]:
                    if neighbor in remaining:
                        remaining.remove(neighbor)
                        component.add(neighbor)
                        new_frontier.append(neighbor)
            frontier = new_frontier
        components.append(list(component))
    return components


def _build_pixel_belong(nation_territory: list[list[int]], pixelno: int) -> np.ndarray:
    pixel_belong = np.full(pixelno, -1, dtype=int)
    for nation_id, pixels in enumerate(nation_territory):
        for pixel in pixels:
            pixel_belong[pixel] = nation_id
    return pixel_belong


def _region_metrics(
    nation_territory: list[list[int]],
    region_pixels: np.ndarray,
    *,
    t_idx: int,
    hhi: np.ndarray,
    regimeno: np.ndarray,
    regimesize: np.ndarray,
) -> list[int]:
    region_set = set(region_pixels)
    territory_sizes = [len([p for p in territory if p in region_set]) for territory in nation_territory]
    pixelno = len(region_pixels)
    if pixelno == 0:
        hhi[t_idx] = 0
        regimeno[t_idx] = 0
        return territory_sizes
    hhi[t_idx] = float(sum((size / pixelno) ** 2 for size in territory_sizes))
    regimeno[t_idx] = int(sum(size > 0 for size in territory_sizes))
    top = sorted(territory_sizes, reverse=True)[: min(10, regimeno[t_idx])]
    if top:
        regimesize[: len(top), t_idx] = np.array(top, dtype=float)
    return territory_sizes


def _apply_transfer(
    nation_territory: list[list[int]],
    nation_orgpixel: list[int],
    lost_by_nation: dict[int, list[int]],
    won_by_nation: dict[int, list[int]],
) -> tuple[list[list[int]], list[int]]:
    for nation, pixels in lost_by_nation.items():
        if not pixels:
            continue
        lost_set = set(pixels)
        nation_territory[nation] = [p for p in nation_territory[nation] if p not in lost_set]

    for nation, pixels in won_by_nation.items():
        if not pixels:
            continue
        nation_territory[nation].extend(pixels)

    surviving = [idx for idx, terr in enumerate(nation_territory) if len(terr) > 0]
    nation_territory = [nation_territory[idx] for idx in surviving]
    nation_orgpixel = [nation_orgpixel[idx] for idx in surviving]
    return nation_territory, nation_orgpixel


def run_simulation(
    map_data: MapData,
    config: ModelConfig,
    *,
    rng: np.random.Generator,
    run_id: int,
    progress_interval: int | None = None,
    progress_callback: Callable[[int, int, int], None] | None = None,
    map_steps: set[int] | None = None,
) -> SimulationResult:
    eadta = map_data.eadta.copy()
    pixelno = len(eadta)
    scenario = config.scenario

    for col in [
        "coastal",
        "medcoast",
        "tmin",
        "tmax",
        "startKK10",
        "startDSMW",
        "startloess",
        "startHSWD",
        "river",
        "medsea",
    ]:
        if col in eadta.columns:
            eadta[col] = eadta[col].fillna(0)

    pixel_coast = eadta["coastal"].to_numpy(dtype=float)
    pixel_medcoast = eadta["medcoast"].to_numpy(dtype=float)

    pixel_eu = _which_equal(eadta["Europe"].to_numpy(dtype=float), 1)
    pixel_cn = _which_equal(eadta["China"].to_numpy(dtype=float), 1)
    pixel_india = _which_equal(eadta["indiancont"].to_numpy(dtype=float), 1)
    pixel_mideast = _which_equal(eadta["mideast"].to_numpy(dtype=float), 1)
    pixel_seasia = _which_equal(eadta["seasia"].to_numpy(dtype=float), 1)
    pixel_europax = _which_equal(eadta["Europe"].to_numpy(dtype=float), 1)

    pixel_africa_et = _which_equal(eadta["AfricaET"].to_numpy(dtype=float), 1)
    pixel_africa_wt = _which_equal(eadta["AfricaWT"].to_numpy(dtype=float), 1)
    pixel_africa_me = _which_equal(eadta["AfricaME"].to_numpy(dtype=float), 1)
    pixel_africa_nh = _which_equal(eadta["AfricaNH"].to_numpy(dtype=float), 1)
    pixel_africa_sh = _which_equal(eadta["AfricaSH"].to_numpy(dtype=float), 1)

    pixel_america_nh = _which_equal(eadta["AmericaNH"].to_numpy(dtype=float), 1)
    pixel_america_sh = _which_equal(eadta["AmericaSH"].to_numpy(dtype=float), 1)
    pixel_america_cl = _which_equal(eadta["AmericaCL"].to_numpy(dtype=float), 1)

    if "medsea" in eadta.columns:
        medsea_region = _which_equal(eadta["medsea"].to_numpy(dtype=float), 1)
    else:
        medsea_region = np.array([], dtype=int)

    pixel_ruggedness = eadta["ElevationS"].to_numpy(dtype=float)

    pixel_temp = eadta["tmin"].to_numpy(dtype=float)
    pixel_cold = np.log(np.maximum(9 - pixel_temp, 1))
    theta_cold = scenario.theta_cold / np.log(9 - _quantile(pixel_temp, 0.1))

    pixel_temp2 = eadta["tmax"].to_numpy(dtype=float)
    pixel_hot = np.log(np.maximum(pixel_temp2 - 21, 1))
    theta_hot = scenario.theta_hot / np.log(_quantile(pixel_temp2, 0.9) - 21)

    theta_rugged = scenario.theta_rugged_90th / _quantile(pixel_ruggedness, 0.9)
    theta_sea = 0.0
    theta0 = 1
    if scenario.theta_rugged_90th == 0 and scenario.theta_cold == 0 and scenario.theta_hot == 0:
        theta0 = 0

    if scenario.yield_data == "YKK10":
        pixel_y = eadta["YKK10"].to_numpy(dtype=float)
    elif scenario.yield_data == "YCSI":
        pixel_y = eadta["YCSI"].to_numpy(dtype=float)
    elif scenario.yield_data == "YGAEZ":
        pixel_y = eadta["YGAEZ"].to_numpy(dtype=float)
    elif scenario.yield_data == "YGAEZ4":
        pixel_y = eadta["YGAEZ4"].to_numpy(dtype=float)
    elif scenario.uniform_y == 1:
        pixel_y = np.full(pixelno, 0.5, dtype=float)
    else:
        raise ValueError("Yield data is not specified and uniform_y is disabled.")

    pixel_y = np.nan_to_num(pixel_y, nan=0.0)

    pixel_steppeeast = eadta["steppeeast"].to_numpy(dtype=float)
    id_steppeeast = _which_equal(pixel_steppeeast, 1)
    pixel_steppe = eadta["steppe_all"].to_numpy(dtype=float)
    id_steppe = _which_equal(pixel_steppe, 1)

    pixel_river = eadta["river"].to_numpy(dtype=float) if scenario.river == 1 else None
    pixel_riverdummy = None
    if scenario.river == 1:
        pixel_riverdummy = (pixel_river > 0).astype(int)

    nation_territory = [list([idx]) for idx in range(pixelno)]
    if config.start_with_regime == 1 and "regimeid" in eadta.columns:
        regime_ids = eadta["regimeid"].to_numpy(dtype=float)
        regime_series = pd.Series(np.arange(pixelno), index=regime_ids)
        nation_territory = [list(vals) for _, vals in regime_series.groupby(level=0)]

    nation_orgpixel = list(range(pixelno))

    chi_cold = (1 / (1 - scenario.ydiscount_cold) - 1) / np.log(max(9 + 4, 1))
    chi_hot = (1 / (1 - scenario.ydiscount_hot) - 1) / np.log(max(29.5 - 21, 1))
    pixel_y = pixel_y / (1 + chi_cold * pixel_cold + chi_hot * pixel_hot)
    pixel_y_unrest = pixel_y.copy()

    if scenario.uniform_y == 1:
        alpha = 1
    else:
        anchor_pixels = np.union1d(pixel_eu, pixel_cn)
        alpha = 1 / _quantile(pixel_y_unrest[anchor_pixels], 0.95)

    pixel_y = pixel_y + config.ymin

    t_max = config.t_max
    hhi = {
        "EU": np.zeros(t_max),
        "CN": np.zeros(t_max),
        "Med": np.zeros(t_max),
        "Mideast": np.zeros(t_max),
        "India": np.zeros(t_max),
        "SEAsia": np.zeros(t_max),
        "EuropaX": np.zeros(t_max),
    }
    regimeno = {
        "EU": np.zeros(t_max),
        "CN": np.zeros(t_max),
        "Med": np.zeros(t_max),
        "Mideast": np.zeros(t_max),
        "India": np.zeros(t_max),
        "SEAsia": np.zeros(t_max),
        "EuropaX": np.zeros(t_max),
        "AfricaET": np.zeros(t_max),
        "AfricaWT": np.zeros(t_max),
        "AfricaME": np.zeros(t_max),
        "AfricaNH": np.zeros(t_max),
        "AfricaSH": np.zeros(t_max),
        "AmericaNH": np.zeros(t_max),
        "AmericaSH": np.zeros(t_max),
        "AmericaCL": np.zeros(t_max),
    }
    regimesize = {
        key: np.zeros((10, t_max)) for key in regimeno.keys()
    }
    unions = {"EU": np.zeros(t_max), "CN": np.zeros(t_max)}
    intersects = {"EU": np.zeros(t_max), "CN": np.zeros(t_max)}

    maize1 = _which_equal(eadta["Maize1"].to_numpy(dtype=float), 1)
    maize2 = _which_equal(
        (eadta["Maize2"].to_numpy(dtype=float) == 1)
        & (eadta["Maize1"].to_numpy(dtype=float) != 1),
        True,
    )
    maize3 = _which_equal(
        (eadta["Maize3"].to_numpy(dtype=float) == 1)
        & (eadta["Maize1"].to_numpy(dtype=float) != 1)
        & (eadta["Maize2"].to_numpy(dtype=float) != 1),
        True,
    )

    america = _which_equal(
        (eadta["Americas"].to_numpy(dtype=float) == 1)
        & (eadta["Maize1"].to_numpy(dtype=float) != 1)
        & (eadta["Maize2"].to_numpy(dtype=float) != 1)
        & (eadta["Maize3"].to_numpy(dtype=float) != 1),
        True,
    )
    africa = _which_equal(eadta["Africa"].to_numpy(dtype=float), 1)
    australia = _which_equal(eadta["Australia"].to_numpy(dtype=float), 1)
    japan = _which_equal(eadta["Japan"].to_numpy(dtype=float), 1)

    startyr = np.zeros(pixelno, dtype=int)
    startyr[maize1] = config.maize1_start
    startyr[maize2] = config.maize2_start
    startyr[maize3] = config.maize3_start
    startyr[america] = config.america_start
    startyr[africa] = config.africa_start
    startyr[australia] = config.australia_start
    startyr[japan] = config.japan_start

    timer = 0
    nationno = len(nation_territory)
    eu_max = 0.0
    cn_max = 0.0
    eu_max_orgpixel = 0
    cn_max_orgpixel = 0
    map_snapshots: dict[int, np.ndarray] = {}
    map_steps_set = set(map_steps) if map_steps else set()

    increment = None
    increment2 = None
    pixel_y_500 = None
    pixel_y_unrest_500 = None
    pixel_single_set = set(map_data.pixel_single)

    while True:
        if timer >= t_max:
            break
        if nationno <= 1:
            for key in hhi:
                hhi[key][timer:] = 1
            break

        timer += 1
        t_idx = timer - 1
        if progress_interval and progress_callback:
            if timer == 1 or timer % progress_interval == 0 or timer == t_max:
                progress_callback(run_id, timer, t_max)

        ratio_y0 = np.ones(pixelno)
        if scenario.bc1000_ratio == "startKK10":
            ratio_y0 = eadta["startKK10"].to_numpy(dtype=float)
        elif scenario.bc1000_ratio == "startDSMW":
            ratio_y0 = eadta["startDSMW"].to_numpy(dtype=float)
        elif scenario.bc1000_ratio == "startloess":
            ratio_y0 = eadta["startloess"].to_numpy(dtype=float)
        elif scenario.bc1000_ratio == "startHSWD":
            ratio_y0 = eadta["startHSWD"].to_numpy(dtype=float)

        if timer == 1:
            pixel_y_500 = pixel_y.copy()
            pixel_y_unrest_500 = pixel_y_unrest.copy()
            increment = pixel_y * (1 - ratio_y0) / 500
            pixel_y = pixel_y * ratio_y0
            pixel_y = np.maximum(pixel_y, 0.000001)
            increment2 = pixel_y_unrest * (1 - ratio_y0) / 500
            pixel_y_unrest = pixel_y_unrest * ratio_y0
            pixel_y_unrest = np.maximum(pixel_y_unrest, 0.000001)
        else:
            pixel_y = np.minimum(pixel_y + increment * (timer >= startyr), pixel_y_500)
            pixel_y_unrest = np.minimum(
                pixel_y_unrest + increment2 * (timer >= startyr), pixel_y_unrest_500
            )

        if scenario.y_steppe is not None:
            pixel_y[id_steppe] = scenario.y_steppe
            pixel_y_unrest[id_steppe] = scenario.y_steppe
        if config.y_steppeeast is not None:
            pixel_y[id_steppeeast] = config.y_steppeeast
            pixel_y_unrest[id_steppeeast] = config.y_steppeeast

        nation_resource = np.array([sum(pixel_y[pixels]) for pixels in nation_territory])
        pixel_belong = _build_pixel_belong(nation_territory, pixelno)

        pixel_unrestprob = np.clip(alpha * pixel_y_unrest, config.unrest_min, 1)
        pixel_unrestprob[timer < startyr] = config.unrest_min
        unrest_flags = rng.random(pixelno) < pixel_unrestprob
        pixel_unrest = np.where(unrest_flags)[0]
        pixel_unrest = np.array(
            [p for p in pixel_unrest if p not in pixel_single_set], dtype=int
        )

        if len(pixel_unrest) == 0:
            nationno = len(nation_territory)
            continue

        pixel_unrest_neighbors = [map_data.neighbors[p] for p in pixel_unrest]

        noseawar_id: np.ndarray = np.array([], dtype=int)
        if scenario.conflict_mech in {"random", "neighborY"}:
            neighbors_land = [map_data.neighbors_land[p] for p in pixel_unrest]
            neighbors_sea = [
                [n for n in nb if n not in set(land)]
                for nb, land in zip(pixel_unrest_neighbors, neighbors_land)
            ]
            landborder = np.array([len(land) for land in neighbors_land], dtype=float)

            if scenario.conflict_mech == "random":
                land_probs = [
                    np.full(len(land), 1 / len(land)) if land else np.array([])
                    for land in neighbors_land
                ]
            else:
                land_probs = [
                    pixel_y[land] / np.sum(pixel_y[land]) if len(land) else np.array([])
                    for land in neighbors_land
                ]
            land_probs = [
                probs * border / 6 if len(probs) else np.array([])
                for probs, border in zip(land_probs, landborder)
            ]

            if scenario.conflict_mech == "random":
                sea_probs = [
                    np.full(len(sea), 1 / len(sea)) if sea else np.array([])
                    for sea in neighbors_sea
                ]
            else:
                sea_probs = [
                    pixel_y[sea] / np.sum(pixel_y[sea]) if len(sea) else np.array([])
                    for sea in neighbors_sea
                ]
            sea_probs = [
                scenario.alpha_sea * probs * (6 - border) / 6 if len(probs) else np.array([])
                for probs, border in zip(sea_probs, landborder)
            ]

            combined_neighbors = [
                list(land) + list(sea)
                for land, sea in zip(neighbors_land, neighbors_sea)
            ]
            combined_probs = [
                np.concatenate([land_p, sea_p]) if len(land_p) + len(sea_p) > 0 else np.array([])
                for land_p, sea_p in zip(land_probs, sea_probs)
            ]
            combined_probs = [
                probs * min(len(probs), 6) / 6 if len(probs) else np.array([])
                for probs in combined_probs
            ]

            pixel_unrest_target = []
            for nbs, probs in zip(combined_neighbors, combined_probs):
                if len(nbs) == 0:
                    pixel_unrest_target.append(-1)
                    continue
                if len(nbs) == 1:
                    pixel_unrest_target.append(nbs[0])
                else:
                    pixel_unrest_target.append(
                        int(rng.choice(nbs, size=1, p=probs / probs.sum())[0])
                    )
            pixel_unrest_target = np.array(pixel_unrest_target, dtype=int)

            unrest_belongnation = pixel_belong[pixel_unrest]
            valid_target = pixel_unrest_target >= 0
            if not np.all(valid_target):
                pixel_unrest = pixel_unrest[valid_target]
                pixel_unrest_target = pixel_unrest_target[valid_target]
                unrest_belongnation = unrest_belongnation[valid_target]
                landborder = landborder[valid_target]
            noseawar_prob = (6 - landborder) / 6 * (1 - scenario.alpha_sea)
            noseawar_id = np.where(rng.random(len(noseawar_prob)) < noseawar_prob)[0]
        elif scenario.conflict_mech == "expgain":
            unrest_belongnation = pixel_belong[pixel_unrest]
            unrest_neighbors_belong = [
                pixel_belong[nb] if len(nb) else np.array([], dtype=int)
                for nb in pixel_unrest_neighbors
            ]
            enemy_neighbors = [
                [n for n, nation in zip(nbs, belongs) if nation != self_nation]
                for nbs, belongs, self_nation in zip(
                    pixel_unrest_neighbors, unrest_neighbors_belong, unrest_belongnation
                )
            ]
            valid_mask = [len(enemies) > 0 for enemies in enemy_neighbors]
            pixel_unrest = pixel_unrest[valid_mask]
            unrest_belongnation = unrest_belongnation[valid_mask]
            enemy_neighbors = [enemies for enemies, keep in zip(enemy_neighbors, valid_mask) if keep]

            enemy_belong = [pixel_belong[enemies] for enemies in enemy_neighbors]
            enemy_cost = [
                np.array(
                    [
                        _sum_cost_y(
                            nation_idx,
                            target,
                            nation_territory,
                            map_data.neighbors_land,
                            pixel_y,
                            scenario.sigma,
                        )
                        for nation_idx, target in zip(enemy_nation, enemies)
                    ]
                )
                for enemy_nation, enemies in zip(enemy_belong, enemy_neighbors)
            ]
            self_cost = [
                np.array(
                    [
                        _sum_cost_y(
                            self_nation,
                            target,
                            nation_territory,
                            map_data.neighbors_land,
                            pixel_y,
                            scenario.sigma,
                        )
                        for target in enemies
                    ]
                )
                for self_nation, enemies in zip(unrest_belongnation, enemy_neighbors)
            ]

            enemy_land_neighbors = [
                [n for n in enemies if n in set(map_data.neighbors_land[pixel])]
                for pixel, enemies in zip(pixel_unrest, enemy_neighbors)
            ]
            enemy_sea_vector = [
                np.array([1 if n not in set(land) else 0 for n in enemies])
                for enemies, land in zip(enemy_neighbors, enemy_land_neighbors)
            ]

            pixel_unrest_expgain = [
                (1 - sea_vec * (1 - scenario.alpha_sea))
                * (attack_y / (attack_y + defense_y))
                * pixel_y[np.array(enemies)]
                / (1 + _theta(np.array(enemies), theta_rugged=theta_rugged, theta_cold=theta_cold,
                              theta_hot=theta_hot, theta_steppe=scenario.theta_steppe,
                              pixel_ruggedness=pixel_ruggedness, pixel_cold=pixel_cold,
                              pixel_hot=pixel_hot, pixel_steppe=pixel_steppe))
                for enemies, sea_vec, attack_y, defense_y in zip(
                    enemy_neighbors, enemy_sea_vector, self_cost, enemy_cost
                )
            ]

            pixel_unrest_target = np.array(
                [enemies[int(np.argmax(gain))] for enemies, gain in zip(enemy_neighbors, pixel_unrest_expgain)],
                dtype=int,
            )

            seawar_id = [
                idx
                for idx, (target, land) in enumerate(zip(pixel_unrest_target, enemy_land_neighbors))
                if target not in land
            ]
            noseawar_id = np.array(
                [idx for idx in seawar_id if rng.random() < (1 - scenario.alpha_sea)],
                dtype=int,
            )
        else:
            raise ValueError(f"Unknown conflict mechanism: {scenario.conflict_mech}")

        if len(pixel_unrest) == 0:
            nationno = len(nation_territory)
            continue

        unrest_target_belong = pixel_belong[pixel_unrest_target]
        warid = np.array(
            [
                idx
                for idx, (a, d) in enumerate(zip(unrest_belongnation, unrest_target_belong))
                if a != d
            ],
            dtype=int,
        )
        if len(noseawar_id) > 0:
            warid = np.array([idx for idx in warid if idx not in set(noseawar_id)], dtype=int)

        if len(warid) > 0:
            pixel_attack = pixel_unrest[warid]
            pixel_attack_target = pixel_unrest_target[warid]
            nation_attack = unrest_belongnation[warid]
            nation_defense = unrest_target_belong[warid]

            defense_connect = [
                _pixel_land_connected(attacker, nation_territory[defender], map_data.neighbors_land)
                for attacker, defender in zip(pixel_attack, nation_defense)
            ]
            attack_connect = [
                _pixel_land_connected(target, nation_territory[attacker], map_data.neighbors_land)
                for target, attacker in zip(pixel_attack_target, nation_attack)
            ]
            seawar = np.array(
                [len(x) == 0 and len(y) == 0 for x, y in zip(defense_connect, attack_connect)],
                dtype=bool,
            )

            defense_connect_y = np.array([sum(pixel_y[x]) for x in defense_connect])
            defense_unconnect_y = nation_resource[nation_defense] - defense_connect_y
            defense_y_est = defense_connect_y + scenario.sigma * defense_unconnect_y

            attack_connect_y = np.array([sum(pixel_y[x]) for x in attack_connect])
            attack_unconnect_y = nation_resource[nation_attack] - attack_connect_y
            attack_y_est = attack_connect_y + scenario.sigma * attack_unconnect_y

            oppos_y_est_sum = np.zeros(nationno)
            for nation_id, val in zip(nation_defense, attack_y_est):
                oppos_y_est_sum[nation_id] += val
            for nation_id, val in zip(nation_attack, defense_y_est):
                oppos_y_est_sum[nation_id] += val

            attack_y_ratio = defense_y_est / np.maximum(oppos_y_est_sum[nation_attack], 1e-9)
            attack_y_cost = attack_y_est * attack_y_ratio
            defense_y_ratio = attack_y_est / np.maximum(oppos_y_est_sum[nation_defense], 1e-9)
            defense_y_cost = defense_y_est * defense_y_ratio

            theta_attack = _theta(
                pixel_attack,
                theta_rugged=theta_rugged,
                theta_cold=theta_cold,
                theta_hot=theta_hot,
                theta_steppe=scenario.theta_steppe,
                pixel_ruggedness=pixel_ruggedness,
                pixel_cold=pixel_cold,
                pixel_hot=pixel_hot,
                pixel_steppe=pixel_steppe,
            )
            theta_defense = _theta(
                pixel_attack_target,
                theta_rugged=theta_rugged,
                theta_cold=theta_cold,
                theta_hot=theta_hot,
                theta_steppe=scenario.theta_steppe,
                pixel_ruggedness=pixel_ruggedness,
                pixel_cold=pixel_cold,
                pixel_hot=pixel_hot,
                pixel_steppe=pixel_steppe,
            )
            theta_max = np.maximum(theta_attack, theta_defense)

            if scenario.river == 1 and pixel_riverdummy is not None:
                anyriver = (pixel_riverdummy[pixel_attack] + pixel_riverdummy[pixel_attack_target]) > 0
                riverconnected = (
                    (pixel_riverdummy[pixel_attack] > 0) & (pixel_riverdummy[pixel_attack_target] > 0)
                )
                theta_max = theta_max - riverconnected * theta_max + (~riverconnected) * anyriver * 2

            if theta0 == 0:
                theta_max = theta_max * 0
            theta_max = theta_max + seawar * theta_sea

            attack_stepborder = np.array(
                [
                    scenario.psi_steppe
                    if (
                        nation_orgpixel[attacker] in set(id_steppeeast)
                        and len(set(nation_territory[attacker]).intersection(id_steppeeast)) > 0
                    )
                    else 1
                    for attacker in nation_attack
                ],
                dtype=float,
            )
            defense_stepborder = np.array(
                [
                    scenario.psi_steppe
                    if (
                        nation_orgpixel[defender] in set(id_steppeeast)
                        and len(set(nation_territory[defender]).intersection(id_steppeeast)) > 0
                    )
                    else 1
                    for defender in nation_defense
                ],
                dtype=float,
            )

            denom = attack_stepborder * attack_y_cost + defense_stepborder * defense_y_cost
            winprob = attack_stepborder * attack_y_cost / denom / (1 + theta_max)
            loseprob = defense_stepborder * defense_y_cost / denom / (1 + theta_max)
            tieprob = np.maximum(1 - winprob - loseprob, 0)

            draw = rng.random(len(winprob))
            outcome_win = draw < winprob
            outcome_lose = (draw >= winprob) & (draw < winprob + loseprob)
            outcome_tie = ~(outcome_win | outcome_lose)

            win_fid = np.where(outcome_win)[0]
            lose_fid = np.where(outcome_lose)[0]

            if len(win_fid) > 0:
                dup_mask = _duplicated_random(pixel_attack_target[win_fid], rng)
                win_ufid = win_fid[~dup_mask]
                if len(win_ufid) != len(win_fid):
                    dup_ids = np.setdiff1d(win_fid, win_ufid)
                    outcome_win[dup_ids] = False
                    outcome_tie[dup_ids] = True
                    win_fid = np.where(outcome_win)[0]
                    lose_fid = np.where(outcome_lose)[0]

            if len(lose_fid) > 0 and len(win_fid) > 0:
                combined = np.concatenate([pixel_attack[lose_fid], pixel_attack_target[win_fid]])
                dup_mask = _duplicated_random(combined, rng)
                lose_keep = dup_mask[: len(lose_fid)] == 0
                win_keep = dup_mask[len(lose_fid) :] == 0

                lose_ufid = lose_fid[lose_keep]
                win_ufid = win_fid[win_keep]
                if len(lose_ufid) != len(lose_fid) or len(win_ufid) != len(win_fid):
                    drop_lose = np.setdiff1d(lose_fid, lose_ufid)
                    drop_win = np.setdiff1d(win_fid, win_ufid)
                    outcome_lose[drop_lose] = False
                    outcome_win[drop_win] = False
                    outcome_tie[drop_lose] = True
                    outcome_tie[drop_win] = True
                    win_fid = np.where(outcome_win)[0]
                    lose_fid = np.where(outcome_lose)[0]

            if len(win_fid) > 0 or len(lose_fid) > 0:
                lost_by_nation: dict[int, list[int]] = {}
                won_by_nation: dict[int, list[int]] = {}
                for idx in win_fid:
                    defender = nation_defense[idx]
                    attacker = nation_attack[idx]
                    target = pixel_attack_target[idx]
                    lost_by_nation.setdefault(defender, []).append(int(target))
                    won_by_nation.setdefault(attacker, []).append(int(target))
                for idx in lose_fid:
                    attacker = nation_attack[idx]
                    defender = nation_defense[idx]
                    lost_pixel = pixel_attack[idx]
                    lost_by_nation.setdefault(attacker, []).append(int(lost_pixel))
                    won_by_nation.setdefault(defender, []).append(int(lost_pixel))

                nation_territory, nation_orgpixel = _apply_transfer(
                    nation_territory, nation_orgpixel, lost_by_nation, won_by_nation
                )

        if theta0 != 0:
            pixel_belong = _build_pixel_belong(nation_territory, pixelno)
            nation_pixelborder = []
            nation_pixelborder_in = []
            nation_dellandborder = []
            nation_coastL = []
            nation_medcoastL = []
            for territory in nation_territory:
                terr_set = set(territory)
                border_neighbors = set()
                for pix in territory:
                    border_neighbors.update(map_data.neighbors[pix])
                pixelborder = list(border_neighbors.difference(terr_set))
                nation_pixelborder.append(pixelborder)

                border_in_neighbors = set()
                for pix in pixelborder:
                    border_in_neighbors.update(map_data.neighbors[pix])
                pixelborder_in = list(border_in_neighbors.intersection(terr_set))
                nation_pixelborder_in.append(pixelborder_in)

                del_landborder = list(terr_set.difference(pixelborder_in))
                nation_dellandborder.append(del_landborder)
                nation_coastL.append(float(np.sum(pixel_coast[del_landborder])) if del_landborder else 0.0)
                nation_medcoastL.append(
                    float(np.sum(pixel_medcoast[del_landborder])) if del_landborder else 0.0
                )

            nation_borderlength = np.array([len(x) for x in nation_pixelborder_in], dtype=float)
            nation_allborderL = nation_borderlength + np.array(nation_coastL) - np.array(nation_medcoastL)
            nation_size2L = [idx for idx, terr in enumerate(nation_territory) if len(terr) > 1]

            pixel_seperate = {}
            for idx in nation_size2L:
                pixelvec = nation_pixelborder_in[idx]
                if not pixelvec:
                    continue
                probs = config.beta * _theta(
                    np.array(pixelvec),
                    theta_rugged=theta_rugged,
                    theta_cold=theta_cold,
                    theta_hot=theta_hot,
                    theta_steppe=scenario.theta_steppe,
                    pixel_ruggedness=pixel_ruggedness,
                    pixel_cold=pixel_cold,
                    pixel_hot=pixel_hot,
                    pixel_steppe=pixel_steppe,
                ) * nation_allborderL[pixel_belong[pixelvec]]
                probs = np.minimum(probs, 1)
                separate = [p for p, keep in zip(pixelvec, rng.random(len(pixelvec)) < probs) if keep]
                if separate:
                    pixel_seperate[idx] = separate

            for idx, separated in pixel_seperate.items():
                nation_territory[idx] = [p for p in nation_territory[idx] if p not in set(separated)]
                for pixel in separated:
                    nation_territory.append([pixel])
                    nation_orgpixel.append(nation_orgpixel[idx])

            surviving = [idx for idx, terr in enumerate(nation_territory) if len(terr) > 0]
            nation_territory = [nation_territory[idx] for idx in surviving]
            nation_orgpixel = [nation_orgpixel[idx] for idx in surviving]

        if config.shockprob_regime > 0:
            breaks = rng.random(len(nation_territory)) < config.shockprob_regime
            broken = [idx for idx, val in enumerate(breaks) if val]
            if broken:
                new_territory = []
                new_orgpixel = []
                for idx, terr in enumerate(nation_territory):
                    if idx in broken:
                        for pixel in terr:
                            new_territory.append([pixel])
                            new_orgpixel.append(nation_orgpixel[idx])
                    else:
                        new_territory.append(terr)
                        new_orgpixel.append(nation_orgpixel[idx])
                nation_territory = new_territory
                nation_orgpixel = new_orgpixel

        if config.shockprob_general > 0 and rng.random() < config.shockprob_general:
            nation_territory = [[idx] for idx in range(pixelno)]
            nation_orgpixel = list(range(pixelno))

        split_territory = []
        split_orgpixel = []
        for idx, terr in enumerate(nation_territory):
            components = _split_disconnect(terr, map_data.neighbors)
            if len(components) > 1:
                for comp in components:
                    split_territory.append(comp)
                    split_orgpixel.append(nation_orgpixel[idx])
            else:
                split_territory.append(terr)
                split_orgpixel.append(nation_orgpixel[idx])
        nation_territory = split_territory
        nation_orgpixel = split_orgpixel

        nationno = len(nation_territory)

        territory_eu = _region_metrics(
            nation_territory, pixel_eu, t_idx=t_idx, hhi=hhi["EU"],
            regimeno=regimeno["EU"], regimesize=regimesize["EU"]
        )
        territory_cn = _region_metrics(
            nation_territory, pixel_cn, t_idx=t_idx, hhi=hhi["CN"],
            regimeno=regimeno["CN"], regimesize=regimesize["CN"]
        )
        _region_metrics(
            nation_territory, medsea_region, t_idx=t_idx, hhi=hhi["Med"],
            regimeno=regimeno["Med"], regimesize=regimesize["Med"]
        )
        _region_metrics(
            nation_territory, pixel_mideast, t_idx=t_idx, hhi=hhi["Mideast"],
            regimeno=regimeno["Mideast"], regimesize=regimesize["Mideast"]
        )
        _region_metrics(
            nation_territory, pixel_india, t_idx=t_idx, hhi=hhi["India"],
            regimeno=regimeno["India"], regimesize=regimesize["India"]
        )
        _region_metrics(
            nation_territory, pixel_seasia, t_idx=t_idx, hhi=hhi["SEAsia"],
            regimeno=regimeno["SEAsia"], regimesize=regimesize["SEAsia"]
        )
        _region_metrics(
            nation_territory, pixel_europax, t_idx=t_idx, hhi=hhi["EuropaX"],
            regimeno=regimeno["EuropaX"], regimesize=regimesize["EuropaX"]
        )

        _region_metrics(
            nation_territory, pixel_africa_et, t_idx=t_idx, hhi=np.zeros(t_max),
            regimeno=regimeno["AfricaET"], regimesize=regimesize["AfricaET"]
        )
        _region_metrics(
            nation_territory, pixel_africa_wt, t_idx=t_idx, hhi=np.zeros(t_max),
            regimeno=regimeno["AfricaWT"], regimesize=regimesize["AfricaWT"]
        )
        _region_metrics(
            nation_territory, pixel_africa_me, t_idx=t_idx, hhi=np.zeros(t_max),
            regimeno=regimeno["AfricaME"], regimesize=regimesize["AfricaME"]
        )
        _region_metrics(
            nation_territory, pixel_africa_nh, t_idx=t_idx, hhi=np.zeros(t_max),
            regimeno=regimeno["AfricaNH"], regimesize=regimesize["AfricaNH"]
        )
        _region_metrics(
            nation_territory, pixel_africa_sh, t_idx=t_idx, hhi=np.zeros(t_max),
            regimeno=regimeno["AfricaSH"], regimesize=regimesize["AfricaSH"]
        )
        _region_metrics(
            nation_territory, pixel_america_nh, t_idx=t_idx, hhi=np.zeros(t_max),
            regimeno=regimeno["AmericaNH"], regimesize=regimesize["AmericaNH"]
        )
        _region_metrics(
            nation_territory, pixel_america_sh, t_idx=t_idx, hhi=np.zeros(t_max),
            regimeno=regimeno["AmericaSH"], regimesize=regimesize["AmericaSH"]
        )
        _region_metrics(
            nation_territory, pixel_america_cl, t_idx=t_idx, hhi=np.zeros(t_max),
            regimeno=regimeno["AmericaCL"], regimesize=regimesize["AmericaCL"]
        )

        eu_max = float(max(territory_eu) if territory_eu else 0)
        cn_max = float(max(territory_cn) if territory_cn else 0)
        eu_max_id = int(np.argmax(territory_eu)) if territory_eu else 0
        cn_max_id = int(np.argmax(territory_cn)) if territory_cn else 0

        intersects["CN"][t_idx] = len([p for p in nation_territory[cn_max_id] if p in set(pixel_cn)])
        unions["CN"][t_idx] = len(set(pixel_cn).union(nation_territory[cn_max_id]))
        intersects["EU"][t_idx] = len([p for p in nation_territory[eu_max_id] if p in set(pixel_eu)])
        unions["EU"][t_idx] = len(set(pixel_eu).union(nation_territory[eu_max_id]))

        eu_max_orgpixel = nation_orgpixel[eu_max_id]
        cn_max_orgpixel = nation_orgpixel[cn_max_id]

        if timer in map_steps_set:
            map_snapshots[timer] = _build_pixel_belong(nation_territory, pixelno)

    return SimulationResult(
        metrics=hhi,
        regimesizes=regimesize,
        regimeno=regimeno,
        unions=unions,
        intersects=intersects,
        map_snapshots=map_snapshots,
        eu_max_orgpixel=eu_max_orgpixel,
        cn_max_orgpixel=cn_max_orgpixel,
        eu_max=eu_max,
        cn_max=cn_max,
    )


def _sum_cost_y(
    nation_id: int,
    target_pixel: int,
    nation_territory: list[list[int]],
    neighbors_land: list[list[int]],
    pixel_y: np.ndarray,
    sigma: float,
) -> float:
    connect = _pixel_land_connected(target_pixel, nation_territory[nation_id], neighbors_land)
    connect_y = sum(pixel_y[connect])
    unconnect_y = sum(pixel_y[nation_territory[nation_id]]) - connect_y
    return connect_y + sigma * unconnect_y
