from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
from zipfile import ZipFile

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry.base import BaseGeometry


@dataclass
class MapData:
    map_gdf: gpd.GeoDataFrame
    eadta: pd.DataFrame
    neighbors_land: list[list[int]]
    neighbors: list[list[int]]
    coast_fids: np.ndarray
    pixel_single: np.ndarray
    eurasia_bounds: tuple[float, float, float, float] | None


def _shapefile_from_zip(zip_path: Path) -> str:
    with ZipFile(zip_path) as zf:
        shp_files = [name for name in zf.namelist() if name.lower().endswith(".shp")]
    if not shp_files:
        raise FileNotFoundError(f"No .shp file found in {zip_path}")
    return f"zip://{zip_path}!{shp_files[0]}"


def _neighbors_from_geometries(
    geometries: Iterable[BaseGeometry],
    *,
    queen: bool,
) -> list[list[int]]:
    geoms = list(geometries)
    sindex = gpd.GeoSeries(geoms).sindex
    neighbors: list[list[int]] = [[] for _ in range(len(geoms))]
    for idx, geom in enumerate(geoms):
        if geom is None or geom.is_empty:
            continue
        candidate_ids = list(sindex.intersection(geom.bounds))
        for cand in candidate_ids:
            if cand == idx:
                continue
            other = geoms[cand]
            if other is None or other.is_empty:
                continue
            if not geom.touches(other):
                continue
            if not queen:
                shared = geom.boundary.intersection(other.boundary)
                if shared.is_empty or shared.length == 0:
                    continue
            neighbors[idx].append(cand)
        neighbors[idx] = sorted(set(neighbors[idx]))
    return neighbors


def load_map_data(map_dir: Path, max_sea_dist: int) -> MapData:
    map_dir = map_dir.expanduser()
    world_geom_path = _shapefile_from_zip(map_dir / "World5_Nodta.zip")
    world_data_path = _shapefile_from_zip(map_dir / "World5.zip")

    map_gdf = gpd.read_file(world_geom_path)
    eadta_gdf = gpd.read_file(world_data_path)
    eadta = eadta_gdf.drop(columns=["geometry"], errors="ignore")

    eurasia_bounds = None
    eurasia_zip = map_dir / "Euroasia_Dis.zip"
    if eurasia_zip.exists():
        eurasia_path = _shapefile_from_zip(eurasia_zip)
        eurasia_gdf = gpd.read_file(eurasia_path)
        eurasia_bounds = tuple(float(x) for x in eurasia_gdf.total_bounds)

    if len(map_gdf) != len(eadta):
        raise ValueError("Map geometry and attribute data row counts do not match.")

    neighbors_land = _neighbors_from_geometries(map_gdf.geometry, queen=False)

    coastnear_table = pd.read_csv(map_dir / "DistMatrix_World.csv")
    if "Unnamed: 0" in coastnear_table.columns:
        coastnear_table = coastnear_table.drop(columns=["Unnamed: 0"])
    coastnear_table = coastnear_table[coastnear_table["Dist"] <= max_sea_dist].copy()
    coastnear_table[["coastid", "NbList"]] = coastnear_table[["coastid", "NbList"]] - 1

    coastnear_list = (
        coastnear_table.groupby("coastid")["NbList"].apply(list).to_dict()
    )
    coast_fids = np.array(sorted(coastnear_list.keys()), dtype=int)

    neighbors = [list(nbs) for nbs in neighbors_land]
    for coast_fid, sea_neighbors in coastnear_list.items():
        merged = set(neighbors[coast_fid])
        merged.update(sea_neighbors)
        neighbors[coast_fid] = sorted(merged)

    pixel_single = np.array(
        [idx for idx, nbs in enumerate(neighbors) if len(nbs) == 0],
        dtype=int,
    )

    return MapData(
        map_gdf=map_gdf,
        eadta=eadta,
        neighbors_land=neighbors_land,
        neighbors=neighbors,
        coast_fids=coast_fids,
        pixel_single=pixel_single,
        eurasia_bounds=eurasia_bounds,
    )
