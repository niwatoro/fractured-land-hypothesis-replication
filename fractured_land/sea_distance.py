from __future__ import annotations

from collections import deque
from typing import Iterable


def sea_distance(
    start: int,
    territory: Iterable[int],
    sea_neighbors: list[list[int]],
    coast_ids: set[int],
) -> tuple[list[int], list[int]]:
    territory_list = list(territory)
    territory_set = set(territory_list)
    graph = [list(territory_set.intersection(sea_neighbors[pix])) for pix in territory_list]

    distances = [float("inf")] * len(graph)
    visited = [False] * len(graph)
    start_idx = territory_list.index(start)
    distances[start_idx] = 0

    while True:
        shortest_distance = float("inf")
        shortest_index = -1
        for idx, dist in enumerate(distances):
            if dist < shortest_distance and not visited[idx]:
                shortest_distance = dist
                shortest_index = idx

        if shortest_index == -1:
            other_coast = coast_ids.difference({start}).intersection(territory_set)
            out_dist = []
            out_nb = []
            for node in other_coast:
                node_idx = territory_list.index(node)
                neighbors = graph[node_idx]
                if neighbors:
                    min_dist = min(distances[territory_list.index(n)] for n in neighbors)
                else:
                    min_dist = float("inf")
                out_dist.append(min_dist + 1)
                out_nb.append(node)
            return out_dist, out_nb

        if graph[shortest_index]:
            for node in graph[shortest_index]:
                node_idx = territory_list.index(node)
                if distances[node_idx] > distances[shortest_index] + 1:
                    distances[node_idx] = distances[shortest_index] + 1
            visited[shortest_index] = True
        else:
            visited[shortest_index] = True


def compute_distance_matrix(
    coast_ids: list[int],
    coast_neighbor_lists: list[list[int]],
    sea_neighbors: list[list[int]],
) -> list[tuple[int, int, int]]:
    coast_id_set = set(coast_ids)
    results = []
    for start, territory in zip(coast_ids, coast_neighbor_lists):
        dist, nb_ids = sea_distance(start, territory, sea_neighbors, coast_id_set)
        for nb, d in zip(nb_ids, dist):
            results.append((start, nb, d))
    return results
