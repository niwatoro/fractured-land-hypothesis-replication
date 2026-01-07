from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ScenarioSettings:
    theta_cold: float
    theta_hot: float
    ydiscount_cold: float
    ydiscount_hot: float
    theta_rugged_90th: float
    sigma: float
    alpha_sea: float
    yield_data: str
    bc1000_ratio: str
    uniform_y: int
    conflict_mech: str
    river: int
    theta_steppe: float
    psi_steppe: float
    y_steppe: float | None


@dataclass(frozen=True)
class ModelConfig:
    scenario_name: str
    scenario: ScenarioSettings
    beta: float = 0.00001
    shockprob_regime: float = 0.0
    shockprob_general: float = 0.0
    max_sea_dist: int = 6
    start_with_regime: int = 0
    ymin: float = 0.001
    unrest_min: float = 0.000001
    t_max: int = 500
    draw_map: int = 0
    t_plot: int = 100
    draw_eurasia: int = 1
    num_cores: int = 3
    try_no: int = 3
    africa_start: int = 0
    america_start: int = 0
    australia_start: int = 0
    japan_start: int = 0
    maize1_start: int = 0
    maize2_start: int = 0
    maize3_start: int = 0
    y_steppeeast: float | None = None


def scenario_settings(name: str) -> ScenarioSettings:
    scenarios = {
        "Preferred": ScenarioSettings(
            theta_cold=2,
            theta_hot=2,
            ydiscount_cold=0.8,
            ydiscount_hot=0.8,
            theta_rugged_90th=2,
            sigma=0.333,
            alpha_sea=0.1,
            yield_data="YGAEZ",
            bc1000_ratio="startloess",
            uniform_y=0,
            conflict_mech="neighborY",
            river=1,
            theta_steppe=-1,
            psi_steppe=3,
            y_steppe=1116.475,
        ),
        "Adj. for hot and cold climate": ScenarioSettings(
            theta_cold=2,
            theta_hot=2,
            ydiscount_cold=0.8,
            ydiscount_hot=0.8,
            theta_rugged_90th=2,
            sigma=0.333,
            alpha_sea=0.1,
            yield_data="YGAEZ",
            bc1000_ratio="startDSMW",
            uniform_y=0,
            conflict_mech="neighborY",
            river=0,
            theta_steppe=0,
            psi_steppe=1,
            y_steppe=None,
        ),
        "Nomadic Pastoralism": ScenarioSettings(
            theta_cold=2,
            theta_hot=2,
            ydiscount_cold=0,
            ydiscount_hot=0,
            theta_rugged_90th=2,
            sigma=0.333,
            alpha_sea=0.1,
            yield_data="YGAEZ",
            bc1000_ratio="startDSMW",
            uniform_y=0,
            conflict_mech="neighborY",
            river=0,
            theta_steppe=0,
            psi_steppe=1,
            y_steppe=1116.475,
        ),
        "Silk Road+Steppe Expansionism": ScenarioSettings(
            theta_cold=2,
            theta_hot=2,
            ydiscount_cold=0,
            ydiscount_hot=0,
            theta_rugged_90th=2,
            sigma=0.333,
            alpha_sea=0.1,
            yield_data="YGAEZ",
            bc1000_ratio="startDSMW",
            uniform_y=0,
            conflict_mech="neighborY",
            river=0,
            theta_steppe=-1,
            psi_steppe=3,
            y_steppe=None,
        ),
        "Loess": ScenarioSettings(
            theta_cold=2,
            theta_hot=2,
            ydiscount_cold=0,
            ydiscount_hot=0,
            theta_rugged_90th=2,
            sigma=0.333,
            alpha_sea=0.1,
            yield_data="YGAEZ",
            bc1000_ratio="startloess",
            uniform_y=0,
            conflict_mech="neighborY",
            river=0,
            theta_steppe=0,
            psi_steppe=1,
            y_steppe=None,
        ),
        "Rivers": ScenarioSettings(
            theta_cold=2,
            theta_hot=2,
            ydiscount_cold=0,
            ydiscount_hot=0,
            theta_rugged_90th=2,
            sigma=0.333,
            alpha_sea=0.1,
            yield_data="YGAEZ",
            bc1000_ratio="startDSMW",
            uniform_y=0,
            conflict_mech="neighborY",
            river=1,
            theta_steppe=0,
            psi_steppe=1,
            y_steppe=None,
        ),
        "Alt. Conflict Mech: Random": ScenarioSettings(
            theta_cold=2,
            theta_hot=2,
            ydiscount_cold=0,
            ydiscount_hot=0,
            theta_rugged_90th=2,
            sigma=0.333,
            alpha_sea=0.1,
            yield_data="YGAEZ",
            bc1000_ratio="startDSMW",
            uniform_y=0,
            conflict_mech="random",
            river=0,
            theta_steppe=0,
            psi_steppe=1,
            y_steppe=None,
        ),
        "Alt. Conflict Mech: Expected Gain": ScenarioSettings(
            theta_cold=2,
            theta_hot=2,
            ydiscount_cold=0,
            ydiscount_hot=0,
            theta_rugged_90th=2,
            sigma=0.333,
            alpha_sea=0.1,
            yield_data="YGAEZ",
            bc1000_ratio="startDSMW",
            uniform_y=0,
            conflict_mech="expgain",
            river=0,
            theta_steppe=0,
            psi_steppe=1,
            y_steppe=None,
        ),
        "Alt. r:KK10": ScenarioSettings(
            theta_cold=2,
            theta_hot=2,
            ydiscount_cold=0,
            ydiscount_hot=0,
            theta_rugged_90th=2,
            sigma=0.333,
            alpha_sea=0.1,
            yield_data="YGAEZ",
            bc1000_ratio="startKK10",
            uniform_y=0,
            conflict_mech="neighborY",
            river=0,
            theta_steppe=0,
            psi_steppe=1,
            y_steppe=None,
        ),
        "Alt. Y CSI": ScenarioSettings(
            theta_cold=2,
            theta_hot=2,
            ydiscount_cold=0,
            ydiscount_hot=0,
            theta_rugged_90th=2,
            sigma=0.333,
            alpha_sea=0.1,
            yield_data="YCSI",
            bc1000_ratio="startDSMW",
            uniform_y=0,
            conflict_mech="neighborY",
            river=0,
            theta_steppe=0,
            psi_steppe=1,
            y_steppe=None,
        ),
        "No Obstacles+Uniform Resource": ScenarioSettings(
            theta_cold=0,
            theta_hot=0,
            ydiscount_cold=0,
            ydiscount_hot=0,
            theta_rugged_90th=0,
            sigma=1,
            alpha_sea=1,
            yield_data="",
            bc1000_ratio="",
            uniform_y=1,
            conflict_mech="neighborY",
            river=0,
            theta_steppe=0,
            psi_steppe=1,
            y_steppe=None,
        ),
        "Uniform Resource": ScenarioSettings(
            theta_cold=2,
            theta_hot=2,
            ydiscount_cold=0,
            ydiscount_hot=0,
            theta_rugged_90th=2,
            sigma=0.333,
            alpha_sea=0.1,
            yield_data="",
            bc1000_ratio="",
            uniform_y=1,
            conflict_mech="neighborY",
            river=0,
            theta_steppe=0,
            psi_steppe=1,
            y_steppe=None,
        ),
        "No Obstacles": ScenarioSettings(
            theta_cold=0,
            theta_hot=0,
            ydiscount_cold=0,
            ydiscount_hot=0,
            theta_rugged_90th=0,
            sigma=1,
            alpha_sea=1,
            yield_data="YGAEZ",
            bc1000_ratio="startDSMW",
            uniform_y=0,
            conflict_mech="neighborY",
            river=0,
            theta_steppe=0,
            psi_steppe=1,
            y_steppe=None,
        ),
        "No Climatic Obstacles": ScenarioSettings(
            theta_cold=0,
            theta_hot=0,
            ydiscount_cold=0,
            ydiscount_hot=0,
            theta_rugged_90th=2,
            sigma=0.333,
            alpha_sea=0.1,
            yield_data="YGAEZ",
            bc1000_ratio="startDSMW",
            uniform_y=0,
            conflict_mech="neighborY",
            river=0,
            theta_steppe=0,
            psi_steppe=1,
            y_steppe=None,
        ),
        "Baseline": ScenarioSettings(
            theta_cold=2,
            theta_hot=2,
            ydiscount_cold=0,
            ydiscount_hot=0,
            theta_rugged_90th=2,
            sigma=0.333,
            alpha_sea=0.1,
            yield_data="YGAEZ",
            bc1000_ratio="startDSMW",
            uniform_y=0,
            conflict_mech="neighborY",
            river=0,
            theta_steppe=0,
            psi_steppe=1,
            y_steppe=None,
        ),
    }
    if name not in scenarios:
        raise ValueError(f"Unknown scenario: {name}")
    return scenarios[name]


def build_config(scenario_name: str) -> ModelConfig:
    return ModelConfig(scenario_name=scenario_name, scenario=scenario_settings(scenario_name))
