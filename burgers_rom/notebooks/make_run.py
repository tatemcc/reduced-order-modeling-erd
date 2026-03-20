from __future__ import annotations

from pathlib import Path
import json

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

from burgers_rom.config import RunConfig, DataConfig, PODConfig, DerivConfig, SINDyConfig, RolloutConfig, PlotConfig
from burgers_rom.run_pipeline import run

cfg = RunConfig(
    name="debug_single",
    seed=0,
    data=DataConfig(
        equation="burgers",
        structure="grid",
        resolution="low",
        split="train",
        n_trajectories=1,
        trajectory_ids=None,
        lookback=200,
        download_if_missing=False,
        time_stride=1,
        time_limit=None,
    ),
    pod=PODConfig(
        rank=5,
        energy_fraction=None,
        center=True,
    ),
    deriv=DerivConfig(),
    sindy=SINDyConfig(
        poly_order=2,
        include_bias=False,
        optimizer="stlsq",
        optimizer_kwargs={"threshold": 0.01},
        constrain_energy=False,
    ),  
    rollout=RolloutConfig(
        horizon_steps=190,
    ),
    plots=PlotConfig(
        sympy_labels=True,
    )
    # outputs_dir="outputs",
)

result = run(cfg, write_artifacts=True)
print(f"{result.run_id}, {result.rundir}, {result.aggregates}")
