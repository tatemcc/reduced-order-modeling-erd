# import all configuration containers
from burgers_rom.config import  DataConfig, PODConfig, DerivConfig, SINDyConfig, RolloutConfig, RunConfig

CONFIG_FILE = "debug_config.yaml"

DataConfig(
    equation="burgers",
    structure="grid",
    resolution="low",
    split="train",
    n_trajectories=1,
    trajectory_ids=None,
    download_if_missing=False,
    lookback=1,
    rollout=1,
    squeeze_lookback_dim=True,
    time_stride=1,
    time_limit=None,
)
PODConfig()
DerivConfig()
SINDyConfig()
RolloutConfig()
RunConfig()

