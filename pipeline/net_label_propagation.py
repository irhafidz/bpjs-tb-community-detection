"""
net_label_propagation.py  -  Label propagation (Raghavan, Albert & Kumara, 2007)
================================================================================
igraph `community_label_propagation`, weighted. Each node repeatedly takes the
label carrying the largest edge weight among its neighbours. No objective
function and no resolution parameter; results vary between runs, so seed
stability matters most for this algorithm.

Needs the networks from net_01_build.ipynb.
Output: <RES_DIR>/two_relations/label_propagation/
Run:    %run net_label_propagation.py
"""
from net_common import Config, run_algorithm, seeded

ALGO = "label_propagation"


def detect(g, cfg, seed):
    seeded(seed)
    return g.community_label_propagation(weights="weight").membership, {}


def run(cfg=None, networks=None):
    """networks: None/"all" = every network, or a list such as ["R1R2_weighted"]."""
    return run_algorithm(cfg or Config(), ALGO, detect, networks, stochastic=True, sweep=False)


if __name__ == "__main__":
    run()
