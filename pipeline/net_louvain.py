"""
net_louvain.py  -  Louvain (Blondel et al., 2008)
==================================================
igraph `community_multilevel`, weighted, modularity at resolution gamma.
Stochastic (node order), so seed stability is reported; resolution sweep included.

Needs the networks from net_01_build.ipynb.
Output: <RES_DIR>/two_relations/louvain/
Run:    %run net_louvain.py      or   from net_louvain import run; run(Config())
"""
from net_common import Config, run_algorithm, seeded

ALGO = "louvain"


def detect(g, cfg, seed):
    seeded(seed)
    return g.community_multilevel(weights="weight", resolution=cfg.resolution).membership, {}


def run(cfg=None, networks=None):
    """networks: None/"all" = every network, or a list such as ["R1R2_weighted"]."""
    return run_algorithm(cfg or Config(), ALGO, detect, networks, stochastic=True, sweep=True)


if __name__ == "__main__":
    run()
