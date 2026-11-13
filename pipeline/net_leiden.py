"""
net_leiden.py  -  Leiden (Traag, Waltman & van Eck, 2019)
==========================================================
igraph `community_leiden(objective_function="modularity")`, weighted,
iterated until no further improvement (n_iterations=-1). Unlike Louvain it
guarantees that every community is internally connected.
Stochastic, so seed stability is reported; resolution sweep included.

Needs the networks from net_01_build.ipynb.
Output: <RES_DIR>/two_relations/leiden/
Run:    %run net_leiden.py      or   from net_leiden import run; run(Config())
"""
from net_common import Config, run_algorithm, seeded

ALGO = "leiden"


def detect(g, cfg, seed):
    seeded(seed)
    m = g.community_leiden(objective_function="modularity", weights="weight",
                           resolution=cfg.resolution, n_iterations=-1).membership
    return m, {}


def run(cfg=None, networks=None):
    """networks: None/"all" = every network, or a list such as ["R1R2_weighted"]."""
    return run_algorithm(cfg or Config(), ALGO, detect, networks, stochastic=True, sweep=True)


if __name__ == "__main__":
    run()
