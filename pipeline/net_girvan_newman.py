"""
net_girvan_newman.py  -  Girvan-Newman (Girvan & Newman, 2002)
==============================================================
Divisive: repeatedly removes the edge with the highest edge betweenness
(igraph `community_edge_betweenness`). Betweenness uses edge lengths 1/weight,
so strong ties (household 1.0) are short and removed last. The dendrogram is
cut where weighted modularity is highest. Deterministic (no seed).

Cost: betweenness is recomputed after every removal, O(m^2 n). Benchmark:
144 nodes / 3,158 edges = 9 s; the largest R2 component (~3,000 nodes) would
take hours. So GN runs per connected component, only on components with at
most cfg.gn_max_edges edges, within cfg.gn_time_budget_s seconds per network
(smallest components first). Larger components stay unsplit; coverage is
reported as gn_nodes_covered_pct. A larger machine (Colab Pro, high RAM)
helps memory, not this single-threaded run time.

Needs the networks from net_01_build.ipynb.
Output: <RES_DIR>/two_relations/girvan_newman/
Run:    %run net_girvan_newman.py
"""
import time

import igraph as ig
import numpy as np

from net_common import Config, run_algorithm

ALGO = "girvan_newman"


def global_q_contribution(sg: ig.Graph, memb, total_w: float) -> float:
    """This component's share of the WHOLE network's weighted modularity:
    sum over groups of  W_in/W - (S/2W)^2  (W = total edge weight of the whole network).
    Components are independent, so maximising each share maximises the global Q."""
    memb = np.asarray(memb)
    el = np.array(sg.get_edgelist()).reshape(-1, 2)
    w = np.asarray(sg.es["weight"], float)
    same = memb[el[:, 0]] == memb[el[:, 1]]
    w_in = np.bincount(memb[el[same, 0]], weights=w[same], minlength=memb.max() + 1)
    s = np.bincount(memb, weights=np.asarray(sg.strength(weights="weight"), float), minlength=memb.max() + 1)
    return float((w_in / total_w - (s / (2 * total_w)) ** 2).sum())


def girvan_newman(g: ig.Graph, max_edges: int, budget_s: float = float("inf"),
                  cut: str = "component") -> tuple[list[int], dict]:
    """Girvan-Newman (2002) per connected component.
    Edge betweenness uses lengths 1/weight (strong ties = short).
    cut="component": each dendrogram is cut where the modularity OF THAT COMPONENT is highest
                     (finer partition: a small component is judged on its own).
    cut="global":    each dendrogram is cut where the component's contribution to the modularity
                     of the WHOLE network is highest (same objective as Louvain/Leiden).
    Components with more than max_edges edges are left as one community."""
    t0 = time.time()
    total_w = float(sum(g.es["weight"]))
    memb = np.zeros(g.vcount(), dtype=int)
    nxt, done_nodes, skipped, over_budget = 0, 0, [], 0
    # smallest components first, so the time budget covers as many as possible
    comps = sorted(g.connected_components(), key=len)
    for i, vs in enumerate(comps, 1):
        if i % 50 == 0 or i == len(comps):
            print(f"    component {i:,}/{len(comps):,} (size {len(vs):,}), {time.time() - t0:,.0f} s")
        sg = g.subgraph(vs)
        m = [0] * sg.vcount()
        too_late = time.time() - t0 > budget_s
        if sg.vcount() > 2 and sg.ecount() <= max_edges and too_late:
            over_budget += 1
            skipped.append((sg.vcount(), sg.ecount()))
        elif sg.vcount() > 2 and sg.ecount() <= max_edges:
            dend = sg.community_edge_betweenness(
                directed=False, weights=[1.0 / w for w in sg.es["weight"]])
            best_q = -np.inf
            for n in range(1, sg.vcount() + 1):
                cand = dend.as_clustering(n).membership
                q = (global_q_contribution(sg, cand, total_w) if cut == "global"
                     else sg.modularity(cand, weights="weight"))
                if q > best_q + 1e-12:
                    best_q, m = q, cand
            done_nodes += sg.vcount()
        elif sg.vcount() <= 2:
            done_nodes += sg.vcount()
        else:
            skipped.append((sg.vcount(), sg.ecount()))
        memb[vs] = np.asarray(m) + nxt
        nxt = memb[vs].max() + 1
    info = {"gn_nodes_covered_pct": 100 * done_nodes / max(g.vcount(), 1),
            "gn_components_not_split": len(skipped),
            "gn_largest_not_split": max(skipped)[0] if skipped else 0,
            "gn_not_split_time_budget": over_budget}
    return memb.tolist(), info


def detect(g, cfg, seed=None):
    print(f"  components with <= {cfg.gn_max_edges:,} edges, budget {cfg.gn_time_budget_s:,} s, cut = {cfg.gn_cut}")
    return girvan_newman(g, cfg.gn_max_edges, cfg.gn_time_budget_s, cfg.gn_cut)


def run(cfg=None, networks=None):
    """networks: None = cfg.gn_networks (main + household); "all" = every network;
    or a list such as ["R1R2_weighted"]."""
    cfg = cfg or Config()
    if networks is None:
        networks = "all" if "all" in cfg.gn_networks else list(cfg.gn_networks)
    # results of the global cut go to their own folder, so both versions are kept
    algo = ALGO if cfg.gn_cut == "component" else f"{ALGO}_globalcut"
    return run_algorithm(cfg, algo, detect, networks, stochastic=False, sweep=False)


if __name__ == "__main__":
    run()
