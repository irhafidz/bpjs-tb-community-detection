"""
net_common.py  -  shared code for the two-relation TB patient networks
=======================================================================
Builds the patient-to-patient networks once and gives every community
detection algorithm the same graph. Algorithms live in their own files:

    net_louvain.py            Louvain            (Blondel et al., 2008)
    net_leiden.py             Leiden             (Traag et al., 2019)
    net_label_propagation.py  Label propagation  (Raghavan et al., 2007)
    net_girvan_newman.py      Girvan-Newman      (Girvan & Newman, 2002)

Relations (Prof Iin, 11 Nov 2026)
  R1 household   same PSTV02 (family number)                         weight 1.0
  R2 clinic_dx   FKTP visit at the same clinic proxy (district FKP06,
                 optionally + type FKP08 / ownership FKP07) with the same
                 3-char ICD-10 (FKP14A), clinic and code from one visit      weight 0.5
  MAIN network   R1 + R2 weighted (a pair in both = 1.5)

Networks written to <RES_DIR>/two_relations/networks/:
  R1_household, R2_clinic_dx, R2_sens_Z20_Z11, R1R2_weighted (MAIN), R1R2_weighted_sens

The data have no facility ID, so "same klinik pratama" is a proxy.
Put this file in bpjs-tb-graph/pipeline/ next to common.py and config.py.
"""
from __future__ import annotations

import csv
import glob
import json
import os
import random
import re
import time
from dataclasses import dataclass, field, asdict, replace

import igraph as ig
import numpy as np
import pandas as pd
import scipy.sparse as sp


# ----------------------------------------------------------------------------
# 0. CONFIG
# ----------------------------------------------------------------------------
@dataclass
class Config:
    # "common" = read through the repo's common.load_tables() (same data and cache as 01-12)
    # "files"  = read CSV files from data_dir with the glob patterns in `tables`
    source: str = "common"
    data_dir: str = "data/csv"
    out_dir: str | None = None   # None = <RES_DIR from config.py>/two_relations

    # glob pattern per table, relative to data_dir (all matching files are stacked)
    # >>> EDIT THESE to match your file names <<<
    tables: dict = field(default_factory=lambda: {
        "membership": "*kepesertaan*.csv",
        "fktp":       "*fktp*kapitasi*.csv",
        "noncap":     "*non*kapitasi*.csv",     # optional
        "fkrtl":      "*fkrtl*.csv",            # optional (only for TB cohort)
        "fkrtl_dx":   "*diagnosis*sekunder*.csv",  # optional (only for TB cohort)
    })

    year: int | None = 2024        # keep visits in this year only (None = all years)
    tb_codes: tuple = ("A15", "A16", "A17", "A18", "A19")

    # Which patients are nodes?
    #   "tb"  = patients with any A15-A19 code (any table)  [default]
    #   "all" = every member in Membership
    cohort: str = "tb"

    # R1: also add the sampled household members of TB patients as nodes
    #     (contact-tracing view). False = only cohort patients are nodes.
    r1_include_contacts: bool = True

    # R2 options
    r2_sources: tuple = ("fktp",)  # add "noncap" to also use non-capitation FKTP visits
    r2_exclude_z: bool = True      # Z00-Z99 are neither diagnosis nor symptom
    r2_icd_filter: tuple | None = None  # e.g. ("A15","A16") to restrict codes; None = all
    # facility proxy besides district (no facility ID in the data):
    #   () = district only; ("FKP08",) = + FKTP type; ("FKP08","FKP07") = + type + ownership
    r2_facility_cols: tuple = ()

    # relation weights (Prof Iin): household = 1, same clinic + same symptom = 0.5
    w_household: float = 1.0
    w_clinic: float = 0.5
    # "binary": a clinic tie is worth w_clinic however many codes are shared
    # "count" : w_clinic x number of shared (facility, ICD) keys
    r2_weight: str = "binary"
    # Z codes kept even when r2_exclude_z=True, matched on the FULL code (FKP15/PNK14).
    # Main analysis: () = all Z codes dropped.
    r2_keep_z_full: tuple = ()
    # Sensitivity analysis: rerun R2 keeping only these Z codes (prefix match on FKP15).
    #   Z20 = contact with / exposure to communicable diseases
    #   Z11 = special screening for infectious and parasitic diseases
    # In this dataset FKP15 holds 3-char code + name (no sub-codes), so Z20.1 / Z11.1
    # cannot be identified; the 3-char codes are used instead.
    run_sensitivity: bool = True
    r2_sensitivity_keep_z: tuple = ("Z20", "Z11")
    max_edges: int = 50_000_000    # safety stop for clique explosion (upper bound)

    # community detection (shared by all algorithm files)
    seed: int = 42
    resolution: float = 1.0
    n_stability_runs: int = 20
    resolution_grid: tuple = (0.5, 0.75, 1.0, 1.5, 2.0)   # Louvain and Leiden only

    # algorithms compared in net_06_compare (each one is run by its own file)
    algorithms: tuple = ("louvain", "leiden", "label_propagation", "girvan_newman", "girvan_newman_globalcut")
    # Girvan-Newman recomputes edge betweenness after every removal (O(m^2 n)),
    # so it is run per connected component and only on components with at most
    # this many edges; larger components are left unsplit and reported.
    gn_max_edges: int = 3000
    # GN is slow, so it runs only on these networks ("all" = every network) and
    # stops splitting new components once this many seconds have passed.
    gn_networks: tuple = ("R1_household", "R1R2_weighted")
    gn_time_budget_s: int = 1800
    # "component": cut each GN dendrogram at that component's own modularity maximum;
    # "global": at the maximum of the whole network's modularity (output folder girvan_newman_globalcut/)
    gn_cut: str = "component"




# columns used per table (codebook Lampiran A-E)
COLS = {
    "membership": ["PSTV01", "PSTV02", "PSTV04", "PSTV05", "PSTV09", "PSTV09_NEW", "PSTV10", "PSTV15",
                   "PSTV17", "PSTV18"],
    "fktp":       ["PSTV01", "FKP02", "FKP03", "FKP05", "FKP05_NEW", "FKP06", "DISTRICT",
                   "FKP07", "FKP08", "FKP09", "FKP14A", "FKP15"],
    "noncap":     ["PSTV01", "PNK02", "PNK03", "PNK06", "PNK06_NEW", "PNK07", "PNK08", "PNK09",
                   "PNK10", "PNK13A", "PNK14"],
    "fkrtl":      ["PSTV01", "FKL02", "FKL03", "FKL15A", "FKL17A"],
    "fkrtl_dx":   ["FKL02", "FKL24A"],
}
REQUIRED = {"membership", "fktp"}


# ----------------------------------------------------------------------------
# 1. LOAD + CLEAN
# ----------------------------------------------------------------------------
def _sniff_sep(path: str) -> str:
    with open(path, "r", encoding="utf-8-sig", errors="replace") as f:
        sample = f.read(65536)
    try:
        return csv.Sniffer().sniff(sample, delimiters=",;|\t").delimiter
    except csv.Error:
        return ","


def load_table(cfg: Config, name: str) -> pd.DataFrame | None:
    pattern = os.path.join(cfg.data_dir, cfg.tables[name])
    files = sorted(glob.glob(pattern))
    if not files:
        if name in REQUIRED:
            raise FileNotFoundError(f"[{name}] no file matches {pattern}")
        print(f"[{name}] no file matches {pattern} -> skipped")
        return None
    want = set(COLS[name])
    parts = []
    for fp in files:
        df = pd.read_csv(fp, sep=_sniff_sep(fp), dtype=str, encoding="utf-8-sig",
                         usecols=lambda c: c.strip().upper() in want,
                         low_memory=False)
        df.columns = [c.strip().upper() for c in df.columns]
        parts.append(df)
    df = pd.concat(parts, ignore_index=True)
    missing = want - set(df.columns)
    if missing:
        print(f"[{name}] missing columns: {sorted(missing)}")
    for c in df.columns:
        df[c] = df[c].str.strip().replace({"": np.nan, "nan": np.nan, "NaN": np.nan})
    print(f"[{name}] {len(files)} file(s), {len(df):,} rows")
    return df


def clean_icd3(s: pd.Series) -> pd.Series:
    s = s.astype("string").str.upper().str.replace(r"[^A-Z0-9]", "", regex=True).str[:3]
    return s.where(s.str.fullmatch(r"[A-Z]\d\d", na=False))


def clean_district(s: pd.Series) -> pd.Series:
    """BPS code if the value holds >= 4 digits (handles '3273', '3273.0');
    otherwise the upper-cased text as given (e.g. 'KOTA BANDUNG')."""
    raw = s.astype("string").str.strip().str.upper()
    dig = raw.str.replace(r"\.0+$", "", regex=True).str.replace(r"\D", "", regex=True)
    out = dig.where(dig.str.len() >= 4)
    return out.fillna(raw.where(raw.str.len() > 0))


def district_key(df: pd.DataFrame, prov_cols: tuple, dist_col: str) -> pd.Series:
    """'PROVINCE | DISTRICT' key. District names repeat across provinces (15 names
    in this sample), so the name alone would merge different districts.
    Uses the harmonised DISTRICT column from common.py when present (FKTP),
    otherwise province (*_NEW = 38 provinces, else old code) + district."""
    if dist_col == "FKP06" and "DISTRICT" in df and df["DISTRICT"].notna().any():
        return df["DISTRICT"].astype("string").str.strip().str.upper().replace("", pd.NA)
    d = clean_district(df[dist_col])
    for pc in prov_cols:
        if pc in df and df[pc].notna().any():
            p = df[pc].astype("string").str.strip().str.upper()
            return (p + " | " + d).where(d.notna() & p.notna(), d)
    return d


def year_of(s: pd.Series) -> pd.Series:
    d = pd.to_datetime(s, errors="coerce", format="ISO8601")
    if d.isna().mean() > 0.5:
        d = pd.to_datetime(s, errors="coerce", dayfirst=True)
    return d.dt.year


def _keep_year(df: pd.DataFrame, date_col: str, cfg: Config) -> pd.DataFrame:
    if cfg.year is None or date_col not in df:
        return df
    return df[year_of(df[date_col]) == cfg.year]


COMMON_NAMES = {"member": "membership", "fktp": "fktp", "fkrtl": "fkrtl",
                "secdx": "fkrtl_dx", "noncap": "noncap"}


def _to_text(df: pd.DataFrame, name: str) -> pd.DataFrame | None:
    if df is None or len(df) == 0:
        print(f"[{name}] empty -> skipped")
        return None
    df = df.copy()
    df.columns = [str(c).strip().upper() for c in df.columns]
    keep = [c for c in COLS[name] if c in df.columns]
    missing = sorted(set(COLS[name]) - set(keep))
    if missing:
        print(f"[{name}] missing columns: {missing}")
    df = df[keep]
    for c in keep:
        df[c] = (df[c].astype("string").str.strip()
                 .replace({"": pd.NA, "nan": pd.NA, "NaN": pd.NA, "NaT": pd.NA, "<NA>": pd.NA}))
    print(f"[{name}] {len(df):,} rows")
    return df


def load_from_common() -> dict:
    """Use the repo's common.load_tables() (cache built by 01_load_and_profile.py)."""
    import common
    T = common.load_tables(harmonise=True)
    T = {k: v for k, v in T.items() if isinstance(v, pd.DataFrame)}   # skip info dicts etc.
    return {new: _to_text(T.get(old), new) for old, new in COMMON_NAMES.items()}


def default_out_dir() -> str:
    """<RES_DIR from config.py>/two_relations (falls back to ./results/...)."""
    import importlib
    for mod in ("config", "common"):
        try:
            m = importlib.import_module(mod)
        except ImportError:
            continue
        for a in ("RES_DIR", "RESULTS", "RES", "RESULTS_DIR"):
            if hasattr(m, a):
                return os.path.join(str(getattr(m, a)), "two_relations")
    return os.path.join("results", "two_relations")


def load_all(cfg: Config) -> dict:
    if cfg.out_dir is None:
        cfg.out_dir = default_out_dir()
    if cfg.source == "common":
        t = load_from_common()
        need = {"membership": ["PSTV01", "PSTV02"], "fktp": ["PSTV01", "FKP06", "FKP14A"]}
        for k, cols in need.items():
            if t.get(k) is None or not set(cols) <= set(t[k].columns):
                raise KeyError(f"[{k}] needs columns {cols}; check the column names "
                               "returned by common.load_tables()")
    else:
        t = {k: load_table(cfg, k) for k in cfg.tables}
    m = t["membership"]
    m = m.dropna(subset=["PSTV01"]).drop_duplicates("PSTV01")
    t["membership"] = m
    return t


# ----------------------------------------------------------------------------
# 2. COHORT (node set)
# ----------------------------------------------------------------------------
def tb_patients(t: dict, cfg: Config) -> set:
    tb = set(cfg.tb_codes)
    ids = set()
    if t.get("fktp") is not None:
        f = _keep_year(t["fktp"], "FKP03", cfg)
        ids |= set(f.loc[clean_icd3(f["FKP14A"]).isin(tb), "PSTV01"].dropna())
    if t.get("noncap") is not None:
        f = _keep_year(t["noncap"], "PNK03", cfg)
        ids |= set(f.loc[clean_icd3(f["PNK13A"]).isin(tb), "PSTV01"].dropna())
    if t.get("fkrtl") is not None:
        h = _keep_year(t["fkrtl"], "FKL03", cfg)
        for c in ("FKL15A", "FKL17A"):
            if c in h:
                ids |= set(h.loc[clean_icd3(h[c]).isin(tb), "PSTV01"].dropna())
        if t.get("fkrtl_dx") is not None:
            d = t["fkrtl_dx"]
            vis = set(d.loc[clean_icd3(d["FKL24A"]).isin(tb), "FKL02"].dropna())
            ids |= set(h.loc[h["FKL02"].isin(vis), "PSTV01"].dropna())
    return ids


def build_patients(t: dict, cfg: Config) -> pd.DataFrame:
    m = t["membership"].copy()
    tb = tb_patients(t, cfg)
    m["is_tb"] = m["PSTV01"].isin(tb)
    m["died"] = m.get("PSTV18", pd.Series(index=m.index)).notna() | \
        m.get("PSTV17", pd.Series("", index=m.index)).fillna("").str.upper().str.contains("MENINGGAL")
    m["district_res"] = (district_key(m, ("PSTV09_NEW", "PSTV09"), "PSTV10")
                         if "PSTV10" in m else pd.NA)
    m["in_cohort"] = m["is_tb"] if cfg.cohort == "tb" else True
    print(f"members {len(m):,} | TB-coded {m['is_tb'].sum():,} | cohort {m['in_cohort'].sum():,}")
    if cfg.cohort == "tb" and len(tb - set(m["PSTV01"])):
        print(f"  note: {len(tb - set(m['PSTV01'])):,} TB patients not in Membership (dropped)")
    return m.reset_index(drop=True)


# ----------------------------------------------------------------------------
# 3. EDGES FROM SHARED KEYS  (patient x key incidence -> B B^T)
# ----------------------------------------------------------------------------
def edges_from_shared_key(pk: pd.DataFrame, node_col: str, key_col: str,
                          max_edges: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    """pk: one row per (node, key). Returns edges (src, dst, weight) and group sizes.
    Two nodes are tied if they share >= 1 key; weight = number of shared keys."""
    pk = pk[[node_col, key_col]].dropna().drop_duplicates()
    groups = pk.groupby(key_col)[node_col].size().rename("n_patients").reset_index()
    groups["pairs"] = groups["n_patients"] * (groups["n_patients"] - 1) // 2
    groups = groups.sort_values("n_patients", ascending=False).reset_index(drop=True)
    bound = int(groups["pairs"].sum())
    print(f"  keys {len(groups):,} | upper bound on pairs {bound:,}")
    if bound > max_edges:
        raise MemoryError(
            f"{bound:,} candidate pairs > max_edges={max_edges:,}. Largest keys:\n"
            f"{groups.head(10).to_string()}\n"
            "Restrict r2_icd_filter, keep r2_exclude_z=True, or raise max_edges.")
    if pk.empty:
        return pd.DataFrame(columns=["src", "dst", "weight"]), groups
    nodes = pd.Index(pk[node_col].unique())
    keys = pd.Index(pk[key_col].unique())
    B = sp.csr_matrix((np.ones(len(pk), dtype=np.int32),
                       (nodes.get_indexer(pk[node_col]), keys.get_indexer(pk[key_col]))),
                      shape=(len(nodes), len(keys)))
    W = sp.triu(B @ B.T, k=1).tocoo()
    e = pd.DataFrame({"src": nodes[W.row], "dst": nodes[W.col], "weight": W.data})
    return e, groups


def build_r1_household(pat: pd.DataFrame, cfg: Config):
    """R1: same PSTV02 (family number)."""
    if cfg.cohort == "tb" and cfg.r1_include_contacts:
        hh_tb = set(pat.loc[pat["in_cohort"], "PSTV02"].dropna())
        nodes = pat[pat["in_cohort"] | pat["PSTV02"].isin(hh_tb)]
    else:
        nodes = pat[pat["in_cohort"]]
    nodes = nodes.dropna(subset=["PSTV02"])
    print(f"R1 household: {len(nodes):,} candidate nodes, "
          f"{nodes['PSTV02'].nunique():,} households")
    e, g = edges_from_shared_key(nodes, "PSTV01", "PSTV02", cfg.max_edges)
    e["weight"] = cfg.w_household * (e["weight"] > 0)   # one household per patient -> 1 x w
    e["relation"] = "household"
    g = g.rename(columns={"PSTV02": "household"})
    return nodes.copy(), e, g


def r2_patient_keys(t: dict, pat: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    """One row per (patient, district|ICD) from clinic visits."""
    src = {"fktp": ("FKP03", "FKP06", "FKP14A", "FKP15"),
           "noncap": ("PNK03", "PNK07", "PNK13A", "PNK14")}
    parts = []
    for name in cfg.r2_sources:
        df = t.get(name)
        if df is None:
            continue
        dcol, distcol, icdcol, fullcol = src[name]
        print(f"  [{name}] rows {len(df):,}", end="")
        df = _keep_year(df, dcol, cfg)
        print(f" | year {cfg.year}: {len(df):,}", end="")
        if distcol not in df or icdcol not in df:
            print(f" | missing {distcol} or {icdcol} -> skipped")
            continue
        prov = ("FKP05_NEW", "FKP05") if name == "fktp" else ("PNK06_NEW", "PNK06")
        dist = district_key(df, prov, distcol)
        print(f" | district non-missing: {dist.notna().sum():,}"
              f" | ICD-10 valid: {clean_icd3(df[icdcol]).notna().sum():,}")
        full = (df[fullcol].astype("string").str.upper().str.replace(r"[^A-Z0-9]", "", regex=True)
                if fullcol in df else pd.Series(pd.NA, index=df.index, dtype="string"))
        fac = dist
        nc_map = {"FKP07": "PNK08", "FKP08": "PNK09", "FKP09": "PNK10"}
        for c in cfg.r2_facility_cols:
            cc = c if name == "fktp" else nc_map.get(c, c)
            if cc in df:
                fac = fac + "/" + df[cc].astype("string").str.strip().str.upper().fillna("NA")
            else:
                print(f"  [{name}] facility column {cc} missing -> not used")
        parts.append(pd.DataFrame({"PSTV01": df["PSTV01"],
                                   "district": dist, "facility": fac,
                                   "icd3": clean_icd3(df[icdcol]), "icd_full": full}))
    if not parts:
        return pd.DataFrame(columns=["PSTV01", "district", "icd3", "key"])
    v = pd.concat(parts, ignore_index=True).dropna(subset=["PSTV01", "district", "icd3"])
    print(f"  rows with patient + district + ICD: {len(v):,}", end="")
    v = v[v["PSTV01"].isin(set(pat.loc[pat["in_cohort"], "PSTV01"]))]
    print(f" | in TB cohort: {len(v):,}", end="")
    isz = v["icd3"].str.startswith("Z")
    keepz = pd.Series(False, index=v.index)
    for code in cfg.r2_keep_z_full:          # prefix match: "Z20" matches Z20, Z201, Z20CONTACT...
        norm = re.sub(r"[^A-Z0-9]", "", code.upper())
        hit = isz & v["icd_full"].fillna("").str.startswith(norm)
        v.loc[hit, "icd3"] = code             # key uses the full code, e.g. 3273|Z20.1
        keepz |= hit
    if cfg.r2_keep_z_full:
        print(f"  kept Z rows: {dict(v.loc[keepz, 'icd3'].value_counts())}")
    if cfg.r2_exclude_z:
        v = v[~isz | keepz]
    print(f" | after Z rule: {len(v):,}")
    if cfg.r2_icd_filter:
        v = v[v["icd3"].isin(set(cfg.r2_icd_filter))]
    v["key"] = v["facility"] + "|" + v["icd3"]
    v = v.drop(columns=["icd_full", "facility"])
    return v.drop_duplicates(["PSTV01", "key"]).reset_index(drop=True)


def build_r2_clinic_dx(t: dict, pat: pd.DataFrame, cfg: Config):
    """R2: same clinic district AND same 3-char ICD-10 (on a clinic visit)."""
    pk = r2_patient_keys(t, pat, cfg)
    print(f"R2 clinic_dx: {pk['PSTV01'].nunique():,} cohort patients with clinic visits, "
          f"{len(pk):,} patient-(district,ICD) pairs")
    e, g = edges_from_shared_key(pk, "PSTV01", "key", cfg.max_edges)
    e = e.rename(columns={"weight": "n_shared"})
    e["weight"] = cfg.w_clinic * (1 if cfg.r2_weight == "binary" else e["n_shared"])
    e["relation"] = "clinic_dx"
    if len(g):
        g[["facility", "icd3"]] = g["key"].str.rsplit("|", n=1, expand=True)
        g["district"] = g["facility"].str.split("/", n=1).str[0]
    else:
        for c in ("facility", "district", "icd3"):
            g[c] = pd.Series(dtype=str)
        print("  WARNING: R2 has no keys -> check the row counts above")
    nodes = pat[pat["PSTV01"].isin(set(pk["PSTV01"]))].copy()
    return nodes, e, g, pk


# ----------------------------------------------------------------------------
# 4. GRAPH + COMMUNITY TABLES (shared by all algorithms)
# ----------------------------------------------------------------------------
def to_graph(nodes: pd.DataFrame, edges: pd.DataFrame) -> ig.Graph:
    """Graph on nodes with >= 1 edge (isolates are reported, not clustered)."""
    used = pd.Index(pd.unique(np.concatenate([edges["src"].values, edges["dst"].values])))
    g = ig.Graph(n=len(used), edges=list(zip(used.get_indexer(edges["src"]),
                                             used.get_indexer(edges["dst"]))),
                 directed=False)
    g.vs["name"] = list(used)
    g.es["weight"] = edges["weight"].astype(float).tolist()
    return g


def seeded(seed: int) -> None:
    """Make igraph's stochastic algorithms reproducible."""
    ig.set_random_number_generator(random)
    random.seed(seed)


def node_table(g: ig.Graph, memb: list[int], nodes: pd.DataFrame,
               others: dict | None = None) -> pd.DataFrame:
    df = pd.DataFrame({"PSTV01": g.vs["name"], "community": memb,
                       "degree": g.degree(), "strength": g.strength(weights="weight")})
    for name, m in (others or {}).items():   # one column per algorithm
        df[f"community_{name}"] = m
    keep = [c for c in ["PSTV01", "PSTV02", "is_tb", "died", "district_res", "PSTV05", "PSTV15"]
            if c in nodes]
    return df.merge(nodes[keep], on="PSTV01", how="left")


def community_table(nt: pd.DataFrame, pk: pd.DataFrame | None = None) -> pd.DataFrame:
    agg = nt.groupby("community").agg(
        size=("PSTV01", "size"),
        households=("PSTV02", "nunique"),
        tb_patients=("is_tb", "sum"),
        deaths=("died", "sum"),
        top_district_res=("district_res",
                          lambda s: s.mode().iat[0] if s.notna().any() else pd.NA),
    )
    agg["pct_tb"] = (100 * agg["tb_patients"] / agg["size"]).round(1)
    if pk is not None:  # R2: most shared (district, ICD) keys inside each community
        k = pk.merge(nt[["PSTV01", "community"]], on="PSTV01")
        top = (k.groupby(["community", "key"]).size().rename("n").reset_index()
                 .sort_values(["community", "n"], ascending=[True, False])
                 .groupby("community").head(3)
                 .assign(txt=lambda d: d["key"] + " (" + d["n"].astype(str) + ")")
                 .groupby("community")["txt"].agg("; ".join))
        agg["top_district|icd"] = top
    return agg.sort_values("size", ascending=False).reset_index()


def combine_edges(*edge_tables: pd.DataFrame) -> pd.DataFrame:
    """Union of relations; weights summed per pair (household 1 + clinic 0.5 = 1.5)."""
    e = pd.concat([x[["src", "dst", "weight", "relation"]] for x in edge_tables], ignore_index=True)
    a = e[["src", "dst"]].min(axis=1); b = e[["src", "dst"]].max(axis=1)
    e = e.assign(src=a, dst=b)
    out = (e.groupby(["src", "dst"], as_index=False)
             .agg(weight=("weight", "sum"),
                  relation=("relation", lambda r: "+".join(sorted(set(r))))))
    print(f"combined edges: {len(out):,}  by relation: {out['relation'].value_counts().to_dict()}")
    return out


def sens_label(cfg: Config) -> str:
    return "R2_sens_" + "_".join(re.sub(r"[^A-Z0-9]", "", c.upper()) for c in cfg.r2_sensitivity_keep_z)


def compare_partitions(ra: dict, rb: dict) -> dict:
    """ARI/NMI between two Louvain results on the patients present in both."""
    a = ra["nodes"][["PSTV01", "community"]]
    b = rb["nodes"][["PSTV01", "community"]]
    m = a.merge(b, on="PSTV01", suffixes=("_a", "_b"))
    out = {"a": ra["label"], "b": rb["label"], "common_patients": len(m),
           "only_in_a": len(a) - len(m), "only_in_b": len(b) - len(m)}
    if len(m) > 1:
        ca, cb = pd.factorize(m["community_a"])[0], pd.factorize(m["community_b"])[0]
        out["ARI"] = ig.compare_communities(list(ca), list(cb), "adjusted_rand")
        out["NMI"] = ig.compare_communities(list(ca), list(cb), "nmi")
    return out



# ----------------------------------------------------------------------------
# 4b. RUNTIME LOG (for the paper: every build / algorithm run is appended)
# ----------------------------------------------------------------------------
RUN_ID = time.strftime("%Y%m%d-%H%M%S")


def machine_info() -> dict:
    """Hardware and software the timings were measured on."""
    import platform, sys
    info = {"python": platform.python_version(), "igraph": ig.__version__,
            "pandas": pd.__version__, "numpy": np.__version__,
            "os": platform.platform(), "logical_cpus": os.cpu_count(),
            "colab": "google.colab" in sys.modules}
    try:
        with open("/proc/cpuinfo") as f:
            info["cpu"] = next(l.split(":", 1)[1].strip() for l in f if l.startswith("model name"))
    except Exception:
        info["cpu"] = platform.processor() or "unknown"
    try:
        with open("/proc/meminfo") as f:
            kb = int(next(l for l in f if l.startswith("MemTotal")).split()[1])
        info["ram_gb"] = round(kb / 1024 ** 2, 1)
    except Exception:
        info["ram_gb"] = None
    return info


def runtime_dir(cfg: Config) -> str:
    return algo_dir(cfg, "runtime")


def log_runtime(cfg: Config, rows: list[dict]) -> None:
    """Append timing rows to <out_dir>/runtime/runtime_log.csv (never overwritten)."""
    if not rows:
        return
    d = runtime_dir(cfg)
    mi = machine_info()
    df = pd.DataFrame(rows).assign(run_id=RUN_ID, timestamp=time.strftime("%Y-%m-%d %H:%M:%S"),
                                   cpu=mi["cpu"], logical_cpus=mi["logical_cpus"],
                                   ram_gb=mi["ram_gb"], python=mi["python"], igraph=mi["igraph"],
                                   colab=mi["colab"])
    p = f"{d}/runtime_log.csv"
    if os.path.exists(p):                    # union of columns: steps log different fields
        df = pd.concat([pd.read_csv(p), df], ignore_index=True)
    df.to_csv(p, index=False)
    with open(f"{d}/environment_{RUN_ID}.json", "w") as f:
        json.dump(mi, f, indent=2)


class Timer:
    """with Timer() as t: ...  ->  t.wall, t.cpu (seconds)"""
    def __enter__(self):
        self._w, self._c = time.perf_counter(), time.process_time()
        return self

    def __exit__(self, *exc):
        self.wall = time.perf_counter() - self._w
        self.cpu = time.process_time() - self._c


# ----------------------------------------------------------------------------
# 5. BUILD ALL NETWORKS ONCE  (net_01_build)
# ----------------------------------------------------------------------------
NODE_COLS = ["PSTV01", "PSTV02", "is_tb", "died", "district_res", "PSTV05", "PSTV15"]


def net_dir(cfg: Config) -> str:
    if cfg.out_dir is None:
        cfg.out_dir = default_out_dir()
    d = os.path.join(cfg.out_dir, "networks")
    os.makedirs(d, exist_ok=True)
    return d


def algo_dir(cfg: Config, algo: str) -> str:
    if cfg.out_dir is None:
        cfg.out_dir = default_out_dir()
    d = os.path.join(cfg.out_dir, algo)
    os.makedirs(d, exist_ok=True)
    return d


def network_stats(name: str, nodes: pd.DataFrame, edges: pd.DataFrame) -> dict:
    row = {"network": name, "candidate_nodes": len(nodes), "edges": len(edges)}
    if len(edges):
        g = to_graph(nodes, edges)
        comp = g.connected_components()
        row.update(nodes_with_edges=g.vcount(), isolates=len(nodes) - g.vcount(),
                   components=len(comp), largest_component=max(comp.sizes()),
                   total_weight=float(sum(g.es["weight"])), density=g.density())
        if "relation" in edges:
            for rel, n in edges["relation"].value_counts().items():
                row[f"edges_{rel}"] = int(n)
    return row


def build_networks(cfg: Config) -> pd.DataFrame:
    """Build R1, R2, R2 sensitivity and the weighted networks; write them to
    <out_dir>/networks/ with a manifest that the algorithm files read."""
    rt = []
    with Timer() as tm:
        t = load_all(cfg)
    rt.append({"step": "build", "task": "load tables (common.load_tables)", "seconds_wall": tm.wall, "seconds_cpu": tm.cpu})
    with Timer() as tm:
        pat = build_patients(t, cfg)
    rt.append({"step": "build", "task": "TB cohort", "seconds_wall": tm.wall, "seconds_cpu": tm.cpu})
    d = net_dir(cfg)

    with Timer() as tm:
        n1, e1, g1 = build_r1_household(pat, cfg)
    rt.append({"step": "build", "task": "edges R1_household", "network": "R1_household",
               "nodes": len(n1), "edges": len(e1), "seconds_wall": tm.wall, "seconds_cpu": tm.cpu})
    with Timer() as tm:
        n2, e2, g2, pk2 = build_r2_clinic_dx(t, pat, cfg)
    rt.append({"step": "build", "task": "edges R2_clinic_dx", "network": "R2_clinic_dx",
               "nodes": len(n2), "edges": len(e2), "seconds_wall": tm.wall, "seconds_cpu": tm.cpu})
    cfg_s = replace(cfg, r2_keep_z_full=cfg.r2_sensitivity_keep_z)
    sl = sens_label(cfg)
    with Timer() as tm:
        n2s, e2s, g2s, pk2s = build_r2_clinic_dx(t, pat, cfg_s)
    rt.append({"step": "build", "task": f"edges {sl}", "network": sl,
               "nodes": len(n2s), "edges": len(e2s), "seconds_wall": tm.wall, "seconds_cpu": tm.cpu})

    nets = {
        "R1_household":       (n1, e1, None),
        "R2_clinic_dx":       (n2, e2, "R2_patient_keys"),
        sl:                   (n2s, e2s, "R2_sens_patient_keys"),
        "R1R2_weighted":      (pd.concat([n1, n2]).drop_duplicates("PSTV01"),
                               combine_edges(e1, e2), "R2_patient_keys"),
        "R1R2_weighted_sens": (pd.concat([n1, n2s]).drop_duplicates("PSTV01"),
                               combine_edges(e1, e2s), "R2_sens_patient_keys"),
    }
    stats, manifest = [], {}
    for name, (nodes, edges, keys) in nets.items():
        nodes[[c for c in NODE_COLS if c in nodes]].to_csv(f"{d}/{name}_nodes.csv", index=False)
        edges.to_csv(f"{d}/{name}_edges.csv", index=False)
        manifest[name] = {"keys": keys, "main": name == "R1R2_weighted"}
        stats.append(network_stats(name, nodes, edges))
    pk2.to_csv(f"{d}/R2_patient_keys.csv", index=False)
    pk2s.to_csv(f"{d}/R2_sens_patient_keys.csv", index=False)
    g1.to_csv(f"{d}/R1_household_sizes.csv", index=False)
    g2.drop(columns="pairs").to_csv(f"{d}/R2_key_sizes.csv", index=False)
    with open(f"{d}/manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)
    with open(f"{d}/config_used.json", "w") as f:
        json.dump(asdict(cfg), f, indent=2, default=str)
    summary = pd.DataFrame(stats).set_index("network")
    summary.to_csv(f"{d}/network_summary.csv")
    log_runtime(cfg, rt)
    print(f"\nnetworks saved to {os.path.abspath(d)}/")
    return summary


def list_networks(cfg: Config) -> dict:
    p = os.path.join(net_dir(cfg), "manifest.json")
    if not os.path.exists(p):
        raise FileNotFoundError(f"{p} not found: run net_01_build.ipynb first")
    with open(p) as f:
        return json.load(f)


def load_network(cfg: Config, name: str):
    """nodes, edges, patient-key table (or None) of one saved network."""
    d = net_dir(cfg)
    dt = {"PSTV01": str, "PSTV02": str, "src": str, "dst": str, "district_res": str}
    nodes = pd.read_csv(f"{d}/{name}_nodes.csv", dtype=dt)
    edges = pd.read_csv(f"{d}/{name}_edges.csv", dtype=dt)
    keys_name = list_networks(cfg)[name]["keys"]
    keys = pd.read_csv(f"{d}/{keys_name}.csv", dtype=dt) if keys_name else None
    return nodes, edges, keys


# ----------------------------------------------------------------------------
# 6. RUN ONE ALGORITHM ON ALL NETWORKS  (used by net_louvain.py etc.)
# ----------------------------------------------------------------------------
def evaluate(g: ig.Graph, memb: list[int], resolution: float = 1.0) -> dict:
    sizes = pd.Series(memb).value_counts()
    comp = g.connected_components()
    return {
        "nodes_in_graph": g.vcount(), "edges": g.ecount(),
        "components": len(comp), "largest_component": max(comp.sizes()),
        "communities": int(sizes.size), "largest_community": int(sizes.max()),
        "median_community_size": float(sizes.median()),
        "communities_size>=3": int((sizes >= 3).sum()),
        "modularity_Q": g.modularity(memb, weights="weight", resolution=resolution),
        # 1.0 = communities are exactly the connected components (nothing split)
        "vs_components_NMI": (ig.compare_communities(memb, comp.membership, "nmi")
                              if len(comp) > 1 else float("nan")),
    }


def run_algorithm(cfg: Config, algo: str, detect, networks=None,
                  stochastic: bool = True, sweep: bool = False) -> dict:
    """detect(g, cfg, seed) -> (membership, info dict).
    Writes <out_dir>/<algo>/{network}_nodes.csv, _communities.csv,
    _resolution_sweep.csv (if sweep) and summary_metrics.csv."""
    out = algo_dir(cfg, algo)
    manifest = list_networks(cfg)
    names = list(manifest) if networks in (None, "all") else list(networks)
    results, rows = {}, []
    for name in names:
        nodes, edges, keys = load_network(cfg, name)
        print(f"\n=== {algo} | {name} ===")
        if edges.empty:
            print("  no edges -> skipped")
            continue
        g = to_graph(nodes, edges)
        with Timer() as tm:
            memb, info = detect(g, cfg, cfg.seed)
        walls, cpus = [tm.wall], [tm.cpu]          # timed runs of the algorithm itself
        m = {"network": name, "algorithm": algo, **evaluate(g, memb, cfg.resolution),
             "isolates_not_clustered": len(nodes) - g.vcount(), "seconds": round(tm.wall, 3), **info}
        t_stab = t_sweep = 0.0
        if stochastic and cfg.n_stability_runs:
            ari, nmi = [], []
            for s in range(1, cfg.n_stability_runs + 1):
                with Timer() as tm:
                    m2, _ = detect(g, cfg, s)
                walls.append(tm.wall); cpus.append(tm.cpu); t_stab += tm.wall
                ari.append(ig.compare_communities(memb, m2, "adjusted_rand"))
                nmi.append(ig.compare_communities(memb, m2, "nmi"))
            m.update(stability_ARI_mean=float(np.mean(ari)), stability_ARI_sd=float(np.std(ari)),
                     stability_NMI_mean=float(np.mean(nmi)))
        sw = None
        if sweep:
            sw = []
            with Timer() as tm:
                for r in cfg.resolution_grid:
                    m3, _ = detect(g, replace(cfg, resolution=r), cfg.seed)
                    sw.append({"resolution": r, "n_communities": len(set(m3)),
                               "modularity": g.modularity(m3, weights="weight", resolution=r)})
            t_sweep = tm.wall
            sw = pd.DataFrame(sw)
            sw.to_csv(f"{out}/{name}_resolution_sweep.csv", index=False)
        w = np.array(walls)
        m.update(seconds_median=float(np.median(w)), seconds_q1=float(np.percentile(w, 25)),
                 seconds_q3=float(np.percentile(w, 75)), timed_runs=len(w),
                 cpu_seconds_median=float(np.median(cpus)))
        log_runtime(cfg, [{"step": "detect", "task": algo, "algorithm": algo, "network": name,
                           "nodes": g.vcount(), "edges": g.ecount(),
                           "components": m["components"], "largest_component": m["largest_component"],
                           "communities": m["communities"],
                           "seconds_wall": float(np.median(w)), "seconds_wall_q1": m["seconds_q1"],
                           "seconds_wall_q3": m["seconds_q3"], "seconds_wall_min": float(w.min()),
                           "seconds_wall_max": float(w.max()), "seconds_cpu": float(np.median(cpus)),
                           "timed_runs": len(w), "seconds_stability_total": t_stab,
                           "seconds_sweep_total": t_sweep,
                           **{k: v for k, v in info.items() if k.startswith("gn_")},
                           "gn_max_edges": cfg.gn_max_edges if algo.startswith("girvan_newman") else None,
                           "gn_cut": cfg.gn_cut if algo.startswith("girvan_newman") else None}])
        nt = node_table(g, memb, nodes)
        ct = community_table(nt, keys)
        nt.to_csv(f"{out}/{name}_nodes.csv", index=False)
        ct.to_csv(f"{out}/{name}_communities.csv", index=False)
        for k in ("communities", "modularity_Q", "largest_community", "stability_ARI_mean", "seconds_median"):
            if k in m:
                v = m[k]
                print(f"  {k:22s} {v:,.4f}" if isinstance(v, float) else f"  {k:22s} {v:,}")
        rows.append(m)
        results[name] = {"metrics": m, "nodes": nt, "communities": ct, "sweep": sw, "graph": g}
    summary = pd.DataFrame(rows)
    # keep the summary of networks run earlier (e.g. GN on another network)
    p = f"{out}/summary_metrics.csv"
    if os.path.exists(p) and len(summary):
        old = pd.read_csv(p)
        summary = pd.concat([old[~old["network"].isin(summary["network"])], summary], ignore_index=True)
    summary.to_csv(p, index=False)
    with open(f"{out}/config_used.json", "w") as f:
        json.dump(asdict(cfg), f, indent=2, default=str)
    print(f"\nsaved to {os.path.abspath(out)}/")
    results["summary"] = summary
    return results


# ----------------------------------------------------------------------------
# 7. COMPARE ALGORITHMS  (net_06_compare)
# ----------------------------------------------------------------------------
def compare_algorithms(cfg: Config, algorithms=None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Reads each algorithm's results; returns (metrics table, pairwise ARI/NMI)."""
    algorithms = algorithms or cfg.algorithms
    if cfg.out_dir is None:
        cfg.out_dir = default_out_dir()
    metrics, agree = [], []
    for name in list_networks(cfg):
        memb = {}
        for a in algorithms:
            p = os.path.join(cfg.out_dir, a, f"{name}_nodes.csv")
            if os.path.exists(p):
                memb[a] = pd.read_csv(p, dtype={"PSTV01": str}).set_index("PSTV01")["community"]
        algs = list(memb)
        for i, a in enumerate(algs):
            for b in algs[i + 1:]:
                common = memb[a].index.intersection(memb[b].index)
                ca = pd.factorize(memb[a].loc[common])[0].tolist()
                cb = pd.factorize(memb[b].loc[common])[0].tolist()
                agree.append({"network": name, "a": a, "b": b, "patients": len(common),
                              "ARI": ig.compare_communities(ca, cb, "adjusted_rand"),
                              "NMI": ig.compare_communities(ca, cb, "nmi")})
    for a in algorithms:
        p = os.path.join(cfg.out_dir, a, "summary_metrics.csv")
        if os.path.exists(p):
            metrics.append(pd.read_csv(p))
        else:
            print(f"[{a}] no results yet ({p})")
    out = algo_dir(cfg, "compare")
    met = pd.concat(metrics, ignore_index=True) if metrics else pd.DataFrame()
    agr = pd.DataFrame(agree)
    met.to_csv(f"{out}/algorithm_comparison.csv", index=False)
    agr.to_csv(f"{out}/algorithm_agreement.csv", index=False)
    print(f"saved to {os.path.abspath(out)}/")
    return met, agr


# ----------------------------------------------------------------------------
# 8. DEMO DATA (synthetic, for testing only)
# ----------------------------------------------------------------------------
def make_demo_data(path: str, n_households: int = 400, seed: int = 0) -> None:
    rng = np.random.default_rng(seed)
    os.makedirs(path, exist_ok=True)
    districts = ["3171", "3172", "3273", "3578", "5171", "7371"]
    icds = ["A15", "A16", "J06", "J18", "R05", "R50", "E11", "I10", "K30", "Z09", "Z20", "Z11"]
    full_of = {"Z20": ["Z20.1", "Z20.8"], "Z11": ["Z11.1", "Z11.5"], "Z09": ["Z09.8"]}
    mem, fk = [], []
    pid = 0
    for h in range(n_households):
        d = rng.choice(districts)
        for _ in range(rng.integers(1, 6)):
            pid += 1
            p = f"P{pid:06d}"
            mem.append({"PSTV01": p, "PSTV02": f"H{h:05d}", "PSTV05": rng.choice(["1", "2"]),
                        "PSTV10": d, "PSTV15": "120.5",
                        "PSTV18": "2024" if rng.random() < .02 else ""})
            tb = rng.random() < (0.35 if h % 7 == 0 else 0.08)
            for v in range(rng.integers(1, 4)):
                code = rng.choice(["A15", "A16"]) if (tb and v == 0) else rng.choice(icds)
                vd = d if rng.random() < .85 else rng.choice(districts)
                fk.append({"PSTV01": p, "FKP02": f"V{len(fk):07d}",
                           "FKP03": f"2024-{rng.integers(1,13):02d}-{rng.integers(1,28):02d}",
                           "FKP06": vd, "FKP14A": code,
                           "FKP07": rng.choice(["PEMERINTAH", "SWASTA"]),
                           "FKP08": rng.choice(["PUSKESMAS", "KLINIK PRATAMA", "DOKTER PRAKTIK PERORANGAN"]),
                           "FKP15": rng.choice(full_of.get(code, [code + ".9"]))})
    pd.DataFrame(mem).to_csv(f"{path}/demo_kepesertaan.csv", index=False)
    pd.DataFrame(fk).to_csv(f"{path}/demo_fktp_kapitasi.csv", index=False, sep=";")
    print(f"demo data written to {path}/ ({len(mem)} members, {len(fk)} FKTP visits)")
