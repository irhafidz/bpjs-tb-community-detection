"""Shared helpers: loading, cleaning, district harmonisation, referral validation, Leiden, saving.

Every numbered script calls load_tables(). The first call reads the five CSVs (about 30 s) and
caches the cleaned tables in OUT_DIR/cache, so later scripts in a fresh Colab start in seconds.
"""
import csv
import importlib
import pickle
import subprocess
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", category=FutureWarning)


def _ensure(pkg, pip_name=None):
    try:
        return importlib.import_module(pkg)
    except ImportError:
        subprocess.run([sys.executable, "-m", "pip", "install", "-q", pip_name or pkg], check=True)
        return importlib.import_module(pkg)


def _mount_drive():
    if Path("/content").exists() and not Path("/content/drive/MyDrive").exists():
        try:
            from google.colab import drive
            drive.mount("/content/drive")
        except Exception as e:  # running with !python: mount Drive in the notebook first
            print("Mount Google Drive in the notebook first (from google.colab import drive):", e)


_mount_drive()
_ensure("igraph", "python-igraph")
_ensure("leidenalg")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import igraph as ig  # noqa: E402
import leidenalg as la  # noqa: E402

from config import (C, CACHE_DIR, DATA_DIR, FILES, RES_DIR, RESOLUTION, RUNS, SEED,  # noqa: E402
                    TAB_DIR)

pd.set_option("display.width", 160, "display.max_columns", 40)


# ---- runtime log: every script appends one row to results/runtime.csv when it finishes ---------
def _runtime_logger():
    import atexit, os, platform, resource, time
    t0, script = time.time(), Path(sys.argv[0]).stem

    def _write():
        peak_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2   # Linux: KB -> GB
        gpu = "none (CPU only)"
        row = pd.DataFrame([dict(script=script, seconds=round(time.time() - t0, 1), peak_ram_gb=round(peak_gb, 2),
                                 cpus=os.cpu_count(), python=platform.python_version(), gpu_used=gpu,
                                 finished=time.strftime("%Y-%m-%d %H:%M"))])
        path = RES_DIR / "runtime.csv"
        row.to_csv(path, mode="a", header=not path.exists(), index=False)
    try:                                   # inside %run the kernel keeps living, so atexit never fires per script
        get_ipython()                      # noqa: F821
        print("(runtime is logged only when run as: !python <script>.py)")
    except NameError:
        atexit.register(_write)


_runtime_logger()


# ---- reading and cleaning ----------------------------------------------------------------------
def _read_one(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        head = fh.readline()
    try:
        sep = csv.Sniffer().sniff(head, delimiters=",;|\t").delimiter
    except csv.Error:
        sep = ","
    df = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False, encoding_errors="replace", low_memory=False)
    df.columns = df.columns.str.strip().str.upper().str.lstrip("﻿")
    for c in df.columns:
        df[c] = df[c].str.strip()
    return df.replace({"": np.nan}).assign(source_file=Path(path).name)


LOAD_LOG = []


def table_files(key):
    """All files of one table: FILES[key] is a glob pattern (all years) or a list of file names."""
    pat = FILES[key]
    paths = sorted(DATA_DIR.glob(pat)) if isinstance(pat, str) else [DATA_DIR / p for p in pat]
    if not paths:
        raise FileNotFoundError(f"{key}: no file matches {pat!r} in {DATA_DIR}")
    return paths


def read_table(key):
    """Read and stack every file of a table (all years). Exact duplicate rows across files are dropped and counted."""
    parts = [_read_one(p) for p in table_files(key)]
    df = pd.concat(parts, ignore_index=True, sort=False)
    data_cols = [c for c in df.columns if c != "source_file"]
    n0 = len(df)
    df = df.drop_duplicates(subset=data_cols).reset_index(drop=True)
    for p, d in zip(table_files(key), parts):
        LOAD_LOG.append(dict(table=key, file=p.name, rows=len(d), columns=d.shape[1] - 1))
    LOAD_LOG.append(dict(table=key, file="ALL (stacked)", rows=n0, columns=len(data_cols),
                         exact_duplicates_dropped=n0 - len(df)))
    print(f"  {key}: {len(parts)} file(s), {n0:,} rows, {n0 - len(df):,} exact duplicates dropped")
    if key == "member" and C["member_id"] in df:   # the same member in several sample years: keep the latest PSTV16
        yr = pd.to_numeric(df.get("PSTV16"), errors="coerce")
        dup = int(df[C["member_id"]].duplicated().sum())
        df = df.assign(_y=yr).sort_values("_y").drop_duplicates(C["member_id"], keep="last").drop(columns="_y")
        LOAD_LOG.append(dict(table=key, file="members in >1 file (kept latest PSTV16)", rows=dup))
        print(f"  member: {dup:,} repeated member rows across files (kept the latest sample year)")
    return df.reset_index(drop=True)


def to_date(x):
    d = pd.to_datetime(x, errors="coerce", format="ISO8601")
    return pd.to_datetime(x, errors="coerce", dayfirst=True) if d.isna().mean() > .5 else d


def icd3(x):
    y = x.str.upper().str.replace(r"[^A-Z0-9]", "", regex=True).str[:3]
    return y.where(y.str.fullmatch(r"[A-Z]\d\d", na=False))


def U(x):
    return x.str.upper().str.strip()


def num(x):
    return pd.to_numeric(x.str.replace(",", ".", regex=False), errors="coerce")


def has(df, key):
    return C[key] in df.columns


def prov_of(district_key):
    return district_key.split(" | ")[0]


def learn_province_map(f, h):
    """(old province, district) -> post-2022 province.

    Learned from the location columns where *_NEW is filled (FKP05/FKP05_NEW, FKL05/FKL05_NEW).
    *_NEW is not updated on every row, so per (old province, district) an updated value wins.
    The referrer and destination provinces (FKL25_NEW, FKP16_NEW) keep the pre-2022 Papua names in
    the 2024 release, so they are mapped through this table.
    """
    pairs = []
    for df, pk_, dk in [(f, "fktp_province", "fktp_district"), (h, "hosp_province", "hosp_district")]:
        po, pn, d = C[pk_], C[pk_] + "_NEW", C[dk]
        if {po, pn, d} <= set(df.columns):
            pairs.append(pd.DataFrame({"po": U(df[po]), "d": U(df[d]), "pn": U(df[pn])}).dropna())
    pmap = {}
    if pairs:
        p = pd.concat(pairs)
        p["changed"] = p.pn != p.po
        for (po_, d_), g in p.groupby(["po", "d"]):
            g2 = g[g.changed] if g.changed.any() else g
            pmap[(po_, d_)] = g2.pn.mode().iat[0]
    return pmap


def district(df, prov_key, dist_key, pmap, harmonise=True):
    """District key 'PROVINCE | DISTRICT'. Placeholders ('TIDAK DIRUJUK', '-') are not districts.

    harmonise=False uses the *_NEW province exactly as delivered (the 'before' state in the paper).
    """
    if not has(df, dist_key):
        return pd.Series(np.nan, index=df.index)
    d = U(df[C[dist_key]]).replace({"TIDAK DIRUJUK": np.nan, "-": np.nan})
    po = U(df[C[prov_key]]) if has(df, prov_key) else pd.Series("?", index=df.index)
    newc = C[prov_key] + "_NEW"
    pn = U(df[newc]) if newc in df.columns else po
    if harmonise:
        mapped = pd.Series([pmap.get((a, b)) for a, b in zip(po, d)], index=df.index)
        p = mapped.fillna(pn).fillna(po).fillna("?")
    else:
        p = pn.fillna(po).fillna("?")
    return (p + " | " + d).where(d.notna())


def _build(harmonise):
    T = {k: read_table(k) for k in FILES}
    m, f, h, s = T["member"], T["fktp"], T["fkrtl"], T["secdx"]
    f["visit_date"] = to_date(f[C["fktp_date"]])
    h["adm_date"], h["dis_date"] = to_date(h[C["fkrtl_in"]]), to_date(h[C["fkrtl_out"]])
    f["dx3"], h["dx3"], s["dx3"] = icd3(f[C["fktp_dx"]]), icd3(h[C["prim_dx"]]), icd3(s[C["sec_dx"]])
    h["adm_dx3"] = icd3(h[C["adm_dx"]]) if has(h, "adm_dx") else np.nan
    h["cost"] = num(h[C["verified"]]) if has(h, "verified") else np.nan
    h.loc[h.cost < 0, "cost"] = np.nan
    pmap = learn_province_map(f, h)
    for df, pk_, dk, out in [(f, "fktp_province", "fktp_district", "district"),
                             (f, "ref_dest_province", "ref_dest_district", "dest_district"),
                             (h, "hosp_province", "hosp_district", "district"),
                             (h, "referrer_province", "referrer_district", "referrer_district")]:
        df[out] = district(df, pk_, dk, pmap, harmonise)
    T["prov_map"] = pmap
    pd.DataFrame(LOAD_LOG).to_csv(RES_DIR / "load_log.csv", index=False)
    return T


def load_tables(harmonise=True, rebuild=False):
    """Cleaned tables as a dict: member, fktp, noncap, fkrtl, secdx (+ prov_map)."""
    path = CACHE_DIR / f"tables_{'harmonised' if harmonise else 'raw'}.pkl"
    if path.exists() and not rebuild:
        with open(path, "rb") as fh:
            return pickle.load(fh)
    print(f"reading CSVs from {DATA_DIR} (harmonise={harmonise}) ...")
    T = _build(harmonise)
    with open(path, "wb") as fh:
        pickle.dump(T, fh, protocol=pickle.HIGHEST_PROTOCOL)
    print("cached:", path)
    return T


def validate_referrals(f, h):
    """Each FKRTL.FKP02 (referral origin) checked: resolves, same member, origin on/before admission."""
    o = (f[[C["fktp_id"], C["member_id"], "visit_date", "district"]].drop_duplicates(C["fktp_id"])
         .rename(columns={C["fktp_id"]: "origin_id", C["member_id"]: "origin_member",
                          "visit_date": "origin_date", "district": "origin_district"}))
    r = (h.loc[h[C["fkrtl_origin"]].notna(),
               [C["fkrtl_id"], C["member_id"], "adm_date", C["fkrtl_origin"], "district", "referrer_district"]]
         .rename(columns={C["fkrtl_origin"]: "origin_id"}).merge(o, on="origin_id", how="left"))
    r["resolves"] = r.origin_member.notna()
    r["same_member"] = r[C["member_id"]] == r.origin_member
    r["lag_days"] = (r.adm_date - r.origin_date).dt.days
    r["valid"] = r.resolves & r.same_member & (r.lag_days >= 0)
    return r


# ---- TB-coded members ----------------------------------------------------------------------------
def tb_members(T, include_admission=None):
    """Members with >= 1 ICD-10 code A15-A19 in the TB fields (see config.TB_INCLUDE_ADMISSION)."""
    from config import TB_CODES, TB_INCLUDE_ADMISSION
    inc = TB_INCLUDE_ADMISSION if include_admission is None else include_admission
    f, h, s = T["fktp"], T["fkrtl"], T["secdx"]
    hv2p = h.set_index(C["fkrtl_id"])[C["member_id"]]
    out = (set(f.loc[f.dx3.isin(TB_CODES), C["member_id"]]) | set(h.loc[h.dx3.isin(TB_CODES), C["member_id"]])
           | set(s.loc[s.dx3.isin(TB_CODES), C["sec_visit"]].map(hv2p).dropna()))
    if inc:
        out |= set(h.loc[h.adm_dx3.isin(TB_CODES), C["member_id"]])
    return out


# ---- community detection -----------------------------------------------------------------------
def leiden_runs(g, weights="w", runs=RUNS, res=RESOLUTION):
    """Leiden (RB configuration = modularity at res 1) from seeds 0..runs-1; best by quality first."""
    parts = [la.find_partition(g, la.RBConfigurationVertexPartition, weights=weights,
                               resolution_parameter=res, seed=sd, n_iterations=-1) for sd in range(runs)]
    best = max(parts, key=lambda p: p.quality())
    return np.array(best.membership), [np.array(p.membership) for p in parts]


def louvain_seeded(g, weights="w", seed=SEED):
    import random
    random.seed(seed)
    ig.set_random_number_generator(random)
    return np.array(g.community_multilevel(weights=weights).membership)


def n_communities(memb, min_size=3):
    return int((pd.Series(memb).value_counts() >= min_size).sum())


# ---- saving ------------------------------------------------------------------------------------
def save_csv(df, name, index=False):
    df.to_csv(RES_DIR / f"{name}.csv", index=index)
    print(f"  -> results/{name}.csv")


def save_tex(df, name, caption, index=False):
    try:
        tex = df.to_latex(index=index, caption=caption, label=f"tab:{name}", escape=True,
                          float_format=lambda v: f"{v:,.3f}".rstrip("0").rstrip("."))
        (TAB_DIR / f"{name}.tex").write_text(tex)
        print(f"  -> tables/{name}.tex")
    except Exception as e:  # to_latex needs jinja2
        print("  (LaTeX skipped:", e, ")")


def show(df, n=30):
    print(df.head(n).to_string())
