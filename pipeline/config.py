"""Paths, file names, column map and thresholds shared by every script.

Edit only this file. Every value can also be overridden by an environment variable, for example
    import os; os.environ["BPJS_DATA_DIR"] = "/content/drive/MyDrive/.../BPJS_TB"
before %run.
"""
import os
from pathlib import Path

# ---- folders ---------------------------------------------------------------------------------
DATA_DIR = Path(os.environ.get("BPJS_DATA_DIR", "/content/drive/MyDrive/2025_doktor_ITS/DATASET/BPJS_TB"))
OUT_DIR = Path(os.environ.get("BPJS_OUT_DIR", "/content/drive/MyDrive/2025_doktor_ITS/graph_TB_ieee"))
CACHE_DIR = OUT_DIR / "cache"          # cleaned tables (private, never commit)
RES_DIR = OUT_DIR / "results"          # aggregate result CSVs
TAB_DIR = OUT_DIR / "tables"           # LaTeX tables for the paper
FIG_DIR = OUT_DIR / "figs"             # figures for the paper
NEO_FULL_DIR = OUT_DIR / "neo4j_full"  # patient-level graph (private, never commit)
PUBLIC_DIR = OUT_DIR / "public"        # aggregate graph for GitHub and Neo4j Aura
for d in (CACHE_DIR, RES_DIR, TAB_DIR, FIG_DIR, NEO_FULL_DIR, PUBLIC_DIR):
    d.mkdir(parents=True, exist_ok=True)

# Each table is a glob pattern: every matching file in DATA_DIR is read and stacked, so all years / releases
# in the folder are used (e.g. TB2015_..., TB2024_...). Use a list of file names instead to pick files by hand.
FILES = {
    "member": os.environ.get("BPJS_FILES_MEMBER", "*kepesertaan*.csv"),
    "fktp": os.environ.get("BPJS_FILES_FKTP", "*fktpkapitasi*.csv"),
    "noncap": os.environ.get("BPJS_FILES_NONCAP", "*nonkapitasi*.csv"),
    "fkrtl": os.environ.get("BPJS_FILES_FKRTL", "*fkrtl*.csv"),
    "secdx": os.environ.get("BPJS_FILES_SECDX", "*diagnosissekunder*.csv"),
}

# ---- column map (BPJS codebook, Lampiran A-E) --------------------------------------------------
C = dict(
    # membership (Lampiran A)
    member_id="PSTV01", family_id="PSTV02", birth="PSTV03", relation="PSTV04", sex="PSTV05",
    ward_snap="PSTV07", segment_snap="PSTV08", weight="PSTV15", status_snap="PSTV17", death_year="PSTV18",
    # FKTP capitation (Lampiran B)
    fktp_id="FKP02", fktp_date="FKP03", fktp_province="FKP05", fktp_district="FKP06", fktp_kind="FKP08",
    fktp_level="FKP10", fktp_segment="FKP12", fktp_exit="FKP13", fktp_dx="FKP14A",
    ref_dest_province="FKP16", ref_dest_district="FKP17", fktp_type="FKP22",
    # FKRTL hospital (Lampiran C)
    fkrtl_id="FKL02", fkrtl_in="FKL03", fkrtl_out="FKL04", fkrtl_origin="FKP02",
    hosp_province="FKL05", hosp_district="FKL06", hosp_class="FKL09", hosp_level="FKL10",
    hosp_segment="FKL12", hosp_ward="FKL13", hosp_exit="FKL14", adm_dx="FKL15A", prim_dx="FKL17A",
    cbg="FKL19", cbg_desc="FKL19A", severity="FKL23", referrer_province="FKL25", referrer_district="FKL26",
    billed="FKL47", verified="FKL48",
    # secondary diagnoses (Lampiran E)
    sec_visit="FKL02", sec_dx="FKL24A", sec_dx_full="FKL24",
    # non-capitation (Lampiran D)
    noncap_id="PNK02", noncap_service="PNK16",
)

# ICD-10 block "Tuberculosis (A15-A19)" (WHO ICD-10 5th ed.; NHS Classifications Browser). B90 (sequelae) is
# outside the block and is reported separately, never counted as TB-coded.
TB_CODES = ("A15", "A16", "A17", "A18", "A19")
# Diagnosis fields that make a member TB-coded (codebook labels): FKP14A clinic diagnosis, FKL17A hospital primary
# diagnosis, FKL24A hospital secondary diagnosis. FKL15A (admission diagnosis, recorded at arrival) is reported as a
# sensitivity only. Set TB_INCLUDE_ADMISSION = True to count it.
TB_INCLUDE_ADMISSION = os.environ.get("BPJS_TB_INCLUDE_ADMISSION", "0") == "1"
SEED = 42
RUNS = 50                 # Leiden restarts (best of RUNS by quality)
RESOLUTION = 1.0
DX_MIN_PATIENTS = 20      # a diagnosis enters the comorbidity graph if carried by >= 20 patients
MIN_PAIR = 20             # a comorbidity tie needs >= 20 patients with both codes
MIN_FLOW = 5              # an inter-district referral tie needs >= 5 patients
MIN_CELL = 5              # public export: any count below 5 is suppressed (small-cell rule)
