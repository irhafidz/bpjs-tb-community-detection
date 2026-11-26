---
title: "[Community Detection Algorithms on a Two-Relation Patient Network from Tuberculosis Claims of Indonesia's National Health Insurance]"
author: "Irmasari Hafidz*, Irandi Putra Pratomo, Ahmad Rusdiansyah, Nur Aini Rakhmawati *corresponding author"
date: "[dd/11/2026]"
conference: "Submitted at the IEEE 2026 XXX"
---

# bpjs-tb-community-detection

This repository accompanies the paper **"[Community Detection Algorithms on a Two-Relation Patient Network from Tuberculosis Claims of Indonesia's National Health Insurance]"**, submitted to [IEEE conference name, year].

The study builds patient–patient networks from the **BPJS Kesehatan 2024 tuberculosis (TB) sample** and compares four community detection algorithms on the same weighted network: **Louvain, Leiden, label propagation and Girvan–Newman**. Patients are linked by two relations: living in the same household, and visiting a primary-care clinic (FKTP) in the same district with the same diagnosis. The algorithms are compared on modularity, community size, stability across random seeds, agreement with each other, runtime and robustness.

## Repository structure

```
bpjs-tb-community-detection/
├── pipeline/   # config, data loading, network construction, one script per algorithm, notebooks 00–09
├── figures/    # high-resolution figures from the paper
├── results/    # aggregate tables (CSV/TeX); no patient-level data
├── docs/       # interpretation notes per pipeline step
└── README.md   # project introduction
```

## Dataset

- **BPJS Kesehatan 2024 TB sample:** 118,747 members in five tables (membership, FKTP capitation visits, FKTP non-capitation visits, FKRTL hospital visits, FKRTL secondary diagnoses). Members entered the sample in 2022, 2023 or 2024 (the year of a TB diagnosis) and are followed up; the release contains their 2024 claims.
- **Study cohort:** 10,823 members with an ICD-10 A15–A19 code in any diagnosis field (FKP14A, PNK13A, FKL15A, FKL17A, FKL24A) on an encounter dated 2024 (10,839 when encounters that started in late 2023 are included).
- **Main network:** 9,062 patients with at least one tie, 36,123 weighted edges, 464 connected components (largest: 3,078 patients, 15,135 edges).
- **Household ties are rare:** of the 10,710 households in the cohort, only 303 (2.8%) have two or more members in the sample, because BPJS samples individuals; among these, 104 (34.3%) contain two or more TB-coded patients. Household ties are 420 of 36,123 edges (1.2%).
- **Data are not included.** The sample data are available on request from BPJS Kesehatan (Open Data BPJS Kesehatan) under a data use agreement. Participant and family IDs are de-identified by BPJS. This repository holds code, aggregate results and figures only.

## Methodology

1. Load and clean the five BPJS tables; select the TB cohort (A15–A19 on an encounter dated 2024).
2. Build two relations from a patient × key incidence matrix *B* as *BBᵀ*:
   - **R1 household:** same family ID (`PSTV02`), weight **1.0**.
   - **R2 district + diagnosis:** FKTP visit in the same province | district with the same 3-character ICD-10 code, weight **0.5**. The data have no facility ID, so the district is a proxy for the clinic. Z codes are excluded; a sensitivity run keeps Z20 and Z11.
   - **Main network** = R1 + R2 (a pair in both gets 1.5).
3. Run four algorithms in python-igraph on the same weighted network:
   - Louvain (`community_multilevel`)
   - Leiden (`community_leiden`, modularity objective)
   - Label propagation (`community_label_propagation`)
   - Girvan–Newman (`community_edge_betweenness`), run on each connected component (100% of patients), with two dendrogram cut rules: whole-network modularity (global cut) and per-component modularity.
4. Evaluate with weighted modularity *Q*, number and size of communities, seed stability (ARI over 20 seeds), pairwise agreement (ARI, NMI), runtime (median of 21 runs), a degree-preserving null model, a clique-overlap check and edge-weight sensitivity.

## Results

Main network (`R1R2_weighted`). Runtime: Google Colab, AMD EPYC 7B12, single thread, all algorithms on the same CPU.

| Algorithm | Communities | Largest | Q | Stability (ARI) | Runtime (s) |
|---|---:|---:|---:|---:|---:|
| Louvain | 517 | 134 | 0.989 | 0.966 | 0.025 |
| Leiden | 519 | 133 | 0.989 | 0.973 | 0.042 |
| Label propagation | 857 | 54 | 0.923 | 0.913 | 0.071 |
| Girvan–Newman, global cut | 520 | 133 | 0.989 | deterministic | 138.3 |
| Girvan–Newman, per component | 1,416 | 95 | 0.889 | deterministic | 121.0 |

- Louvain, Leiden and Girvan–Newman with the global cut gave nearly the same partition (ARI 0.91–0.97). Louvain and Leiden reached it in under 0.05 s, about 3,300–5,400 times faster than Girvan–Newman.
- For Girvan–Newman, the cut rule changed the result more than the algorithm: cutting each component at its own modularity maximum gave 1,416 communities (Q 0.889), a consequence of the resolution limit of modularity.
- Label propagation found smaller groups and agreed least with the other methods (ARI 0.53–0.56).
- Observed modularity far exceeded a degree-preserving null model (0.989 vs 0.350 ± 0.003; largest component alone 0.966 vs 0.311 ± 0.002). Only 13% of communities coincide with a single district × ICD-10 clique. The Leiden partition was stable across edge weights (ARI 0.97–0.98) and Z-code handling (ARI 0.991).
- Building the network from the five tables took about 21 s, longer than any of the fast algorithms.

## How to run

```bash
pip install -r requirements.txt
cd pipeline
# set DATA_DIR and RES_DIR in config.py, then:
python -c "from net_common import Config, build_networks; build_networks(Config())"
python net_louvain.py; python net_leiden.py; python net_label_propagation.py
python net_girvan_newman.py
```

The notebooks `net_00` to `net_09` reproduce every table and figure in Google Colab. All parameters are fields of `Config` in `net_common.py`, and each run writes `config_used.json` next to its results. For comparable runtimes, run all algorithms in one Colab session (Colab may assign a different CPU per session; the CPU is logged in `runtime/runtime_log.csv`).

## Acknowledgment

[Funding statement, e.g. grant number and funder.] The authors thank BPJS Kesehatan (Open Data BPJS Kesehatan) for providing the sample data.

## License and disclaimer

Code: MIT License. Data © BPJS Kesehatan; not redistributed.

Disclaimer: network ties indicate a shared household or a shared district and diagnosis, not confirmed contact or transmission. Results describe the BPJS 2024 TB sample, not the full population.
