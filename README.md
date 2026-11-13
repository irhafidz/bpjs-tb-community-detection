---
title: "[Final revised paper title]"
author: "Irmasari Hafidz*, Irandi Putra Pratomo, Ahmad Rusdiansyah, Nur Aini Rakhmawati *corresponding author"
date: "[dd/11/2026]"
conference: "Submitted at the IEEE 2026 XXX"
---


## Introduction #bpjs-tb-community-detection
This repository accompanies the research paper **"[Final revised paper title]"**, accepted at the IEEE 2026 XXX , "

The study builds patient–patient networks from the BPJS Kesehatan 2019-2024 cohort tuberculosis (TB) sample data and compares four community detection algorithms on the same weighted network: **Louvain**, **Leiden**, **label propagation** and **Girvan–Newman**. Patients are linked by two relations: living in the same household, and visiting a primary-care clinic (FKTP) in the same district with the same diagnosis. The algorithms are compared on modularity, community size, stability across random seeds, agreement with each other, runtime and robustness.

## Repository Structure

```text
bpjs-tb-community-detection/
├── pipeline/   # config, data loading, network construction, one script per algorithm, notebooks 00–09
├── figures/    # high-resolution figures from the paper
├── results/    # aggregate tables (CSV/TeX); no patient-level data
├── docs/       # interpretation notes per pipeline step
└── README.md   # project introduction
```

## Dataset

- BPJS Kesehatan 2024 TB sample: 118,747 members in five tables (membership, FKTP capitation visits, FKTP non-capitation visits, FKRTL hospital visits, FKRTL secondary diagnoses).
- Study cohort: 10,839 members with an ICD-10 A15–A19 code on a 2024 encounter.
- Main network: 9,062 patients with at least one tie, 36,123 weighted edges, 464 connected components (largest: 3,078 patients).
- **Data are not included.** The sample data are available on request from Open Data BPJS Kesehatan under a data use agreement. Participant and family IDs are de-identified by BPJS. This repository holds code, aggregate results and figures only.

## Methodology

1. Load and clean the five BPJS tables; select the TB cohort (A15–A19 in 2024).
2. Build two relations from a patient × key incidence matrix *B* as *BBᵀ*:
   - **R1 household:** same family ID, weight 1.0.
   - **R2 district + diagnosis:** FKTP visit in the same province | district with the same 3-character ICD-10 code, weight 0.5 (Z codes excluded; sensitivity run keeps Z20/Z11).
   - Main network = R1 + R2 (a pair in both gets 1.5).
3. Run four algorithms in python-igraph on the same weighted network:
   - Louvain (`community_multilevel`)
   - Leiden (`community_leiden`, modularity objective)
   - Label propagation (`community_label_propagation`)
   - Girvan–Newman (`community_edge_betweenness`, run per connected component, with two dendrogram cut rules: global modularity and per-component modularity)
4. Evaluate with weighted modularity *Q*, number and size of communities, seed stability (ARI over 20 seeds), pairwise agreement (ARI, NMI), runtime (21 runs), a degree-preserving null model, a clique-overlap check and edge-weight sensitivity.

## Results

| Algorithm | Communities | Largest | *Q* | Stability (ARI) | Runtime (s) |
|---|---:|---:|---:|---:|---:|
| Louvain | 517 | 134 | 0.989 | 0.966 | 0.040 |
| Leiden | 519 | 133 | 0.989 | 0.973 | 0.116 |
| Label propagation | 857 | 54 | 0.923 | 0.913 | 0.058 |
| Girvan–Newman, global cut | 520 | 133 | 0.989 | deterministic | 206.2 |
| Girvan–Newman, per component | 1,416 | 95 | 0.889 | deterministic | 180.0 |

Louvain, Leiden and Girvan–Newman with the global cut gave nearly the same partition (ARI 0.91–0.97). Louvain and Leiden reached this in well under a second, about 5,000 times faster than Girvan–Newman. Label propagation found smaller groups and agreed least with the other methods (ARI 0.53–0.56). Observed modularity far exceeded a degree-preserving null model (0.989 vs 0.350 ± 0.003), and the Leiden partition was stable across edge weights (ARI 0.97–0.98).

## How to Run

```bash
pip install -r requirements.txt
cd pipeline
# set DATA_DIR and RES_DIR in config.py, then:
python -c "from net_common import Config, build_networks; build_networks(Config())"
python net_louvain.py; python net_leiden.py; python net_label_propagation.py
python net_girvan_newman.py
```

The notebooks `net_00` to `net_09` reproduce every table and figure in Google Colab. All parameters are fields of `Config` in `net_common.py`, and each run writes `config_used.json` next to its results.

## Acknowledgment

[Funding statement, e.g. grant number and funder.] The authors thank Open Data BPJS Kesehatan for providing the sample data.

## License and Disclaimer

Code: MIT License. Data © BPJS Kesehatan; not redistributed.
Disclaimer: Network ties indicate a shared household or a shared district and diagnosis, not confirmed contact or transmission. Results describe the BPJS 2024 TB sample, not the full population.
