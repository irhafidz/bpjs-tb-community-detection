# Pipeline

Python modules and Google Colab notebooks that build the patient networks and run the four community detection algorithms. Run the notebooks in order (00 → 09).

```text
pipeline/
├── config.py                               # paths (DATA_DIR, RES_DIR); edit for your machine
├── common.py                               # loads and cleans the five BPJS tables
├── net_common.py                           # Config, network construction, evaluation, runtime logging
├── net_louvain.py                          # Louvain: detect(g, cfg, seed) + run()
├── net_leiden.py                           # Leiden: detect(g, cfg, seed) + run()
├── net_label_propagation.py                # label propagation: detect(g, cfg, seed) + run()
├── net_girvan_newman.py                    # Girvan–Newman per component, cut = "component" | "global"
│
├── net_00_dataset.ipynb                    # dataset profile and TB cohort
├── net_01_build.ipynb                      # build networks (R1 household, R2 district + diagnosis, R1R2 weighted)
├── net_02_louvain.ipynb                    # Louvain
├── net_03_leiden.ipynb                     # Leiden
├── net_04_label_propagation.ipynb          # label propagation
├── net_05_girvan_newman.ipynb              # Girvan–Newman, both cut rules
├── net_06_compare.ipynb                    # agreement between algorithms (ARI, NMI)
├── net_07_figures.ipynb                    # figures for one algorithm (set ALGO)
├── net_07b_networks_all_algorithms.ipynb   # all algorithms on the same layout
├── net_07c_girvan_newman_rules.ipynb       # global vs per-component cut
├── net_08_runtime.ipynb                    # runtime, 21 runs per algorithm
├── net_09_robustness.ipynb                 # null model, clique check, weight sensitivity
└── README.md
```

Steps 02–05 are independent of each other. Steps 06–08 read their results. Step 09 needs only step 01.
