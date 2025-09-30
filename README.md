# Code Archive and Publication Information: DMD Critical-Period Analyses and Surgical Pathology Papers
This repository provides a stable, citable archive of analysis scripts primarily accompanying my Duchenne muscular dystrophy (DMD) studies, along with additional information and citation details for other published papers (e.g., gastric and lung pathology studies).  

All analysis scripts live in dedicated branches. The main branch contains this README and the license only.

## What this repo is (and is not)
The goal is reproducibility and scholarly transparency. Each branch includes its own README with exact instructions, versions, and data accession IDs. The main page does not reproduce those details and does not mirror datasets.

## How to use this archive
1) Choose the branch that matches the manuscript you want to reproduce.  
2) Follow that branch’s README end‑to‑end.  
3) If you need to report a technical issue, open an Issue with the exact branch name and commit hash.

## Reproducibility notes (global)
Analyses were developed in R (4.x) and/or Python where noted.  
Randomness is controlled with explicit `set.seed(...)` calls in each branch to ensure bit‑for‑bit repeatability.

## How to cite
If you use code from a branch in academic work, please cite the corresponding paper for that branch and also cite this repository (branch name and commit hash).  
Suggested wording:  
> “We used analysis scripts from the ‘<branch‑name>’ branch (commit <hash>) of the DMD Critical‑Period Analyses repository by Yamakado T.  (year), released under the Apache 2.0 license.”

## License
Copyright (c) 2025 Tetsuhiro Yamakado.  
Released under the Apache License, Version 2.0.  
You may use, reproduce, and modify the code under the terms of that license. Third‑party dependencies remain under their own licenses.

## List of Publications
1) Yamakado T, Tanei ZI, Ishikawa Y, Kimura T, Ishida Y, Tanaka S. Sporadic gastric foveolar-type adenoma with morular metaplasia. Pathol Int. 2024;74(8):489-492. doi:10.1111/pin.13462
2) Yamakado T, Sato-Yazawa H, Ishii J, et al. Spontaneous Transformation from Lung Adenocarcinoma to MYC-amplified Large Cell Neuroendocrine Carcinoma. Pathol Int. 2025;75(2):105-113. doi:10.1111/pin.13507

## Contact
Tetsuhiro Yamakado — For questions or bug reports, please open an Issue in this repository.
