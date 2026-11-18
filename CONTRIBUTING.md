# Contributing

## Adding your part

1. Clone the repo and create a branch: `git checkout -b <your-name>/<module>` (e.g. `paulina/topology`).
2. Put your code in `analysis/<module>/scripts/`. Number the scripts in the order they run (`01_build_degree.R`, `02_spearman.R`, …).
3. Replace every hardcoded path with the central config:
   ```r
   source(file.path(Sys.getenv("ADSEX_ROOT", "."), "config", "paths.R"))
   in_dir  <- file.path(DATA_DIR, "coexpr_networks")
   out_dir <- module_results("05_network_topology")
   ```
   Shell scripts use `source "${ADSEX_ROOT:-.}/config/paths.sh"` and the variables `$DATA_DIR` and `$RESULTS_DIR`.
4. Fill in your module's `README.md`: owner, inputs, outputs, how to run, and parameters (thresholds, seeds, number of permutations).
5. If you need new packages, add them to `environment.yml`.
6. Open a pull request into `main`. Someone else reviews it before merging.

## Rules

- **Never commit data or results.** Patient-level data are under an AMP-AD data use agreement. Before every commit, check `git status` to make sure no `.rds`, `.tsv` or counts files are staged.
- A module reads from `data/` or from **another module's results**, and writes **only** to `results/<its own module>/`.
- If your input is produced by another module, name that module in your README. Don't copy files between modules.
- Each script starts with a header block (title, author, date, objective, inputs, outputs), following the existing scripts.
- Fix random seeds (`set.seed`) and record them in the README.
- Write commit messages that say what changed and why.
