# Heterogeneity and mediation of the effect of the fluid strategy in sepsis-induced hypotension

This is the code we used for the paper below.

Abgrall G, Chelly Dagdia Z, Shapiro NI, Lai C, Pham T, Zeitouni K, Monnet X. Heterogeneity and causal mediation of fluid resuscitation strategy in sepsis-induced hypotension: a secondary analysis of the CLOVERS trial.

The analyses use the data of the CLOVERS trial (Shapiro NI et al. N Engl J Med 2023;388:499-510), obtained from the NHLBI Biologic Specimen and Data Repository Information Coordinating Center (BioLINCC). The data are not in this repository. They can be requested from BioLINCC (https://biolincc.nhlbi.nih.gov/studies/petal_clovers/). We thank the CLOVERS investigators and the NHLBI PETAL Network for making the trial data available through BioLINCC, and the patients who took part in the trial. This work was prepared using CLOVERS Research Materials obtained from the NHLBI Biologic Specimen and Data Repository Information Coordinating Center and does not necessarily reflect the opinions or views of CLOVERS or the NHLBI.

The code expects the two BioLINCC files `DATASET.csv` and `DERIVED.csv` in a folder `data/` next to the code. The files written at patient level stay in `data/` and in the BART working folders, which are not tracked. The outputs kept in the notebooks are aggregates.

## Files

`CLOVERS_HTE.ipynb`: heterogeneity of the treatment effect. Table 1, Table E1, the subgroup interactions (Table E2), the Fisher randomisation test and its power (eFigure 2, Table E4), the individual effects predicted by BART with their calibration test, repeated over 100 random splits of the data, eFigure 3 and Figure 2.

`subgroup_power.py`: the three subgroup interactions with their confidence intervals and the power of the interaction tests by simulation (Methods E1). The joint model and the overlap of the subgroups, also computed here, are given in the response to the reviewers.

`CLOVERS_Mediation.ipynb`: mediation through the renal and the pulmonary pathway. The notebook builds the analysis cohorts, describes the mediators, gives the values of Table 1 on fluids and respiratory support and the effect of the strategy on each mediator, writes and runs the CMAverse script, then gives the relationship between the two mediators (Table E3) and the Bonferroni correction.

`cmaverse_models.R`: the 14 CMAverse models of the paper, refitted with the settings and seed of the notebook run. The estimates are the same, and the script also keeps the counterfactual probabilities of each bootstrap draw, from which the risk ratios are computed.

`effect_scales.R`: risk ratios, risk differences and odds ratios with their bootstrap intervals and E-values, from those probabilities (Figure 3, Table E3, Results E1 and E2, eFigure 4).

`dml_fit.R` and `dml_run.R`: double machine learning estimates of the two indirect effects with medoutcon (Results E3). `dml_fit.R` fits one pathway with one cross-fitting seed, and `dml_run.R` runs the ten seeds of each pathway and pools them.

Figure 1 and eFigure 1 are diagrams. Figure 3 and eFigure 4 were drawn from the output of `effect_scales.R`, and their drawing code is not included.

## Running the code

Everything is run from the repository folder, in this order. The first code cell of each notebook sets the path of R 4.4.3 on Windows, to be changed for another installation.

1. `CLOVERS_HTE.ipynb`. The Figure 2 cell at the top reads the BART predictions of section 10 and is run last. The randomisation test alone took 2.4 hours.
2. `subgroup_power.py` (about half an hour).
3. `CLOVERS_Mediation.ipynb`. Section D1 writes the two analysis cohorts to `data/mediation/`, where the R scripts read them. The CMAverse run of section D3 takes about 20 hours.
4. `Rscript cmaverse_models.R` (about 18 hours), then `Rscript effect_scales.R`. `cmaverse_models.R` saves a checkpoint after each model and starts again from it if stopped.
5. `Rscript dml_run.R`, which calls `dml_fit.R` twenty times (about 2.5 hours in all) and writes the pooled estimates to `output/dml_estimates.csv`.

The heterogeneity notebook writes its tables and figures to `output_hte/`, and the scripts write to `output/`.

## Notes from running the code

R is called in two ways. The CMAverse cells of the mediation notebook use rpy2, and BART is run with Rscript from a short script that the heterogeneity notebook writes to its working folder.

The notebooks were run on Windows, where R's UTF-8 output can crash rpy2's console writer, and the kernel with it. The first code cell of the mediation notebook replaces those writers.

The version of sl3 used here calls xgboost functions that xgboost 3 removed, so `dml_fit.R` has its own small xgboost learner (`Lrnr_xgb_es`).

## About the notebooks

The notebooks are shown with the outputs of the runs reported in the paper. They were run between April and September 2026 and were not rerun for this repository. Data checks, diagnostic printouts, tables and figures that are not in the paper, and fallbacks that would only have run without R were removed with their outputs. So was the negative-control mediation model, which the paper no longer reports. It was the last model fitted, so the other estimates do not depend on it. Comments, labels and printed messages were tidied, and file paths were made relative. Apart from printed messages and these removals, the code is the code that was run, and the estimates in the outputs are unchanged. The execution counts are out of order because the notebooks were run over several sessions.

In the code, `headroom` is the escalation cohort of the paper: the patients alive at 24 hours who were not on invasive ventilation at randomisation.

The CMAverse outcome model is logistic, so the ratio effects printed in the mediation notebook are odds ratios. The paper gives risk ratios, computed from the same bootstrap draws by `effect_scales.R`.

## Software

Python 3.11.9 with numpy 2.4.4 (2.4.6 for the runs after 25 July 2026), pandas 2.3.3, scipy 1.16.3, scikit-learn 1.7.2, statsmodels 0.14.6, matplotlib 3.10.9, joblib 1.5.3, tqdm 4.67.3, rpy2 3.5.17 and ipykernel 7.3.0. Sections 6 to 13 of the heterogeneity notebook were run in April 2026 in an earlier environment (Python 3.11.14, matplotlib 3.10.6) whose other versions were not recorded.

R 4.4.3 with CMAverse 0.1.0, mice 3.18.0, EValue 4.1.4, bartCause 1.0-10, dbarts 0.9-32, data.table 1.17.4 and jsonlite 2.0.0. The double machine learning scripts used a separate R library with medoutcon 0.2.5, sl3 1.4.5, xgboost 3.2.1.1, mice 3.19.0, origami 1.0.7, data.table 1.18.2.1, future 1.70.0 and delayed 0.5.0.

CMAverse, sl3 and medoutcon are installed from GitHub, at the commits used here:

```r
remotes::install_github("BS1125/CMAverse@b2ce0598ed362ddfa0db8fff5486ee43d3e73f54")
remotes::install_github("tlverse/sl3@0e8f2365bcbe54010b8120c04a7a2dcfc8119227")
remotes::install_github("nhejazi/medoutcon@307b576790e52c70d05fef8af00e9a88a2147e06")
```

The random seeds are set in the code: 42 in the notebooks and in `cmaverse_models.R`, 20261027 in `subgroup_power.py` and for the single imputation of the double machine learning, and 1001 to 1010 for the double machine learning fits.

## Contact

Questions and corrections are welcome as GitHub issues, or by email to the corresponding author of the paper.

## Licence

The code is under the MIT licence (`LICENSE`). The licence does not cover the CLOVERS data, which are available only from BioLINCC.
