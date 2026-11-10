# Code for "A pre-specified mineralization-inhibitor score separates systemic sclerosis fibroblasts within one study and does not replicate in an independent series"

Glen C. Ritschel (Ritschel Research), Aditi Sinha (Georgia Institute of Technology).
Paper: https://doi.org/10.5281/zenodo.24026978

This repository holds the code behind the numbered results the paper cites, byte-identical to the files that ran, so the sha256 values the paper gives at each point of use can be checked against them. Nothing here was rewritten for the deposit.

## Check it in two minutes

From this folder, with the packages in `requirements.txt`:

    python VERIFY_no_figure_moves_2026-10-26.py

This re-derives the held-out result, -0.2939 across 50 donors, from `PLATE_TABLE_GSE195452_optionB_2026-10-26.csv` alone, under both the current scoring module and the version before the 2026-10-26 fixes, and requires both to match. No download is needed.

    python selftest_inhibition.py
    python selftest_replication_harness.py
    python selftest_donor_group.py
    python selftest_flex_harness.py
    python selftest_stratified.py
    python selftest_cle_arms.py

Each prints PASS. Run them from this folder: they import their sibling modules from it, so `python -I` will not work.

## What is here

- Scoring and harness: `inhibition_score.py`, `replication_harness.py`, `donor_group.py`, `flex_harness.py`, `stratified_null.py`.
- `run_replication_GPL24676.py`: the one-look script for the held-out replication (RESULT_11). It refuses to run without an explicit flag, because the look could be taken only once. `prep_GPL24676_cohort.py` prepares its input.
- `run_step5_glen_score_2026-10-30.py`: the first implementation's per-donor scores for the blinded comparison in Section 7.
- `RUN_*.py` and `PREFLIGHT_*.py`: the runs behind the disease-control and detection results.
- `colab/`: cells run in Google Colab, including the GSE292979 chain from fetch (R0) through cell calling (R1r2), the donor matrix (R2) and the n = 14 sensitivity analysis. These need Colab and network access to GEO and were not re-run for this deposit.
- `gse292979_section5/`: the co-author's own implementation of the Section 5 group comparison on GSE292979, contributed by her: the n = 15 exact batch-stratified test, the n = 14 sensitivity analysis, and the n = 14 matrix needed to re-run it.
- `inhibition_score.py.bak_pre_mararo_2026-10-26`: the scoring module before three specification defects found by the second analyst were fixed. Kept because the verification above compares against it.

Code comments credit fixes to Enock Mararo, the blinded second analyst, by name, with his consent.

## Result to code

| Result | Code |
|---|---|
| RESULT_02 guards, RESULT_05, RESULT_TNAP | `inhibition_score.py` |
| RESULT_03, RESULT_08 | `donor_group.py`, `inhibition_score.py` |
| RESULT_06, RESULT_10 | `stratified_null.py` |
| RESULT_11 held-out replication | `run_replication_GPL24676.py`, `replication_harness.py`, `prep_GPL24676_cohort.py` |
| RESULT_12 | `PREFLIGHT_cd45_GPL24676_v3_FLAT_2026-10-22.py` |
| RESULT_13 | `RUN_cd45_detection_GPL24676_2026-10-22.py`, `colab/COLAB_CELL_modules_2026-10-22.txt` |
| RESULT_14 | `PREFLIGHT_cle_stratified_reachability_v3_2026-10-23.py`, `flex_harness.py` |
| RESULT_15 | `RUN_cle_detection_2026-10-23.py` |
| RESULT_16 | `RUN_cle_one_look_2026-10-23.py`, `colab/COLAB_CELL_cle_ONE_LOOK_2026-10-23.txt` |
| RESULT_18 | `colab/COLAB_CELL_leukocyte_detection_v2_2026-10-24.txt` |
| RESULT_19 | `RUN_leukocyte_detection_v3_2026-10-25.py` |
| RESULT_20 | `RUN_leukocyte_detection_v4_2026-10-25.py` |
| RESULT_21 (`dab9a7930f56a416`) | `RUN_purity_aggregation_diagnostic_v6_2026-10-25.py` |
| RESULT_22 (`0829e27b706a4ae0`) leave-one-out | `colab/COLAB_CELL_enpp1_leave_one_out_2026-10-26.txt` |
| Section 7 first implementation | `run_step5_glen_score_2026-10-30.py` |
| GSE292979 cell calling, donor matrix, n = 14 | `colab/COLAB_CELL_R0*`, `R1r2`, `R2`, `R2n14` |
| Section 5 GSE292979 group comparison, n = 15 exact test and n = 14 sensitivity | `gse292979_section5/` |

Two result numbers were each used twice in the program's working files. The paper cites results by sha256, which resolves them; the hashes above are the ones the paper means.

## What is not here, stated plainly

- **RESULT_01, RESULT_04, RESULT_07** ran in Colab on 2026-10-20 using `inhibition_score.py`; the driving cells were not saved as files.
- **RESULT_09** (sequencing batch exists) ran as a Colab pre-flight cell that was not saved.
- **RESULT_17** is a per-gene table assembled from the detection fractions of earlier results, not a separate run.
- **RESULT_23, the depth check in Section 6.6**: the code is not available. Its figures are reported in the paper but cannot be re-run from this repository.
- **The second analyst's independent implementation** is his and is not in this repository.

## License

CC BY 4.0, matching the paper and this program's other deposits.
