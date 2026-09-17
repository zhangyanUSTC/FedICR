# FedICR: reproduction guide

Code to reproduce the simulation experiments and PISA 2022 application in the FedICR paper. Run all commands from this directory, which contains run.R.

## 1. Installation

Install R with Rscript available in your terminal, then install the required packages:

~~~sh
Rscript scripts/install_packages.R
~~~

S1_SS figures also require Python 3 and Pillow:

~~~sh
python3 -m pip install Pillow
~~~

If R cannot find Pillow, set E_FICR_PYTHON to the Python executable where Pillow is installed.

## 2. Running the experiments

Test S1_M with 100 repetitions per setting:

~~~sh
Rscript run.R --task S1_M --seeds 100 --workers 10
~~~

Run all simulation experiments with 1,000 repetitions per setting:

~~~sh
Rscript run.R --task all --seeds 1000 --workers 10
~~~

For PISA, download the official student questionnaire SAS file and rebuild the inputs using the [data preparation guide](data/README.md) ([中文版](data/README_CN.md)). The raw-to-analysis workflow is:

~~~sh
Rscript scripts/pisa_data_rebuild.R
Rscript pisa/prepare_pv1.R
~~~

Skip either preparation step if its output already exists. Then run the PISA application:

~~~sh
Rscript run.R --task pisa-all --workers 10
~~~

The all task runs simulations only; pisa-all runs the three PISA stages. To run one experiment, replace the task name with a task from Section 3. Simulation parameter settings run sequentially, with repetitions parallelized within each setting.

| Option | Meaning |
| --- | --- |
| --task | Experiment name; see Section 3 |
| --seeds | Repetitions per simulation setting, using seeds 1 through this number; default 1000 |
| --workers | Number of parallel workers; default 10 |
| --clients | Restrict simulations to one client count M |
| --ss-count | Number of strengths for S1_SS or all: starts at 110, step 10; default 20, ending at 300 |
| --output | A new results directory; existing paths are rejected |
| --data | Prepared PISA data directory; default data/ |

Seed and worker counts must be positive integers. PISA does not accept --seeds or --clients. Use --help to display command-line help.

For a shorter SS test, this runs strengths 110, 120 and 130 with M = 10:

~~~sh
Rscript run.R --task S1_SS --clients 10 --ss-count 3 --seeds 20 --workers 10
~~~

Choose workers to fit available CPU and memory. On a cluster, request resources through its scheduler before running.

## 3. Tasks and paper results

Here n is the sample size per client, M is the client count, and N = nM.

| Task | Experiment settings | Paper results |
| --- | --- | --- |
| S1_M | (n,M) = (400,10), (160,25), (100,40); SS = 10 | Figure 1; supplement Table A.1 |
| S1_N | n = 200, M = 20; SS = 5,10,20 | Table 1 |
| S1_SS | (n,M) = (400,10), (100,40), (160,25); SS = 110:10:300 | Figure 2; supplement Figure A.1 |
| S2 | (n,M) = (400,10), (100,40), (160,25); SS = 10 | Table 2: group-wise MSE |
| S3 | (n,M) = (600,20), (300,40), (240,50); SS = 10; candidate K = 2:6 | Figure 3; Table 3: selected-K results |
| pisa-cva | PV1; candidate K = 2:6; 10 repeats per K | PISA group-number selection |
| pisa-compare | PV1; K = 2; 5 initializations | Tables 4–5: coefficients, group sizes and algorithm comparison |
| pisa-sensitivity | Ten matched PV sets; K = 2; 5 initializations per set | Supplement: PISA sensitivity results |

Simulations use true K = p = 4, coefficients SS times the identity matrix, three initialization rounds and an iteration limit of 30. S2 uses two starting types per round. Use 1,000 repetitions per setting for the paper design; smaller runs are for testing.

PPFL is included only in S1_N and S2.

PISA uses an iteration limit of 40 and convergence threshold 1e-4. CVA selects the smallest mean score. Comparison and sensitivity use fixed K = 2; sensitivity independently refits all ten PV sets, including PV1.

## 4. Finding and reading results

The terminal prints the results directory. A full simulation run is organized as follows:

~~~text
output/all_<timestamp>/
  run_info.rds
  S1_M/       output.txt, per-setting RData/PNG, S1M_combined.png
  S1_N/       output.txt, per-setting RData
  S1_SS/
    output.txt
    n400_m10/  per-setting RData, figure data, S1_SS_M10.png
    n100_m40/  per-setting RData, figure data, S1_SS_M40.png
    n160_m25/  per-setting RData, figure data, S1_SS_M25.png
    S1_SS_M25_M40_combined.png
  S2/         output.txt, per-setting RData/table PNG
  S3/         output.txt, per-setting RData, S3_selected_K.png, figure data
~~~

A single simulation task writes these contents directly into its run directory. PISA uses one shared directory with text/RData pairs named Pisa_CVA_*, Pisa_All_Algorithms_* and Pisa_PV_sensitivity_*.

Start with the text reports for numerical summaries and PNGs for figures. Each simulation parameter setting has one RData file containing its repetitions and summaries. run_info.rds records the run settings and environment.

| Task | Main objects in its RData file |
| --- | --- |
| S1_M | report_table: printed summary; iter_list: iteration means; runs: individual repetitions |
| S1_N, S1_SS | report_table: means/SDs; all: additional summaries; runs: individual repetitions |
| S2 | result_table: group-wise MSE summary; runs: individual repetitions |
| S3 | K_out: selected K; s_out: selection/instability summary; result_matrices: candidate s scores; selected_summary: MSE2/ARI; runs: individual repetitions |
| pisa-cva | cva_score: repeated scores; cva_summary: summary by K |
| pisa-compare | cicr: reference fit; z0: reference labels; E_FICR_values, aa_score, summary_table: algorithm/SC summaries |
| pisa-sensitivity | group_summary, coefficient_summary, difference_summary, pairwise_summary, acc_matrix |

Load a file in R, replacing the example path with an actual result file:

~~~r
e <- new.env()
load("output/YOUR_RUN/S1_M/S1_M_RES_n400_m10_gk4_p4_s1000_num3_SS10.Rdata", envir = e)
e$report_table
e$runs[[1]]
~~~

For a single-task run, omit the extra S1_M/ path component. Saved result keys E_FICR, KR and SMA denote E-FICR, CICR and SC, respectively. Use per-repetition records for further calculations rather than rounded summaries.

## 5. Interpreting the results

- MSE uses coefficient-distance label matching. Acc independently maximizes label agreement and is reported as a percentage. References are the true simulated groups/coefficients for simulations and pooled CICR for the PISA algorithm comparison; PISA Acc is not accuracy against observed true student groups.
- MSE2 is the average instance-level squared coefficient error. ARI measures chance-adjusted agreement between estimated and true simulation labels.
- In S2, rows are Oracle, CICR, E-FICR, AA, SC, SC2 and PPFL. The first four columns are group-wise MSE means; the last four are their SDs.
- In S1_M, report_table's ell column is ell/N. Raw losses and result_table's ell remain totals. Lowercase l is the residual fitting loss, not the clustering loss ell.
- The S1_M Ell curve is log(mean_ell_t - mean_ell_T), where the mean is across repetitions and T is the final iteration. Only positive gaps are plotted. The combined figure shows the first 25 MSE/Acc iterations and 20 Ell iterations; all 30 iterations remain in the saved summaries.
- The S1_SS third plot column is log(mean ell), using the mean total loss across repetitions at each SS.
- S3 selects K separately for each method by minimizing normalized grouping instability s (s_n in the paper), equal to 1 minus the Rand index. Mean_min_s and SD_min_s summarize the score at each repetition's selected K. CICR, E-FICR and AA split clients into three folds; SC splits observations within each client. Repetitions with incomplete candidate scores have no selected K. Selection percentages use valid selections as the denominator; check the reported valid counts and failed_candidates for missing results.
- PISA coefficient/group-size reports order Group 1 by the larger mean fitted value. SC best/worst MSE and Acc are selected independently and may come from different clients.
- Sensitivity aligns groups to its PV1 fit; coefficient differences are Group 2 minus Group 1. Its 45 pairwise Acc values each use optimal label matching. Across-PV coefficient summaries are descriptive, not formal multiple-imputation pooling.

S1_M produces its combined figure only when all three M values are run. S1_SS produces a figure for each available M and the M25/M40 combined figure when both have matching SS sets.

S3_selected_K.png shows the selected-K proportions, with labels indicating the percentage selecting the true K = 4.

S3 output.txt ends with a LaTeX table of MSE2 and ARI means and standard deviations for the completed M settings. Oracle and CICR appear only in the first M column. Bold values indicate the best means at the displayed precision, including ties; -- denotes a missing mean and NA an unavailable standard deviation.

## 6. Redrawing figures and interrupted runs

To regenerate figures from saved results without refitting, run this in R from the project directory:

~~~r
suppressPackageStartupMessages(source("simulation/S0_basic.R"))
options(E_FICR.root = normalizePath("."))

source("scripts/render_S1M_ijoc.R")
render_S1M_ijoc("output/YOUR_RUN/S1_M")

source("scripts/paper_figures.R")
paper_figures("output/YOUR_RUN/S1_SS")
paper_figures("output/YOUR_RUN/S3")
~~~

Use each task's results directory, not the parent all directory. S1_M requires three matching result files. Redrawing overwrites the corresponding figures/figure summaries, not raw per-setting results. S1_SS still needs Python/Pillow.

Press Ctrl+C to interrupt a run; ensure workers have exited before restarting. Completed settings stay saved, but there is no automatic resume. Rerun unfinished settings in a new output directory. If only plotting fails, redraw from the saved results. Different software environments may cause small numerical differences.
