# Rebuilding the PISA 2022 data

[中文版](README_CN.md)

The repository provides code to reconstruct the analysis inputs from the official OECD student questionnaire file. Source and derived student-level data are not distributed with the repository. Simulation experiments do not require these files.

## 1. Download and place the source data

Open the [OECD PISA data-file download page](https://webfs.oecd.org/pisa2022/index.html). Under **PISA 2022 Data → SAS Data Files (Compressed)**, download **Student questionnaire data file**. Extract the archive and place the SAS file here:

~~~text
data/CY08MSP_STU_QQQ.SAS7BDAT
~~~

Use the 2022 student questionnaire file, not another year's data or a school, teacher or cognitive-item file. The source is several gigabytes; allow sufficient memory and disk space. The [OECD PISA 2022 Database](https://www.oecd.org/en/data/datasets/pisa-2022-database.html) also links to the downloads.

## 2. Install and rebuild

Run from the project root, where run.R is located:

~~~sh
Rscript scripts/install_packages.R
Rscript scripts/pisa_data_rebuild.R
Rscript pisa/prepare_pv1.R
~~~

The first preparation command creates the ten-PV dataset; the second extracts PV1 without changing the sample.

| Generated file | Contents | Analysis task |
| --- | --- | --- |
| data/pisa_federated_data_pv.rds | Ten matched PV sets; 35 columns per client | pisa-sensitivity |
| data/pisa_federated_data_pv1.rds | PV1 only; 8 columns per client | pisa-cva, pisa-compare |

Both preparation scripts refuse to overwrite existing outputs. Skip a step if its prepared file already exists.

## 3. What the reconstruction does

The script reads CNT, ANXMAT, MISCED, FISCED and all 30 READ/SCIE/MATH PV columns. It jointly removes rows missing any of these fields, adds an intercept column named 1, and splits rows by CNT into a list of country/economy data frames. Within-client student order is preserved. It does not average PVs or standardize the retained variables; final student weights and replicate weights are not used.

CNT identifies the client. For each matched PV set j, mathematics is the response and reading/science are predictors alongside ANXMAT, MISCED, FISCED and the intercept. Keep the generated column order and common student sample; do not filter or sort PV sets independently.

## 4. Check the sample and run PISA

Reconstruction prints a validation summary. The study sample has **75 clients, 451,851 students, 35 columns and 10 PV sets**. PV1 keeps the same clients/students and has 8 columns.

If counts differ, check the source version and missing-value handling; do not alter records merely to match the expected size.

~~~sh
Rscript run.R --task pisa-all --workers 10
~~~

If prepared files are stored elsewhere, add --data followed by that directory. Quote paths containing spaces. See the [main reproduction guide](../README.md) for individual tasks and result interpretation.
