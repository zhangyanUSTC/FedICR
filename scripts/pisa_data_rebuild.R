# Rebuild matched PV data independently of the estimation code.

# Construct the plausible-value field names for the requested PV indices.
pisa_pv_columns <- function(pv_indices = 1:10) {
  pv_indices <- as.integer(pv_indices)
  if (anyNA(pv_indices) || any(!pv_indices %in% 1:10)) {
    stop("pv_indices 必须是 1 到 10 之间的整数。")
  }

  unlist(
    lapply(
      c("READ", "SCIE", "MATH"),
      function(domain) sprintf("PV%d%s", pv_indices, domain)
    ),
    use.names = FALSE
  )
}

# Validate the structure and values of a reconstructed federated dataset.
validate_pisa_federated_data <- function(
    federated_data,
    pv_indices = 1:10,
    base_predictors = c("1", "ANXMAT", "MISCED", "FISCED")
) {
  if (!is.list(federated_data) || length(federated_data) == 0L) {
    stop("federated_data 必须是非空客户端列表。")
  }

  reference_names <- names(federated_data[[1]])
  valid_schema <- vapply(
    federated_data,
    function(client_data) {
      is.data.frame(client_data) &&
        identical(names(client_data), reference_names)
    },
    logical(1)
  )
  if (!all(valid_schema)) {
    stop("所有客户端必须是字段结构一致的 data.frame。")
  }

  required <- c(base_predictors, "CNT", pisa_pv_columns(pv_indices))
  missing_columns <- setdiff(required, reference_names)
  if (length(missing_columns) > 0L) {
    stop("数据缺少字段：", paste(missing_columns, collapse = ", "))
  }

  numeric_columns <- setdiff(required, "CNT")
  n_vec <- vapply(federated_data, nrow, integer(1))
  valid_values <- vapply(
    federated_data,
    function(client_data) {
      nrow(client_data) > 0L &&
        all(vapply(client_data[numeric_columns], is.numeric, logical(1))) &&
        !anyNA(client_data[required]) &&
        all(is.finite(as.matrix(client_data[numeric_columns])))
    },
    logical(1)
  )
  if (!all(valid_values)) {
    stop("客户端数据必须非空，且分析字段必须为无缺失的有限数值。")
  }

  data.frame(
    clients = length(federated_data),
    observations = sum(n_vec),
    min_client_n = min(n_vec),
    median_client_n = stats::median(n_vec),
    max_client_n = max(n_vec),
    columns = length(reference_names),
    pv_sets = length(pv_indices),
    stringsAsFactors = FALSE
  )
}

# Read the official PISA SAS file, clean it, split it by country/economy, and
# save the resulting client list as an RDS file.
rebuild_pisa_federated_data <- function(
    sas_path,
    output_path = file.path("data", "pisa_federated_data_pv.rds"),
    pv_indices = 1:10
) {
  if (!requireNamespace("haven", quietly = TRUE)) {
    stop("数据重建需要 R 包 haven，请先运行 install.packages('haven')。")
  }
  if (!file.exists(sas_path)) {
    stop("找不到 PISA SAS 数据文件：", sas_path)
  }
  if (file.exists(output_path)) {
    stop("Output already exists; choose a new output path: ", output_path)
  }

  required <- c(
    "CNT",
    "ANXMAT",
    "MISCED",
    "FISCED",
    pisa_pv_columns(pv_indices)
  )
  sas_data <- haven::read_sas(sas_path, col_select = tidyselect::all_of(required))
  missing_columns <- setdiff(required, names(sas_data))
  if (length(missing_columns) > 0L) {
    stop("SAS 数据缺少字段：", paste(missing_columns, collapse = ", "))
  }

  data_cleaned <- as.data.frame(sas_data[, required, drop = FALSE])
  data_cleaned <- data_cleaned[
    stats::complete.cases(data_cleaned),
    ,
    drop = FALSE
  ]
  data_cleaned <- data.frame(
    `1` = rep(1, nrow(data_cleaned)),
    data_cleaned,
    check.names = FALSE
  )
  federated_data <- split(data_cleaned, data_cleaned$CNT, drop = TRUE)

  validation <- validate_pisa_federated_data(
    federated_data,
    pv_indices = pv_indices
  )

  output_dir <- dirname(output_path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  saveRDS(federated_data, file = output_path, compress = FALSE)

  message("数据重建完成：", normalizePath(output_path))
  print(validation)
  federated_data
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 2L) {
    stop("Usage: Rscript scripts/pisa_data_rebuild.R [SAS_FILE] [OUTPUT_RDS]")
  }
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
  root <- dirname(dirname(normalizePath(utils::URLdecode(file_arg))))
  sas_path <- if (length(args) >= 1L) args[1L] else
    file.path(root, "data", "CY08MSP_STU_QQQ.SAS7BDAT")
  output_path <- if (length(args) >= 2L) args[2L] else
    file.path(root, "data", "pisa_federated_data_pv.rds")
  invisible(rebuild_pisa_federated_data(sas_path, output_path))
}
