argv <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", grep("^--file=", argv, value = TRUE)[1L])
root <- dirname(normalizePath(utils::URLdecode(gsub("~+~", " ", script, fixed = TRUE))))
source(file.path(root, "scripts", "control.R"))
tryCatch(run_project(root, commandArgs(trailingOnly = TRUE)), error = function(e) {
  message("Replication stopped: ", conditionMessage(e))
  quit(status = 1L)
})
