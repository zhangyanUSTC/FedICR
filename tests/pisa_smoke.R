# Synthetic test; no real student records are used.
root <- normalizePath(".")
source(file.path(root, "scripts", "control.R"))
data_dir <- tempfile("pisa_fixture_")
dir.create(data_dir)
set.seed(81)
d <- lapply(1:6, function(m) {
  n <- 160
  x <- data.frame(rep(1,n), CNT = paste0("Client",m), ANXMAT = rnorm(n),
                  MISCED = rnorm(n), FISCED = rnorm(n), check.names = FALSE)
  names(x)[1] <- "1"
  z <- rep(1:2, length.out = n)
  for (j in 1:10) {
    x[[paste0("PV",j,"READ")]] <- rnorm(n)
    x[[paste0("PV",j,"SCIE")]] <- rnorm(n)
    x[[paste0("PV",j,"MATH")]] <- 5*z + x$ANXMAT + .5*x$MISCED +
      x[[paste0("PV",j,"READ")]] + x[[paste0("PV",j,"SCIE")]] + rnorm(n)
  }
  x
})
names(d) <- paste0("Client",1:6)
saveRDS(d, file.path(data_dir,"pisa_federated_data_pv.rds"))
saveRDS(lapply(d, function(x) x[,c("1","CNT","ANXMAT","MISCED","FISCED","PV1READ","PV1SCIE","PV1MATH")]),
        file.path(data_dir,"pisa_federated_data_pv1.rds"))
out <- tempfile("pisa_smoke_")
run_project(root,c("--task","pisa-all","--data",data_dir,"--workers","2","--output",out))
files <- list.files(out,"[.]Rdata$",full.names=TRUE)
stopifnot(length(files)==3L)
e <- new.env()
load(files[grepl("sensitivity",files)],e)
stopifnot(nrow(e$pairwise_summary)==45L, !"ARI" %in% names(e$pairwise_summary),
          all(diag(e$acc_matrix)==100), all(is.finite(e$beta_array)),
          all(e$group_summary$Group1_n + e$group_summary$Group2_n == 960L))
cat("PASS: three PISA pipelines, raw result files and 45 pairwise Acc values.\n")
