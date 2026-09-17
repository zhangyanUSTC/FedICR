library(Matrix)
library(foreach)
library(doParallel)
library(parallel)
library(clue)
library(mclust)

library(MASS)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(combinat)
library(matrixStats)
library(fossil)

INIT_beta_M <- function(
    X, Y, K_hat, lam,
    max_samples = 100L,
    local_center_factor = 10L
) {
  M <- length(X)
  p <- ncol(X[[1]])
  X <- lapply(X, as.matrix)
  Y <- lapply(Y, as.numeric)
  
  beta_local_list <- lapply(seq_len(M), function(m) {
    X_m <- X[[m]]
    Y_m <- Y[[m]]
    
    n_choose <- min(nrow(X_m), max_samples)
    idx <- sample.int(nrow(X_m), n_choose)
    
    X_sub <- X_m[idx, , drop = FALSE]
    Y_sub <- Y_m[idx]
    
    X_diag <- Matrix::bdiag(
      lapply(seq_len(n_choose), function(i) {
        Matrix::Matrix(X_sub[i, , drop = FALSE], sparse = TRUE)
      })
    )
    
    L <- n_choose * Matrix::Diagonal(n_choose) -
      Matrix::Matrix(
        1,
        nrow = n_choose,
        ncol = n_choose,
        sparse = TRUE
      )
    
    A_TA <- Matrix::kronecker(
      L,
      Matrix::Diagonal(p)
    )
    
    betan <- as.numeric(
      solve(
        crossprod(X_diag) + lam * A_TA,
        crossprod(X_diag, Y_sub)
      )
    )
    
    beta_n <- t(
      matrix(betan, nrow = p, ncol = n_choose)
    )
    
    local_K <- min(
      local_center_factor * K_hat,
      n_choose,
      nrow(unique(as.data.frame(beta_n)))
    )
    
    if (local_K <= 1L) {
      matrix(colMeans(beta_n), nrow = 1L)
    } else {
      kmeans(
        beta_n,
        centers = local_K,
        nstart = 10
      )$centers
    }
  })
  
  beta_pool <- do.call(rbind, beta_local_list)
  
  if (nrow(beta_pool) < K_hat) {
    stop("局部中心总数小于 K_hat，无法进行全局聚类。")
  }

  kmeans(
    beta_pool,
    centers = K_hat,
    nstart = 10
  )$centers
}

INIT_beta_M_parallel <- function(X,Y,K_hat,seed,lam,
                                 max_samples = 100L,
                                 local_center_factor = 10L,
                                 n_cores = NULL) {
  set.seed(seed)
  M <- length(X)
  p <- ncol(X[[1]])
  for (m in seq_len(M)) {
    X[[m]] <- as.matrix(X[[m]])
    Y[[m]] <- as.numeric(Y[[m]])
  }
  if (is.null(n_cores)) n_cores <- getOption("E_FICR.workers", 10L)
  n_cores <- min(n_cores, M)
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  doParallel::registerDoParallel(cl)
  beta_local_list <- foreach::foreach(
    m = seq_len(M),
    .packages = c("Matrix", "stats"),
    .errorhandling = "stop"
  ) %dopar% {
    # Fixed client streams make initialization independent of worker scheduling.
    set.seed(seed + m - 1L)

    X_m <- X[[m]]
    Y_m <- Y[[m]]
    
    n_m <- nrow(X_m)
    n_choose_m <- min(n_m, max_samples)
    
    i_choose <- sample.int(
      n = n_m,
      size = n_choose_m,
      replace = FALSE
    )
    
    X_sub <- X_m[i_choose, , drop = FALSE]
    Y_sub <- Y_m[i_choose]
    
    X_blocks <- lapply(
      seq_len(n_choose_m),
      function(i) {
        Matrix::Matrix(
          X_sub[i, , drop = FALSE],
          sparse = TRUE
        )
      }
    )
    
    X_diag <- Matrix::bdiag(X_blocks)
    
    L_m <- n_choose_m * Matrix::Diagonal(n_choose_m) -
      Matrix::Matrix(
        1,
        nrow = n_choose_m,
        ncol = n_choose_m,
        sparse = TRUE
      )
    
    A_TA_m <- Matrix::kronecker(
      L_m,
      Matrix::Diagonal(p)
    )
    lhs <- crossprod(X_diag) + lam * A_TA_m
    rhs <- crossprod(X_diag, Y_sub)
    
    betan <- as.numeric(
      solve(lhs, rhs)
    )
    
    beta_n <- t(
      matrix(
        betan,
        nrow = p,
        ncol = n_choose_m
      )
    )
    
    n_distinct <- nrow(
      unique(as.data.frame(beta_n))
    )
    
    local_K <- min(
      local_center_factor * K_hat,
      n_choose_m,
      n_distinct
    )
    
    if (local_K <= 1L) {
      matrix(
        colMeans(beta_n),
        nrow = 1L
      )
    } else {
      stats::kmeans(
        beta_n,
        centers = local_K,
        nstart = 10
      )$centers
    }
  }
  
  beta_M5_list <- do.call(
    rbind,
    beta_local_list
  )
  
  if (nrow(beta_M5_list) < K_hat) {
    stop("局部中心总数小于 K_hat，无法进行全局聚类。")
  }
  beta_dis <- stats::kmeans(
    beta_M5_list,
    centers = K_hat,
    nstart = 10
  )$centers
  beta_dis
}

update_Z_list <- function(X, Y, beta_hat) {
  M <- length(X)
  K_hat <- nrow(beta_hat)
  
  z_hat <- vector(
    mode = "list",
    length = M
  )
  names(z_hat) <- names(X)
  
  for (m in seq_len(M)) {
    X_m <- X[[m]]
    Y_m <- as.numeric(Y[[m]])
    n_m <- nrow(X_m)
    
    prediction_m <- X_m %*% t(beta_hat)
    
    residual_sq <- (
      prediction_m - matrix(
        Y_m,
        nrow = n_m,
        ncol = K_hat
      )
    )^2
    
    label_m <- max.col(
      -residual_sq,
      ties.method = "first"
    )
    
    z_m <- matrix(
      0,
      nrow = n_m,
      ncol = K_hat
    )
    
    z_m[
      cbind(seq_len(n_m), label_m)
    ] <- 1
    
    z_hat[[m]] <- z_m
  }
  
  z_hat
}

l_calculate <- function(X, Y, beta_hat, z_hat) {
  sse <- sum(vapply(seq_along(X), function(m) {
    pred <- X[[m]] %*% t(beta_hat)
    labels <- max.col(z_hat[[m]], ties.method = "first")
    sum((Y[[m]] - pred[cbind(seq_along(Y[[m]]), labels)])^2)
  }, numeric(1)))
  sse / sum(lengths(Y))
}

SCORE_hat <- function(X, Y, beta0, beta_hat, z0, z_hat) {
  K <- nrow(beta0)
  K_hat <- nrow(beta_hat)
  N <- sum(lengths(Y))
  l <- sum(vapply(seq_along(X), function(m) {
    pred <- X[[m]] %*% t(beta_hat)
    label <- max.col(z_hat[[m]], ties.method = "first")
    
    sum(
      (Y[[m]] - pred[cbind(seq_along(Y[[m]]), label)])^2
    )
  }, numeric(1))) / N
  
  if (K != K_hat) {
    return(c(mse = NA, acc = NA, loss = l))
  }
  
  cost_matrix <- outer(
    seq_len(K),
    seq_len(K_hat),
    Vectorize(function(i, j) {
      sum((beta0[i, ] - beta_hat[j, ])^2)
    })
  )
  
  best_perm <- as.integer(clue::solve_LSAP(cost_matrix))
  
  mse <- sum(
    cost_matrix[cbind(seq_len(K), best_perm)]
  )
  
  # Optimize Acc label matching independently of coefficient matching.
  overlap <- Reduce(`+`, Map(function(a, b) crossprod(a, b), z0, z_hat))
  label_perm <- as.integer(clue::solve_LSAP(overlap, maximum = TRUE))
  correct <- sum(overlap[cbind(seq_len(K), label_perm)])
  
  acc <- correct / N * 100
  
  c(
    mse = mse,
    acc = acc,
    loss = l
  )
}

CICR_alg <- function(X,Y,beta_init,iters = 40,shred = 1e-4) {
  M <- length(X)
  p <- ncol(X[[1]])
  K_hat <- nrow(beta_init)
  X <- lapply(X, as.matrix)
  Y <- lapply(Y, as.numeric)
  beta_hat <- beta_init
  z_hat <- update_Z_list(X = X,Y = Y,beta_hat = beta_hat)
  
  iter_used <- 1L
  change <- NA_real_
  if (iters >= 2L) {
    for (iter in 2:iters) {
      beta_old <- beta_hat
      for (k in seq_len(K_hat)) {
        sigma_k <- matrix(0,nrow = p,ncol = p)
        XY_k <- numeric(p)
        group_size_k <- 0
        for (m in seq_len(M)) {
          X_m <- X[[m]]
          Y_m <- Y[[m]]
          w_mk <- z_hat[[m]][,k]
          group_size_k <- group_size_k + sum(w_mk)
          
          sigma_k <- sigma_k +
            crossprod(X_m, X_m * w_mk)
          XY_k <- XY_k +
            as.vector(crossprod(X_m, Y_m * w_mk))
        }
        sigma_k_reg <- sigma_k
        beta_hat[k, ] <- as.vector(
          solve(sigma_k_reg, XY_k)
        )
      }
      z_hat <- update_Z_list(X = X,Y = Y,beta_hat = beta_hat)
      iter_used <- iter
      change<-mean(abs(beta_old-beta_hat))/mean(abs(beta_hat))
      if(change<shred) break
    }
  }
  now_out <- c(change,Reduce(`+`, lapply(z_hat, colSums)),
               l_calculate(X, Y, beta_hat, z_hat)
  )
  return(
    list(
      out = now_out,
      t = iter_used,
      beta_hat = beta_hat
    )
  )
}

E_FICR_alg <- function(X,Y,beta_init,z0,beta0,iters = 40,shred = 1e-4) {
  M <- length(X)
  p <- ncol(X[[1]])
  K_hat <- nrow(beta_init)
  
  X <- lapply(X, as.matrix)
  Y <- lapply(Y, as.numeric)
  
  beta_hat <- beta_init
  z_hat <- update_Z_list(X = X,Y = Y,beta_hat = beta_hat)
  
  iter_used <- 1L
  change <- NA_real_
  if (iters >= 2L) {
    for (iter in 2:iters){
      beta_old <- beta_hat
      num_hat <- matrix(0,nrow = M,ncol = K_hat)
      grad_hat <- array(0,dim = c(M, K_hat, p))
      for (m in seq_len(M)) {
        X_m <- X[[m]]
        Y_m <- Y[[m]]
        
        for (k in seq_len(K_hat)) {
          w_mk <- z_hat[[m]][, k]
          num_hat[m, k] <- sum(w_mk)
          residual_mk <- drop(
            X_m %*% beta_hat[k, ]
          ) - Y_m
          grad_hat[m, k, ] <- as.vector(
            crossprod(
              X_m,
              w_mk * residual_mk
            )
          )
        }
      }
      
      max_positions <- apply(
        num_hat,
        2,
        which.max
      )
      for (k in seq_len(K_hat)) {
        total_num_k <- sum(num_hat[, k])
        m_select <- max_positions[k]
        X_select <- X[[m_select]]
        w_select <- z_hat[[m_select]][, k]
        
        hessian_k <- crossprod(
          X_select,
          X_select * w_select
        )
        grad_total_k <- apply(
          grad_hat[, k, , drop = FALSE],
          3,
          sum
        )
        beta_hat[k, ] <- beta_hat[k, ] -
          num_hat[m_select, k] / total_num_k *
          as.vector(
            solve(
              hessian_k,
              grad_total_k
            )
          )
      }
      z_hat <- update_Z_list(X = X,Y = Y,beta_hat = beta_hat)
      iter_used <- iter
      change<-mean(abs(beta_old-beta_hat))/mean(abs(beta_hat))
      if(change<shred) break
    }
  }
  
  score <- SCORE_hat(
    X = X,
    Y = Y,
    beta0 = beta0,
    beta_hat = beta_hat,
    z0 = z0,
    z_hat = z_hat
  )
  now_out <- c(
    iteration = iter_used,
    score
  )
  return(
    list(
      out = now_out,
      t = iter_used,
      change = change,
      beta_hat = beta_hat
    )
  )
}

SC_alg <- function(X, Y, beta_init, z0, beta0,iters, shred, single_m) {
  M <- length(X)
  K_hat <- nrow(beta_init)
  beta_hat <- beta_init
  z_hat <- update_Z_list(X = X,Y = Y,beta_hat = beta_hat)
  iter_used <- 1L
  if (iters >= 2L) {
    for (iter in seq.int(2L, iters)) {
      beta_old <- beta_hat
      X_s <- X[[single_m]]
      Y_s <- Y[[single_m]]
      for (k in seq_len(K_hat)) {
        w <- z_hat[[single_m]][, k]
        
        sigma <- crossprod(X_s, X_s * w)
        XY <- crossprod(X_s, Y_s * w)
        
        beta_hat[k, ] <- as.vector(solve(sigma, XY))
      }
      z_hat <- update_Z_list(X = X,Y = Y,beta_hat = beta_hat)
      
      change <- mean(abs(beta_hat - beta_old)) /
        max(mean(abs(beta_hat)), .Machine$double.eps)
      
      iter_used <- iter
      
      if (change < shred) break
    }
  }
  now_out <- c(
    iteration = iter_used,
    SCORE_hat(X, Y, beta0, beta_hat, z0, z_hat)
  )
  return(
    list(
      out = now_out,
      t = iter_used,
      beta_hat = beta_hat
    )
  )
}

AA_alg <- function(X, Y, beta_init, z0, beta0, iters, shred, num_init){
  M <- length(X)
  p <- ncol(X[[1]])
  K_hat <- nrow(beta_init) / num_init
  beta_M <- vector("list", M)
  client_ids<-seq_along(X)
  for (m in client_ids) {
    best_loss <- Inf
    for (ini in seq_len(num_init)) {
      rows <- ((ini - 1) * K_hat + 1):(ini * K_hat)
      fit <- SC_alg(
        X = X,
        Y = Y,
        beta_init = beta_init[rows, , drop = FALSE],
        z0 = z0,
        beta0 = beta0,
        iters = iters,
        shred = shred,
        single_m = m
      )
      if (!is.null(fit) && !is.null(fit$out)) {
        current_loss <- tail(fit$out, 1)
        
        if (is.finite(current_loss) && current_loss < best_loss) {
          best_loss <- current_loss
          beta_M[[m]] <- fit$beta_hat
        }
      }
    }
  }
  beta_M_valid <- Filter(Negate(is.null), beta_M)
  if (length(beta_M_valid) == 0) {
    stop("所有客户端的局部估计均失败。")
  }
  
  beta_pool <- do.call(rbind, beta_M_valid)
  
  if (nrow(beta_pool) < K_hat) {
    stop("可用的局部参数数量小于 K_hat。")
  }
  
  beta_hat <- kmeans(
    beta_pool,
    centers = K_hat,
    nstart = 20
  )$centers
  z_hat <- update_Z_list(X = X,Y = Y,beta_hat = beta_hat)
  
  if (!is.null(beta0) && !is.null(z0)) {
    out <- c(
      iteration = iters,
      SCORE_hat(X, Y, beta0, beta_hat, z0, z_hat)
    )
  } else {
    out <- c(
      iteration = iters,
      loss = l_calculate(X, Y, beta_hat, z_hat)
    )
  }
  list(
    out = out,
    t = iters,
    beta_hat = beta_hat
  )
}

cva_M <- function(X, Y, beta_init, K_hat,fit_alg = CICR_alg,shred = 1e-4,
                  iters = 40) {
  M <- length(X)
  num_init <- nrow(beta_init) %/% K_hat
  fold <- sample(rep(1:3, length.out = M))
  X_train1 <- X[fold == 1]
  Y_train1 <- Y[fold == 1]
  
  X_train2 <- X[fold == 2]
  Y_train2 <- Y[fold == 2]
  
  X_test <- X[fold == 3]
  Y_test <- Y[fold == 3]
  
  best_fit <- function(X_train, Y_train) {
    best <- NULL
    best_loss <- Inf
    
    for (ini in seq_len(num_init)) {
      rows <- ((ini - 1) * K_hat + 1):(ini * K_hat)
      
      fit <- fit_alg(
        X = X_train,
        Y = Y_train,
        beta_init = beta_init[rows, , drop = FALSE],
        iters = iters,
        shred = shred
      )
      
      if (!is.null(fit)) {
        loss <- fit$out[2+K_hat]
        if (is.finite(loss) && loss < best_loss) {
          best <- fit
          best_loss <- loss
        }
      }
    }
    if (is.null(best)) stop("所有初始化均估计失败。")
    best
  }
  fit1 <- best_fit(X_train1, Y_train1)
  fit2 <- best_fit(X_train2, Y_train2)
  
  update_test_z <- function(beta_hat) {
    Map(function(X_m, Y_m) {
      pred <- X_m %*% t(beta_hat)
      
      labels <- max.col(
        -abs(sweep(pred, 1, Y_m, "-")),
        ties.method = "first"
      )
      
      z_m <- matrix(0, nrow(X_m), K_hat)
      z_m[cbind(seq_len(nrow(X_m)), labels)] <- 1
      z_m
    }, X_test, Y_test)
  }
  
  z_test1 <- update_test_z(fit1$beta_hat)
  z_test2 <- update_test_z(fit2$beta_hat)
  labels1 <- unlist(
    lapply(z_test1, function(z) {
      max.col(z, ties.method = "first")
    }),
    use.names = FALSE
  )
  
  labels2 <- unlist(
    lapply(z_test2, function(z) {
      max.col(z, ties.method = "first")
    }),
    use.names = FALSE
  )
  
  tab <- table(labels1, labels2)
  
  choose2 <- function(x) {
    x * (x - 1) / 2
  }
  
  a <- sum(choose2(tab))
  
  row_pairs <- sum(choose2(rowSums(tab)))
  col_pairs <- sum(choose2(colSums(tab)))
  total_pairs <- choose2(length(labels1))
  
  b <- total_pairs - row_pairs - col_pairs + a
  
  RI <- (a + b) / total_pairs
  
  return(1-RI)
}

cva_n <- function(
    X, Y, beta_init, K_hat,
    fit_alg = CICR_alg,
    shred = 1e-4,
    iters = 40
) {
  M <- length(X)
  num_init <- nrow(beta_init) %/% K_hat
  fold_list <- lapply(X, function(X_m) {
    n_m <- nrow(X_m)
    sample(rep(1:3, length.out = n_m))
  })
  names(fold_list) <- names(X)
  make_split <- function(fold_id) {
    list(
      X = Map(
        function(X_m, fold_m) {
          X_m[fold_m == fold_id, , drop = FALSE]
        },
        X, fold_list
      ),
      Y = Map(
        function(Y_m, fold_m) {
          Y_m[fold_m == fold_id]
        },
        Y, fold_list
      )
    )
  }
  
  train1 <- make_split(1)
  train2 <- make_split(2)
  test   <- make_split(3)
  
  best_fit <- function(X_train, Y_train) {
    best <- NULL
    best_loss <- Inf
    
    for (ini in seq_len(num_init)) {
      rows <- ((ini - 1) * K_hat + 1):(ini * K_hat)
      
      fit <- fit_alg(
        X = X_train,
        Y = Y_train,
        beta_init = beta_init[rows, , drop = FALSE],
        iters = iters,
        shred = shred
      )
      
      if (!is.null(fit)) {
        loss <- fit$out[2+K_hat]
        
        if (is.finite(loss) && loss < best_loss) {
          best <- fit
          best_loss <- loss
        }
      }
    }
    
    if (is.null(best)) {
      stop("所有初始化均估计失败。")
    }
    
    best
  }
  
  fit1 <- best_fit(train1$X, train1$Y)
  fit2 <- best_fit(train2$X, train2$Y)
  
  get_test_labels <- function(beta_hat) {
    Map(function(X_m, Y_m) {
      pred <- X_m %*% t(beta_hat)
      
      max.col(
        -abs(sweep(pred, 1, Y_m, "-")),
        ties.method = "first"
      )
    }, test$X, test$Y)
  }
  
  labels1_list <- get_test_labels(fit1$beta_hat)
  labels2_list <- get_test_labels(fit2$beta_hat)
  
  labels1 <- unlist(labels1_list, use.names = FALSE)
  labels2 <- unlist(labels2_list, use.names = FALSE)
  
  tab <- table(labels1, labels2)
  choose2 <- function(x) x * (x - 1) / 2
  
  same_both <- sum(choose2(tab))
  same_1 <- sum(choose2(rowSums(tab)))
  same_2 <- sum(choose2(colSums(tab)))
  total_pairs <- choose2(length(labels1))
  
  different_both <-
    total_pairs - same_1 - same_2 + same_both
  
  RI <-1- (same_both + different_both) / total_pairs
}
