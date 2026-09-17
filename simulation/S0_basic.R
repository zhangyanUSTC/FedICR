library(MASS)
library(clue)
library(doParallel)
library(ggplot2)
library(gridExtra)
library(foreach)
library(cowplot)
library(combinat)
library(matrixStats)
library(mclust)
library(Matrix)
library(fossil)

# Support both legacy and position-specific ggplot legend boxes.
simulation_legend <- function(plot, position = "bottom") {
  g <- if (inherits(plot, "gtable")) plot else ggplot2::ggplotGrob(plot)
  for (name in c(paste0("guide-box-", position), "guide-box")) {
    indices <- which(g$layout$name == name)
    for (i in indices) {
      if (!inherits(g$grobs[[i]], "zeroGrob")) return(g$grobs[[i]])
    }
  }
  NULL
}

# Normalize legacy result keys in memory; never rewrite saved files.
E_FICR_result_names <- function(x) {
  if (!is.null(names(x))) {
    names(x) <- sub("^(FedICR|FICR|DKR|E-FICR)$", "E_FICR", names(x))
  }
  x
}

print_simulation_table <- function(x, ..., pad = TRUE) {
  display_method <- function(v) sub("^E_FICR$", "E-FICR", v)
  format_number <- function(v) {
    if (!is.numeric(v)) return(if (is.character(v)) display_method(v) else v)
    integer_value <- is.finite(v) & v == trunc(v)
    v <- round(v, 6L)
    v[!is.na(v) & v == 0] <- 0
    out <- formatC(v, format = "f", digits = 6L, drop0trailing = !pad)
    out[integer_value] <- formatC(v[integer_value], format = "f", digits = 0L)
    out
  }
  if (is.data.frame(x)) {
    x[] <- lapply(x, format_number)
  } else {
    x[] <- format_number(x)
  }
  if (!is.null(rownames(x))) rownames(x) <- display_method(rownames(x))
  if (!is.null(names(x))) names(x) <- display_method(names(x))
  print(x, quote = FALSE, ...)
  invisible(NULL)
}


mean_index<-function(res,iters){
  out<-NULL
  for (j in 1:iters) {
    sum<-colMeans(res[res[, 4] == j, , drop = FALSE])
    out<-rbind(out,c(j,sum))
  }
  return(out)
}


INIT_beta_1K<-function(X,Y,K_hat,seed,lam,select_m){
  M<-dim(X)[1];p<-dim(X)[3];n<-dim(X)[2]
  if(n>100){
    n_choose<-100
  }else {n_choose<-n}
  
  combinations <- combn(n_choose, 2)
  C <- matrix(0, nrow = ncol(combinations), ncol = n_choose)
  for (k in 1:ncol(combinations)) {
    i <- combinations[1, k]
    j <- combinations[2, k]
    C[k, i] <- 1
    C[k, j] <- -1
  }
  A_TA<-kronecker(t(C)%*%C,diag(p))
  set.seed(seed)
  i_choose<-sample(1:n, n_choose, replace = FALSE)
  X_diag<-t(X[select_m,i_choose[1],])
  for (i in i_choose[2:n_choose]) {
    X_diag = bdiag(X_diag, t(X[select_m,i,]))
  }
  X_diag <- as.matrix(X_diag)
  betan<-solve(t(X_diag)%*%X_diag+lam*A_TA)%*%t(X_diag)%*%Y[select_m,i_choose]
  beta_n<-t(matrix(betan,ncol = n_choose))
  beta_n<-kmeans(beta_n,centers = K_hat,nstart=10)$centers
  return(beta_n)
}

INIT_beta_M<-function(X,Y,K_hat,seed,lam){
  set.seed(seed)
  M<-dim(X)[1];p<-dim(X)[3];n<-dim(X)[2]
  if(n>100){
    n_choose<-100
  }else {n_choose<-n}
  combinations <- combn(n_choose, 2)
  C <- matrix(0, nrow = ncol(combinations), ncol = n_choose)
  for (k in 1:ncol(combinations)) {
    i <- combinations[1, k]
    j <- combinations[2, k]
    C[k, i] <- 1
    C[k, j] <- -1
  }
  A_TA<-kronecker(t(C)%*%C,diag(p))
  
  beta_M5<-NULL
  
  for (m in 1:M) {
    i_choose<-sample(1:n, n_choose, replace = FALSE)
    X_diag<-t(X[m,i_choose[1],])
    for (i in i_choose[2:n_choose]) {
      X_diag = bdiag(X_diag, t(X[m,i,]))
    }
    X_diag <- as.matrix(X_diag)
    betan<-solve(t(X_diag)%*%X_diag+lam*A_TA)%*%t(X_diag)%*%Y[m,i_choose]
    beta_n<-t(matrix(betan,ncol = n_choose))
    beta_M5<-rbind(beta_M5,kmeans(beta_n,centers =5*K_hat,nstart=10)$centers)
  }
  beta_dis<-kmeans(beta_M5,centers = K_hat,nstart=10)$centers
  return(beta_dis)
}

INIT_beta_pool <- function(X, Y, K_hat, seed, lam = .001,
                           pool_size = 200L) {
  M <- dim(X)[1]; n <- dim(X)[2]; p <- dim(X)[3]
  stopifnot(all(dim(Y) == c(M, n)), K_hat >= 2L, K_hat <= M * n,
            pool_size >= K_hat, pool_size == as.integer(pool_size), lam > 0,
            all(is.finite(X)), all(is.finite(Y)))
  X_all <- do.call(rbind, lapply(seq_len(M), function(m)
    matrix(X[m, , ], nrow = n, ncol = p)))
  Y_all <- as.vector(t(Y))
  sample_size <- min(as.integer(pool_size), M * n)
  set.seed(seed)
  selected <- sample.int(M * n, sample_size, replace = FALSE)
  X_sample <- X_all[selected, , drop = FALSE]
  Y_sample <- Y_all[selected]
  
  X_diag <- Matrix::sparseMatrix(
    i = rep(seq_len(sample_size), each = p),
    j = seq_len(sample_size * p), x = as.vector(t(X_sample)),
    dims = c(sample_size, sample_size * p))
  graph_laplacian <- sample_size * Matrix::Diagonal(sample_size) -
    Matrix::Matrix(1, sample_size, sample_size, sparse = TRUE)
  penalty <- kronecker(graph_laplacian, Matrix::Diagonal(p))
  beta_vector <- Matrix::solve(crossprod(X_diag) + lam * penalty,
                               crossprod(X_diag, Y_sample))
  beta_sample <- matrix(as.numeric(beta_vector), nrow = sample_size,
                        ncol = p, byrow = TRUE)
  kmeans(beta_sample, centers = K_hat, nstart = 10)$centers
}

INIT_random<-function(X,Y,beta0,K_hat,seed){
  set.seed(seed)
  M<-dim(X)[1];p<-dim(X)[3]
  beta_init<-NULL
  for(k in 1:(K_hat)){
    beta_init<-rbind(beta_init,c(runif(p,min=-beta0[1,1],max=beta0[1,1])))
  }
  return(beta_init)
}

SCORE_hat<-function(X,Y,beta0,beta_hat,z0,z_hat){
  M<-dim(X)[1];n<-dim(X)[2];p<-dim(X)[3]
  K<-dim(beta0)[1];K_hat<-dim(beta_hat)[1]
  l<-0
  for (m in 1:M) {
    l<-l+sum((Y[m,]-diag(z_hat[m,,]%*%beta_hat%*%t(X[m,,])))^2)
  }
  l<-l/M/n
  if(K!=K_hat){acc<-NA;mse<-NA;ell<-NA;mse_group<-rep(NA_real_, K);
  acc_group<-rep(NA_real_, K)} else {
    cost_matrix <- matrix(NA, nrow = K, ncol = K)
    for (i in 1:K) {
      for (j in 1:K) {
        cost_matrix[i, j] <- sum((beta0[i, ] - beta_hat[j, ])^2)
      }
    }
    best_perm <- solve_LSAP(cost_matrix)
    mse <- sum(diag(cost_matrix[cbind(1:K, best_perm)]))
    mse_group<-cost_matrix[cbind(1:K, best_perm)]
    # Optimize Acc label matching independently of coefficient matching.
    agreement <- matrix(0, nrow = K, ncol = K)
    for (m in seq_len(M)) {
      agreement <- agreement + crossprod(z0[m,,], z_hat[m,,])
    }
    label_perm <- solve_LSAP(agreement, maximum = TRUE)
    acc<-sum(z0*z_hat[,,label_perm, drop = FALSE])/M/n*100
    acc_group<-apply(z0*z_hat[,,label_perm, drop = FALSE], c(3), sum)
    
    ell<-0
    for (m in 1:M) {
      ell<-ell+sum(diag(z_hat[m,,best_perm]%*%beta0%*%t(X[m,,])-z0[m,,]%*%beta0%*%t(X[m,,]))^2)
    }
  }
  mse2 <- 0
  clusters0 <- integer(M * n)
  clusters1 <- integer(M * n)
  for (m in seq_len(M)) {
    mse2 <- mse2 + sum((z0[m,,] %*% beta0 - z_hat[m,,] %*% beta_hat)^2)
    ids <- (m - 1L) * n + seq_len(n)
    clusters0[ids] <- max.col(z0[m,,], ties.method = "first")
    clusters1[ids] <- max.col(z_hat[m,,], ties.method = "first")
  }
  mse2 <- mse2 / M / n
  ARI <- mclust::adjustedRandIndex(clusters1, clusters0)
  return(c(mse,acc,mse_group,acc_group,mse2,ARI,l,ell))
}

l_calculate<-function(X,Y,beta_hat,z_hat){
  M<-dim(X)[1];n<-dim(X)[2]
  l<-0
  for (m in 1:M) {
    l<-l+sum((Y[m,]-diag(z_hat[m,,]%*%beta_hat%*%t(X[m,,])))^2)
  }
  l<-l/M/n
  return(l)
}

RI_calculate<-function(z1,z2){
  M<-dim(z1)[1]
  clusters0 <- NULL;clusters1 <- NULL
  for (m in 1:M) {
    clusters0<-c(clusters0,apply(z1[m,,], 1, which.max))
    clusters1<-c(clusters1,apply(z2[m,,], 1, which.max))
  }
  ARI<-rand.index(clusters1, clusters0)
  return(ARI)
}

update_Z <- function(X, Y, beta_hat) {
  M <- dim(X)[1]; n <- dim(X)[2]; K_hat <- nrow(beta_hat)
  z_hat <- array(0, dim = c(M, n, K_hat))
  for (m in 1:M) {
    residuals <- abs(Y[m, ] - X[m, , ] %*% t(beta_hat)) 
    min_indices <- apply(residuals, 1, which.min) 
    z_hat[cbind(rep(m, n), 1:n, min_indices)] <- 1
  }
  return(z_hat)
}

get_best_res <- function(old_res, new_res, metric_col = 17) {
  if (is.null(new_res) || ncol(new_res) <= 3) return(old_res)
  if (is.null(old_res) || ncol(old_res) <= 3) return(new_res)
  if (old_res[nrow(old_res), metric_col] > new_res[nrow(new_res), metric_col]) {
    return(new_res)
  }
  return(old_res)
}

get_best_res_2 <- function(old_res, new_res, metric_col = 17) {
  if (is.null(new_res) || ncol(new_res) <= 3) return(old_res)
  if (is.null(old_res) || ncol(old_res) <= 3) return(new_res)
  if (old_res[metric_col] > new_res[metric_col]) {
    return(new_res)
  }
  return(old_res)
}


generate_X <- function(M, n, beta, Sigma, batch_size = 10000) {
  p <- ncol(beta)
  N <- M * n
  
  pairs <- combn(nrow(beta), 2)
  beta_diff <- beta[pairs[1, ], , drop = FALSE] -
    beta[pairs[2, ], , drop = FALSE]
  diff_norm <- sqrt(rowSums(beta_diff^2))
  
  X_mat <- matrix(NA_real_, N, p)
  n_saved <- 0
  
  while (n_saved < N) {
    
    x <- MASS::mvrnorm(
      n = max(batch_size, N - n_saved),
      mu = rep(0, p),
      Sigma = Sigma
    )
    x <- matrix(x, ncol = p)
    
    x_norm <- sqrt(rowSums(x^2))
    
    margin <- abs(x %*% t(beta_diff)) /
      outer(x_norm, diff_norm)
    
    keep <- which(
      x_norm >= 0.1 &
        apply(margin, 1, min) >= 0.01
    )
    
    if (length(keep) > 0) {
      take <- min(length(keep), N - n_saved)
      
      X_mat[(n_saved + 1):(n_saved + take), ] <-
        x[keep[1:take], , drop = FALSE]
      
      n_saved <- n_saved + take
    }
  }
  
  X <- array(NA_real_, dim = c(M, n, p))
  
  for (j in 1:p) {
    X[, , j] <- matrix(
      X_mat[, j],
      nrow = M,
      ncol = n,
      byrow = TRUE
    )
  }
  
  return(X)
}


E_FICR_alg <- function(X,Y,beta_init,z0,beta0,iters,shred,save_history = FALSE,save_estimate=FALSE) {
  tryCatch({
    M<-dim(X)[1];n<-dim(X)[2];p<-dim(X)[3];K_hat<-dim(beta_init)[1]
    beta_hat <- beta_init
    z_hat <- update_Z(X, Y, beta_hat)
    if (save_history) {
      out <- c(1, SCORE_hat(X, Y, beta0, beta_hat, z0, z_hat))
      if (iters == 1L) out <- matrix(out, nrow = 1L)
    }
    t <- 1
    if (iters >= 2) {
      for (t in 2:iters) {
        num_hat <- array(0, dim = c(M, K_hat))
        grad_hat <- array(0, dim = c(M, K_hat, p))
        for (m in 1:M) {
          for (k in 1:K_hat) {
            w <- z_hat[m,,k]
            num_hat[m,k] <- sum(w)
            grad_hat[m,k,] <- crossprod(X[m,,], X[m,,] * w) %*% beta_hat[k,] -
              crossprod(X[m,,], Y[m,] * w)
          }
        }
        max_positions <- apply(num_hat, 2, which.max)
        for (k in 1:K_hat) {
          m_select<-max_positions[k]
          w <- z_hat[m_select,,k]
          beta_hat[k,] <- beta_hat[k,] - num_hat[m_select,k] / sum(num_hat[,k]) *
            solve(crossprod(X[m_select,,], X[m_select,,] * w),
                  colSums(grad_hat[,k,]))
        }
        z_hat <- update_Z(X, Y, beta_hat)
        if (save_history) {
          out <- rbind(out,c(t, SCORE_hat(X, Y, beta0, beta_hat, z0, z_hat)))
        }
        if(mean(abs(beta_init-beta_hat))/mean(abs(beta_hat))<shred) break
        beta_init<-beta_hat
      }
    }
    if (save_history) {
      if (t < iters) {
        last_score <- out[nrow(out), -1]
        for (tt in (t + 1):iters) {
          out <- rbind(out, c(tt, last_score))
        }
      }
    } else {out <- c(t,SCORE_hat(X, Y, beta0, beta_hat, z0, z_hat))}
    
    if(save_estimate){
      return(list(
        beta_hat=beta_hat,
        z_hat=z_hat,
        l=l_calculate(X,Y,beta_hat,z_hat),
        out=out,
        t=t
      ))
    }
    return(list(
      out = out,
      t = t
    ))
  }, error = function(e) {
    return(NULL)
  })
}

CICR_alg <- function(X,Y,beta_init,z0,beta0,iters,shred,save_history = FALSE,save_estimate=FALSE) {
  tryCatch({
    M<-dim(X)[1];n<-dim(X)[2];p<-dim(X)[3];K_hat<-dim(beta_init)[1]
    beta_hat <- beta_init
    z_hat <- update_Z(X, Y, beta_hat)
    if (save_history) {
      out <- c(1, SCORE_hat(X, Y, beta0, beta_hat, z0, z_hat))
      if (iters == 1L) out <- matrix(out, nrow = 1L)
    }
    t <- 1
    if (iters >= 2) {
      for (t in 2:iters) {
        for(k in 1:K_hat){
          sigma<-array(0,dim = c(p, p))
          XY<-0
          for (m in 1:M) {
            w <- z_hat[m,,k]
            sigma <- sigma + crossprod(X[m,,], X[m,,] * w)
            XY <- XY + crossprod(X[m,,], Y[m,] * w)
          }
          beta_hat[k,] <- solve(sigma, XY)
        }
        z_hat <- update_Z(X, Y, beta_hat)
        if (save_history) {
          out <- rbind(out,c(t, SCORE_hat(X, Y, beta0, beta_hat, z0, z_hat)))
        }
        if(mean(abs(beta_init-beta_hat))/mean(abs(beta_hat))<shred) break
        beta_init<-beta_hat
      }
    }
    if (save_history) {
      if (t < iters) {
        last_score <- out[nrow(out), -1]
        for (tt in (t + 1):iters) {
          out <- rbind(out, c(tt, last_score))
        }
      }
    } else {out <- c(t,SCORE_hat(X, Y, beta0, beta_hat, z0, z_hat))}
    
    if(save_estimate){
      return(list(
        beta_hat=beta_hat,
        z_hat=z_hat,
        l=l_calculate(X,Y,beta_hat,z_hat)
      ))
    }
    return(list(
      out = out,
      t = t
    ))
  }, error = function(e) {
    return(NULL)
  })
}

SC_alg<-function(X,Y,beta_init,z0,beta0,iters,shred,single_m,save_history = FALSE,save_estimate=FALSE) {
  tryCatch({
    M<-dim(X)[1];n<-dim(X)[2];p<-dim(X)[3];K_hat<-dim(beta_init)[1]
    beta_hat <- beta_init
    residuals <- abs(Y[single_m, ] - X[single_m, , ] %*% t(beta_hat))
    z_single <- matrix(0, nrow = n, ncol = K_hat)
    z_single[cbind(seq_len(n), apply(residuals, 1, which.min))] <- 1
    if (save_history) {
      z_hat <- update_Z(X, Y, beta_hat)
      out <- c(1, SCORE_hat(X, Y, beta0, beta_hat, z0, z_hat))
      if (iters == 1L) out <- matrix(out, nrow = 1L)
      beta_iter<-array(0 ,dim = c(iters, K_hat, p));beta_iter[1,,]<-beta_init
    }
    t <- 1
    if (iters >= 2) {
      for (t in 2:iters) {
        for(k in 1:K_hat){
          w <- z_single[,k]
          sigma <- crossprod(X[single_m,,], X[single_m,,] * w)
          XY <- crossprod(X[single_m,,], Y[single_m,] * w)
          beta_hat[k,] <- solve(sigma, XY)
        }
        residuals <- abs(Y[single_m, ] - X[single_m, , ] %*% t(beta_hat))
        z_single[,] <- 0
        z_single[cbind(seq_len(n), apply(residuals, 1, which.min))] <- 1
        if (save_history) {
          z_hat <- update_Z(X, Y, beta_hat)
          out <- rbind(out,c(t, SCORE_hat(X, Y, beta0, beta_hat, z0, z_hat)))
          beta_iter[t,,]<-beta_hat
        }
        if(mean(abs(beta_init-beta_hat))/mean(abs(beta_hat))<shred) break
        beta_init<-beta_hat
      }
    }
    if (save_history) {
      if (t < iters) {
        last_score <- out[nrow(out), -1]
        for (tt in (t + 1):iters) {
          out <- rbind(out, c(tt, last_score))
          beta_iter[tt,,]<-beta_hat
        }
      }
    }
    z_hat <- update_Z(X, Y, beta_hat)
    if (!save_history) {out <- c(t,SCORE_hat(X, Y, beta0, beta_hat, z0, z_hat))}
    
    if(save_history){
      history <- list(
        out = out,
        t = t,
        beta_iter=beta_iter
      )
      if (save_estimate) {
        history$beta_hat <- beta_hat
        history$z_hat <- z_hat
        history$l <- l_calculate(X, Y, beta_hat, z_hat)
      }
      return(history)
    }
    
    if(save_estimate){
      return(list(
        beta_hat=beta_hat,
        z_hat=z_hat,
        l=l_calculate(X,Y,beta_hat,z_hat),
        out=out,
        t=t
      ))
    }
    return(list(
      out = out,
      t = t,
      beta_hat=beta_hat
    ))
  }, error = function(e) {
    return(NULL)
  })
}

AA_alg<-function(X,Y,beta_init,z0,beta0,iters,shred,num_init,seed,save_history = FALSE,save_estimate=FALSE) {
  tryCatch({
    M<-dim(X)[1];n<-dim(X)[2];p<-dim(X)[3];K_hat<-dim(beta_init)[1]/num_init
    z_hat<-array(0 ,dim = c(M, n, K_hat))
    if (save_history) {
      beta_M<-array(0 ,dim = c(M,iters, K_hat, p))
    } else {beta_M<-array(0 ,dim = c(M,K_hat, p))}
    for (m in 1:M ) {
      SC_out<-NULL
      SC_loss<-Inf
      for (ini in 1:num_init) {
        fit <- SC_alg(
          X, Y,
          beta_init = beta_init[((ini - 1) * K_hat + 1):(ini * K_hat), ],
          z0 = z0, beta0 = beta0, iters = iters, shred = shred,
          single_m = m, save_history = save_history,
          save_estimate = save_estimate
        )
        if (is.null(fit)) next
        if (save_history) {
          beta_m <- fit$beta_iter
          candidate_out <- fit$out[iters, ]
        } else {
          beta_m <- fit$beta_hat
          candidate_out <- fit$out
        }
        candidate_loss <- candidate_out[length(candidate_out) - 1L]
        if (is.finite(candidate_loss) && candidate_loss < SC_loss) {
          SC_out <- candidate_out
          SC_loss <- candidate_loss
          if (save_history) beta_M[m, , , ] <- beta_m
          else beta_M[m, , ] <- beta_m
        }
      }
    }
    set.seed(seed)
    if (save_history) {
      beta_hat<-array(0 ,dim = c(iters, K_hat, p))
      out<-NULL
      for(t in 1:iters){
        beta_M_t<-beta_M[,t,,]
        beta_M_t<-beta_M_t[apply(beta_M_t,1, function(sub_matrix) !all(sub_matrix == 0)),,]
        beta_M_t<-matrix(beta_M_t,ncol = p)
        beta_hat[t,,]<-kmeans(beta_M_t,centers = K_hat,nstart = 20)$centers
        z_hat <- update_Z(X, Y, beta_hat[t,,])
        out<-rbind(out,c(t,SCORE_hat(X,Y,beta0,beta_hat[t,,],z0,z_hat)))
      }
    } else {
      valid_clients <- apply(beta_M, 1, function(sub_matrix) !all(sub_matrix == 0))
      beta_M_t <- matrix(beta_M[valid_clients,,,drop = FALSE], ncol = p)
      if (nrow(beta_M_t) < K_hat) stop("Too few valid local AA estimates")
      beta_hat<-kmeans(beta_M_t,centers = K_hat,nstart = 20)$centers
      z_hat <- update_Z(X, Y, beta_hat)
      out<- c(iters,SCORE_hat(X,Y,beta0,beta_hat,z0,z_hat))
    }
    if(save_estimate){
      final_beta <- if (save_history) beta_hat[iters,,] else beta_hat
      return(list(
        beta_hat=final_beta,
        z_hat=z_hat,
        l=l_calculate(X,Y,final_beta,z_hat)
      ))
    }
    return(list(
      out = out,
      t = iters
    ))
  }, error = function(e) {
    return(NULL)
  })
}

Oracle_alg<-function(X,Y,z0,beta0){
  M<-dim(X)[1];n<-dim(X)[2];p<-dim(X)[3];K<-dim(beta0)[1]
  beta_hat<-matrix(NA,ncol = p,nrow = K)
  for(k in 1:K){
    sigma<-array(0,dim = c(p, p))
    su<-0
    for(m in 1:M){
      w <- z0[m,,k]
      sigma <- sigma + crossprod(X[m,,], X[m,,] * w)
      su <- su + crossprod(X[m,,], Y[m,] * w)
    }
    beta_hat[k,] <- solve(sigma, su)
  }
  z_hat <- update_Z(X, Y, beta_hat)
  out<-SCORE_hat(X,Y,beta0,beta_hat,z0,z_hat)
  return(out)
}

plot_mse <- function(iter_my, iter_non, iter_sin,iter_ave,oracle, name, iters, na.rm = FALSE) {
  df_index <- data.frame(
    Index = c(1:iters, 1:iters, 1:iters, 1:iters),
    Values = c(iter_my, iter_non, iter_sin,iter_ave),
    algo = c(rep("E-FICR", iters), rep("CICR", iters),rep("SC", iters),rep("AA", iters))
  )
  df_index$algo <- factor(df_index$algo, levels = c("CICR", "E-FICR", "AA","SC"))
  plot <- ggplot(df_index, aes(x = Index, y = Values, color = algo, shape = algo)) +
    geom_line(linewidth = .7, na.rm = na.rm) +
    geom_point(size = 1.1, na.rm = na.rm) +
    geom_hline(yintercept = oracle, color = "#fbb45d", linetype = "dashed", linewidth = 1) +
    labs(x="t",y = name, color = "Algorithm:", shape = "Algorithm:") +  
    scale_color_manual(values = c("#699ed4", "#ef8183", '#7A70B5','#81B3A9')) + 
    scale_shape_manual(values = c(17, 15, 20,18)) +  
    theme_minimal(base_size = 14, base_family = "sans")+
    guides(color = guide_legend(override.aes = list(size = 5)))+
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      axis.title.x = element_blank(),
      legend.text = element_text(size = 14),
      legend.key.size = unit(1.5, "cm")
    ) +
    guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
  return(plot)
}

PPFL_alg <- function(X, Y,beta_init, z0, beta0, K_hat, iters, lr_Q = 0.01, lr_c = 0.05,lambda=0.001,save_history = FALSE) {
  tryCatch({
    M <- dim(X)[1]; n <- dim(X)[2]; p <- dim(X)[3]
    Q <- t(beta_init)
    C <- matrix(1/K_hat, nrow = M, ncol = K_hat)
    if (save_history) {out <- NULL}
    
    for (t in 1:iters) {
      grad_Q_total <- matrix(0, nrow = p, ncol = K_hat)
      
      C_mean <- colMeans(C)
      
      for (m in 1:M) {
        X_m <- matrix(X[m,,], nrow = n, ncol = p)
        Y_m <- matrix(Y[m,], nrow = n, ncol = 1)
        
        for (inner in 1:3) {
          c_m <- C[m, ]
          
          pred <- (X_m %*% Q) %*% c_m
          grad_c_data <- - t(X_m %*% Q) %*% (Y_m - pred) / n
          
          grad_c_penalty <- 2 * lambda * (c_m - C_mean)
          
          G_c_i <- as.vector(grad_c_data + grad_c_penalty)
          
          numerator <- c_m * exp(-lr_c * G_c_i)
          C[m, ] <- numerator / sum(numerator)
        }
        
        c_m <- C[m, ]
        pred <- (X_m %*% Q) %*% c_m
        grad_Q_m <- - (t(X_m) %*% (Y_m - pred) %*% t(c_m)) / n
        grad_Q_total <- grad_Q_total + grad_Q_m
      }
      
      Q <- Q - lr_Q * (grad_Q_total / M)
      
      if (save_history) {
        beta_hat_K <- t(Q) 
        z_hat <- update_Z(X, Y, beta_hat_K)
        score <- SCORE_hat(X, Y, beta0, beta_hat_K, z0, z_hat)
        out <- rbind(out, c(t, score))
      } else {
        beta_hat_K <- t(Q) 
        z_hat <- update_Z(X, Y, beta_hat_K)
        score <- SCORE_hat(X, Y, beta0, beta_hat_K, z0, z_hat)
        out <- c(t, score)
      }
    }
    return(list(out = out, t = t, beta_hat_K = beta_hat_K))
  }, error = function(e) {
    message("PPFL Algorithm Error: ", e)
    return(NULL)
  })
}

S1_main <- function(beta0, num_init = 5, seed = 1, iters = 40,
                    K_hat = 3, n = 300, M = 10, shred = 0.0001,
                    balance = TRUE, rho = 0.3) {
  K_true <- dim(beta0)[1]
  p <- dim(beta0)[2]
  
  set.seed(seed)
  mean_vector <- rep(0, p)
  cov_matrix <- rho^abs(outer(1:p, 1:p, "-"))
  
  X <- generate_X(
    M = M,
    n = n,
    beta = beta0,
    Sigma = cov_matrix
  )
  Y <- matrix(0, nrow = M, ncol = n)
  ep <- matrix(rnorm(M * n, sd = 1), nrow = M, ncol = n)
  
  z0 <- array(0, dim = c(M, n, K_true))
  if (balance == TRUE) {
    for (m in 1:M) {
      for (k in 1:K_true) {
        group_index <- ((k - 1) * n / K_true + 1):(k * n / K_true)
        Y[m, group_index] <-
          X[m, group_index, ] %*% beta0[k, ] + ep[m, group_index]
        z0[m, group_index, k] <- 1
      }
    }
  }
  
  E_FICR_out <- NULL
  CICR_out <- NULL
  SC_out <- NULL
  AA_out <- NULL
  SC_out2 <- NULL
  beta_init <- NULL
  
  for (ini in 1:num_init) {
    beta_init <- rbind(
      beta_init,
      INIT_beta_1K(
        X, Y, K_hat,
        seed = 1000 * ini + seed,
        lam = 0.001,
        select_m = ini
      )
    )
  }
  
  e <- AA_alg(
    X, Y,
    beta_init = beta_init,
    z0 = z0,
    beta0 = beta0,
    iters = iters,
    shred = shred,
    num_init = num_init,
    seed = seed,
    save_history = TRUE
  )
  AA_out <- cbind(seed, K_hat, e$t, e$out, 0)
  
  for (ini in 1:num_init) {
    init_index <- ((ini - 1) * K_hat + 1):(ini * K_hat)
    
    a <- E_FICR_alg(
      X, Y,
      beta_init = beta_init[init_index, ],
      z0 = z0,
      beta0 = beta0,
      iters = iters,
      shred = shred,
      save_history = TRUE
    )
    b <- CICR_alg(
      X, Y,
      beta_init = beta_init[init_index, ],
      z0 = z0,
      beta0 = beta0,
      iters = iters,
      shred = shred,
      save_history = TRUE
    )
    d <- SC_alg(
      X, Y,
      beta_init = beta_init[init_index, ],
      z0 = z0,
      beta0 = beta0,
      iters = iters,
      shred = shred,
      single_m = 1,
      save_history = TRUE
    )
    f <- SC_alg(
      X, Y,
      beta_init = beta_init[init_index, ],
      z0 = z0,
      beta0 = beta0,
      iters = iters,
      shred = shred,
      single_m = M,
      save_history = TRUE
    )
    
    a <- cbind(seed, K_hat, a$t, a$out, ini)
    b <- cbind(seed, K_hat, b$t, b$out, ini)
    d <- cbind(seed, K_hat, d$t, d$out, ini)
    f <- cbind(seed, K_hat, f$t, f$out, ini)
    
    E_FICR_out <- get_best_res(E_FICR_out, a)
    CICR_out <- get_best_res(CICR_out, b)
    SC_out <- get_best_res(SC_out, d)
    SC_out2 <- get_best_res(SC_out2, f)
  }
  
  res_col <- c(
    "seed", "k_hat", "t_max", "t", "mse", "acc",
    "mse_g1", "mse_g2", "mse_g3", "mse_g4",
    "acc_g1", "acc_g2", "acc_g3", "acc_g4",
    "mse2", "ARI", "l", "ell", "init"
  )
  for (obj_name in c(
    "E_FICR_out", "CICR_out", "SC_out", "SC_out2", "AA_out"
  )) {
    temp_obj <- get(obj_name)
    if (!is.null(temp_obj)) {
      colnames(temp_obj) <- res_col
      assign(obj_name, temp_obj)
    }
  }
  
  Oracle_out <- Oracle_alg(X, Y, z0, beta0 = beta0)
  if (K_hat == K_true) {
    Oracle_out <- cbind(seed, K_hat, t(Oracle_out))
    colnames(Oracle_out) <- c(
      "seed", "k_hat", "mse", "acc",
      "mse_g1", "mse_g2", "mse_g3", "mse_g4",
      "acc_g1", "acc_g2", "acc_g3", "acc_g4",
      "mse2", "ARI", "l", "ell"
    )
  } else {
    Oracle_out <- NULL
  }
  
  return(list(
    Oracle_out = Oracle_out,
    CICR_out = CICR_out,
    E_FICR_out = E_FICR_out,
    AA_out = AA_out,
    SC_out = SC_out,
    SC_out2 = SC_out2
  ))
}

S12_main<-function(beta0,num_init=5,seed=1,iters=40,K_hat=3,n=300,M=10,shred=0.0001,balance=TRUE,rho=0.3,run_sc2=TRUE,run_ppfl=TRUE){
  K_true<-dim(beta0)[1];p<-dim(beta0)[2]
  set.seed(seed)
  mean_vector <- rep(0,p)
  cov_matrix <- rho^abs(outer(1:p, 1:p, "-"))
  
  X <- generate_X(
    M = M,
    n = n,
    beta = beta0,
    Sigma = cov_matrix
  )
  Y<- matrix(0, nrow = M, ncol = n)
  ep<- matrix(rnorm(M*n,sd=1), nrow = M, ncol = n)
  
  z0<-array(0 ,dim = c(M, n, K_true))
  if(balance == TRUE){
    for (m in 1:M) {
      for(k in 1:K_true){
        Y[m,((k-1)*n/K_true+1):(k*n/K_true)]<-X[m,((k-1)*n/K_true+1):(k*n/K_true),]%*%beta0[k,]+ep[m,((k-1)*n/K_true+1):(k*n/K_true)]
        z0[m,((k-1)*n/K_true+1):(k*n/K_true),k]<-1
      }
    }
  }

  E_FICR_out<-NULL;CICR_out<-NULL;SC_out<-NULL;AA_out<-NULL;SC_out2<-NULL
  PPFL_out <- NULL
  beta_init<-NULL
  for (ini in 1:num_init) {
    beta_init<-rbind(beta_init, INIT_beta_1K(X,Y,K_hat,seed=1000*ini+seed,lam=0.001,select_m=ini))
  }
  
  e<-AA_alg(X,Y,beta_init=beta_init,z0=z0,beta0=beta0,iters=iters,shred=shred,num_init=num_init,seed=seed)
  AA_out<-t(as.matrix(c(seed,K_hat,e$t,e$out,0)))
  
  for (ini in 1:num_init) {
    a<-E_FICR_alg(X,Y,beta_init=beta_init[((ini-1)*K_hat+1):(ini*K_hat),],z0=z0,beta0=beta0,iters=iters,shred=shred)
    b<-CICR_alg(X,Y,beta_init=beta_init[((ini-1)*K_hat+1):(ini*K_hat),],z0=z0,beta0=beta0,iters=iters,shred=shred)
    d<-SC_alg(X,Y,beta_init=beta_init[((ini-1)*K_hat+1):(ini*K_hat),],z0=z0,beta0=beta0,iters=iters,shred=shred,single_m = 1)
    if (run_sc2) f<-SC_alg(X,Y,beta_init=beta_init[((ini-1)*K_hat+1):(ini*K_hat),],z0=z0,beta0=beta0,iters=iters,shred=shred,single_m = M)
    a<-t(as.matrix(c(seed,K_hat,a$t,a$out,ini)))
    b<-t(as.matrix(c(seed,K_hat,b$t,b$out,ini)))
    d<-t(as.matrix(c(seed,K_hat,d$t,d$out,ini)))
    if (run_sc2) f<-t(as.matrix(c(seed,K_hat,f$t,f$out,ini)))
    if (run_ppfl) {
      p_ppfl <- PPFL_alg(X, Y,beta_init=beta_init[((ini-1)*K_hat+1):(ini*K_hat),], z0=z0, beta0=beta0, K_hat=K_hat, iters=iters,
                         lr_Q = 0.05, lr_c = 0.05)
      p_ppfl_res <- if(!is.null(p_ppfl)) t(as.matrix(c(seed,K_hat,p_ppfl$t,p_ppfl$out,ini))) else NULL
    }
    
    E_FICR_out  <- get_best_res_2(E_FICR_out, a)
    CICR_out   <- get_best_res_2(CICR_out, b)
    SC_out  <- get_best_res_2(SC_out, d)
    if (run_sc2) SC_out2 <- get_best_res_2(SC_out2, f)
    if (run_ppfl) PPFL_out <- get_best_res_2(PPFL_out, p_ppfl_res)
  }
  res_col <- c("seed","k_hat","t_max","t","mse","acc",
               "mse_g1","mse_g2","mse_g3","mse_g4",
               "acc_g1","acc_g2","acc_g3","acc_g4",
               "mse2","ARI","l","ell","init")
  for (obj_name in c("E_FICR_out","CICR_out","SC_out","SC_out2","AA_out","PPFL_out")) {
    temp_obj <- get(obj_name)
    if (!is.null(temp_obj)) {
      colnames(temp_obj) <- res_col
      assign(obj_name, temp_obj)
    }
  }
  Oracle_out<-Oracle_alg(X,Y,z0,beta0=beta0)
  if(K_hat==K_true){
    Oracle_out<-cbind(seed,K_hat,t(Oracle_out));colnames(Oracle_out)<-c("seed","k_hat","mse","acc","mse_g1","mse_g2","mse_g3","mse_g4","acc_g1","acc_g2","acc_g3","acc_g4","mse2","ARI","l","ell")
  } else Oracle_out<-NULL
  ans <- list(Oracle_out=Oracle_out,CICR_out=CICR_out,E_FICR_out=E_FICR_out,
              AA_out=AA_out,SC_out=SC_out,SC_out2=SC_out2,PPFL_out=PPFL_out)
  if (!run_sc2) ans$SC_out2 <- NULL
  if (!run_ppfl) ans$PPFL_out <- NULL
  return(ans)
}

unba_result_columns <- function(K = 4L, oracle = FALSE) {
  group_mse_cols <- paste0("mse_g", seq_len(K))
  group_acc_cols <- paste0("acc_g", seq_len(K))
  
  if (oracle) {
    c(
      "seed", "k_hat", "mse", "acc",
      group_mse_cols, group_acc_cols,
      "mse2", "ARI", "l", "ell"
    )
  } else {
    c(
      "seed", "k_hat", "t_max", "t", "mse", "acc",
      group_mse_cols, group_acc_cols,
      "mse2", "ARI", "l", "ell", "init"
    )
  }
}

unba_na_result <- function(seed, K, oracle = FALSE) {
  columns <- unba_result_columns(K, oracle)
  result <- matrix(
    NA_real_, nrow = 1L, ncol = length(columns),
    dimnames = list(NULL, columns)
  )
  result[1L, c("seed", "k_hat")] <- c(seed, K)
  result
}

S2_main<-function(beta0,num_init=5,seed=1,iters=40,K_hat=3,n=300,M=10,shred=0.0001,balance=TRUE,rho=0.3){
  K_true<-dim(beta0)[1];p<-dim(beta0)[2]
  set.seed(seed)
  mean_vector <- rep(0,p)
  cov_matrix <- rho^abs(outer(1:p, 1:p, "-"))
  data <- mvrnorm(n=M*n, mu = mean_vector, Sigma = cov_matrix)
  
  X<- array(data, dim = c(M, n, p))
  Y<- matrix(0, nrow = M, ncol = n)
  ep<- matrix(rnorm(M*n,sd=1), nrow = M, ncol = n)
  
  z0<-array(0 ,dim = c(M, n, K_true))
  if(balance == TRUE){
    for (m in 1:M) {
      for(k in 1:K_true){
        Y[m,((k-1)*n/K_true+1):(k*n/K_true)]<-X[m,((k-1)*n/K_true+1):(k*n/K_true),]%*%beta0[k,]+ep[m,((k-1)*n/K_true+1):(k*n/K_true)]
        z0[m,((k-1)*n/K_true+1):(k*n/K_true),k]<-1
      }
    }
  } else { 
    for (m in 1:floor(M/2)) {
      radio<-c(0,cumsum(c(0.45,0.45,0,0.1))) 
      for(k in 1:K_true){
        if(round(n*radio[k]+1)-round(n*radio[k+1])<=0){
          X[m,round(n*radio[k]+1):round(n*radio[k+1]),]<-X[m,round(n*radio[k]+1):round(n*radio[k+1]),]*sqrt(k)
          Y[m,round(n*radio[k]+1):round(n*radio[k+1])]<-X[m,round(n*radio[k]+1):round(n*radio[k+1]),]%*%beta0[k,]+ep[m,round(n*radio[k]+1):round(n*radio[k+1])]
          z0[m,round((n*radio[k]+1):round(n*radio[k+1])),k]<-1
        }
      }
    }
    for (m in floor(M/2+1):M) {
      radio<-c(0,cumsum(c(0.45,0.45,0.1,0)))
      for(k in 1:K_true){
        if(round(n*radio[k]+1)-round(n*radio[k+1])<=0){
          X[m,round(n*radio[k]+1):round(n*radio[k+1]),]<-X[m,round(n*radio[k]+1):round(n*radio[k+1]),]*sqrt(k)
          Y[m,round(n*radio[k]+1):round(n*radio[k+1])]<-X[m,round(n*radio[k]+1):round(n*radio[k+1]),]%*%beta0[k,]+ep[m,round(n*radio[k]+1):round(n*radio[k+1])]
          z0[m,round((n*radio[k]+1):round(n*radio[k+1])),k]<-1
        }
      }
    }
  }

  errors <- data.frame(seed = integer(), method = character(), init = integer(),
                       message = character())
  record_error <- function(method, ini, message) {
    errors <<- rbind(errors, data.frame(seed = seed, method = method,
                                        init = ini, message = message))
  }
  attempt <- function(method, ini, expr) {
    ans <- tryCatch(expr, error = function(e) {
      record_error(method, ini, conditionMessage(e)); NULL
    })
    if (is.null(ans)) record_error(method, ini, "No result (algorithm returned NULL)")
    ans
  }
  pack <- function(a, ini) {
    if (is.null(a) || length(a$out) != 15L || any(!is.finite(a$out))) return(NULL)
    matrix(c(seed, K_hat, a$t, a$out, ini), nrow = 1,
           dimnames = list(NULL, unba_result_columns()))
  }
  best_loss <- function(old, candidate) {
    if (is.null(candidate)) return(old)
    if (is.null(old) || candidate[1, "l"] < old[1, "l"]) candidate else old
  }
  beta_init<-NULL
  for (ini in 1:num_init) {
    beta_init<-rbind(beta_init, INIT_beta_M(X,Y,K_hat,seed=seed+1000*ini,lam=0.001))
    beta_init<-rbind(beta_init, INIT_random(X,Y,beta0=beta0,K_hat=K_hat,seed=1000*ini+seed))
  }
  e <- attempt("AA", NA_integer_, AA_alg(X, Y, beta_init, z0, beta0,
                                         iters, shred, num_init = 2*num_init, seed = seed))
  AA_out <- pack(e, 2*num_init)
  best <- setNames(vector("list", 5L), c("E_FICR", "KR", "SMA", "SMA2", "PPFL"))
  for (ini in seq_len(2*num_init)) {
    b0 <- beta_init[(ini - 1L) * K_hat + seq_len(K_hat), , drop = FALSE]
    fits <- list(
      E_FICR = attempt("E_FICR", ini, E_FICR_alg(X, Y, b0, z0, beta0, iters, shred)),
      KR = attempt("KR", ini, CICR_alg(X, Y, b0, z0, beta0, iters, shred))
      , SMA = attempt("SMA", ini, SC_alg(X, Y, b0, z0, beta0, iters, shred, single_m = 1))
      , SMA2 = attempt("SMA2", ini, SC_alg(X, Y, b0, z0, beta0, iters, shred, single_m = M))
      , PPFL = attempt("PPFL", ini, PPFL_alg(X, Y, b0, z0, beta0, K_hat, iters,
                                             lr_Q = .05, lr_c = .05))
    )
    for (method in names(fits)) {
      candidate <- pack(fits[[method]], ini)
      if (!is.null(fits[[method]]) && is.null(candidate))
        record_error(method, ini, "Invalid or non-finite output")
      best[method] <- list(best_loss(best[[method]], candidate))
    }
  }
  oracle <- attempt("Oracle", NA_integer_, Oracle_alg(X, Y, z0, beta0))
  Oracle_out <- if (length(oracle) == 14L && all(is.finite(oracle)))
    matrix(c(seed, K_hat, oracle), nrow = 1L,
           dimnames = list(NULL, unba_result_columns(oracle = TRUE))) else
             unba_na_result(seed, K_hat, oracle = TRUE)
  fallback <- function(x) if (is.null(x)) unba_na_result(seed, K_hat) else x
  list(Oracle_out = Oracle_out, KR_out = fallback(best$KR),
       E_FICR_out = fallback(best$E_FICR), AA_out = fallback(AA_out),
       SMA_out = fallback(best$SMA), SMA_out2 = fallback(best$SMA2),
       PPFL_out = fallback(best$PPFL), errors = errors)
}


# CVA returns normalized instability s_n = 1 - RI; smaller is better.
cva_M <- function(X, Y, beta_init, z0, beta0, K_hat,
                  fit_alg = CICR_alg, seed = 1, shred, iters) {
  tryCatch({
    M <- dim(X)[1]
    num_init <- nrow(beta_init) / K_hat
    if (M < 3L) stop("cva_M requires at least three clients")
    
    set.seed(seed)
    fold <- sample(rep_len(1:3, M))
    split_data <- function(ids) {
      list(
        X = X[ids, , , drop = FALSE],
        Y = Y[ids, , drop = FALSE],
        z0 = z0[ids, , , drop = FALSE]
      )
    }
    train1 <- split_data(fold == 1L)
    train2 <- split_data(fold == 2L)
    test <- split_data(fold == 3L)
    
    best_fit <- function(dat) {
      best <- NULL
      for (ini in seq_len(num_init)) {
        rows <- ((ini - 1L) * K_hat + 1L):(ini * K_hat)
        fit <- fit_alg(
          dat$X, dat$Y, beta_init = beta_init[rows, , drop = FALSE],
          z0 = dat$z0, beta0 = beta0, shred = shred, iters = iters,
          save_estimate = TRUE
        )
        if (!is.null(fit) && is.finite(fit$l) &&
            (is.null(best) || fit$l < best$l)) best <- fit
      }
      best
    }
    
    fit1 <- best_fit(train1)
    fit2 <- best_fit(train2)
    if (is.null(fit1) || is.null(fit2)) return(NULL)
    
    z_test1 <- update_Z(test$X, test$Y, fit1$beta_hat)
    z_test2 <- update_Z(test$X, test$Y, fit2$beta_hat)
    1 - RI_calculate(z_test1, z_test2)
  }, error = function(e) NULL)
}

# For SC, split observations within each client.
cva_n <- function(X, Y, beta_init, z0, beta0, K_hat,
                  fit_alg = SC_alg, seed = 1, shred, iters,
                  single_m = 1L) {
  tryCatch({
    M <- dim(X)[1]
    n <- dim(X)[2]
    p <- dim(X)[3]
    K_true <- nrow(beta0)
    num_init <- nrow(beta_init) / K_hat
    if (n %% 3L != 0L) stop("cva_n requires n to be divisible by three")
    if (single_m < 1L || single_m > M) stop("single_m is outside 1, ..., M")
    
    set.seed(seed)
    fold <- t(replicate(M, sample(rep(1:3, each = n / 3L))))
    n_fold <- n / 3L
    make_split <- function(which_fold) {
      X_part <- array(0, dim = c(M, n_fold, p))
      Y_part <- matrix(0, nrow = M, ncol = n_fold)
      z_part <- array(0, dim = c(M, n_fold, K_true))
      for (m in seq_len(M)) {
        ids <- fold[m, ] == which_fold
        X_part[m, , ] <- X[m, ids, ]
        Y_part[m, ] <- Y[m, ids]
        z_part[m, , ] <- z0[m, ids, ]
      }
      list(X = X_part, Y = Y_part, z0 = z_part)
    }
    train1 <- make_split(1L)
    train2 <- make_split(2L)
    test <- make_split(3L)
    
    best_fit <- function(dat) {
      best <- NULL
      for (ini in seq_len(num_init)) {
        rows <- ((ini - 1L) * K_hat + 1L):(ini * K_hat)
        fit <- fit_alg(
          dat$X, dat$Y, beta_init = beta_init[rows, , drop = FALSE],
          z0 = dat$z0, beta0 = beta0, shred = shred, iters = iters,
          single_m = single_m, save_estimate = TRUE
        )
        if (!is.null(fit) && is.finite(fit$l) &&
            (is.null(best) || fit$l < best$l)) best <- fit
      }
      best
    }
    
    fit1 <- best_fit(train1)
    fit2 <- best_fit(train2)
    if (is.null(fit1) || is.null(fit2)) return(NULL)
    
    z_test1 <- update_Z(test$X, test$Y, fit1$beta_hat)
    z_test2 <- update_Z(test$X, test$Y, fit2$beta_hat)
    1 - RI_calculate(z_test1, z_test2)
  }, error = function(e) NULL)
}

cva_M_forAA_alg <- function(X, Y, beta_init, z0, beta0, K_hat,
                            fit_alg = AA_alg, seed = 1, shred, iters) {
  tryCatch({
    M <- dim(X)[1]
    num_init <- nrow(beta_init) / K_hat
    if (M < 3L) stop("cva_M_forAA_alg requires at least three clients")
    
    set.seed(seed)
    fold <- sample(rep_len(1:3, M))
    split_data <- function(ids) {
      list(
        X = X[ids, , , drop = FALSE],
        Y = Y[ids, , drop = FALSE],
        z0 = z0[ids, , , drop = FALSE]
      )
    }
    train1 <- split_data(fold == 1L)
    train2 <- split_data(fold == 2L)
    test <- split_data(fold == 3L)
    
    fit_once <- function(dat) {
      fit_alg(
        dat$X, dat$Y, beta_init = beta_init, z0 = dat$z0,
        beta0 = beta0, shred = shred, iters = iters,
        num_init = num_init, seed = seed, save_estimate = TRUE
      )
    }
    fit1 <- fit_once(train1)
    fit2 <- fit_once(train2)
    if (is.null(fit1) || is.null(fit2)) return(NULL)
    
    z_test1 <- update_Z(test$X, test$Y, fit1$beta_hat)
    z_test2 <- update_Z(test$X, test$Y, fit2$beta_hat)
    1 - RI_calculate(z_test1, z_test2)
  }, error = function(e) NULL)
}

S3_main <- function(beta0, num_init = 5, seed = 1, iters = 40,
                    n = 300, M = 10, shred = 0.0001, balance = TRUE,
                    rho = 0.3, K_choose = c(2, 3, 4),
                    evaluate_selected = FALSE) {
  K_true <- nrow(beta0)
  p <- ncol(beta0)
  if (!is.numeric(K_choose) || !length(K_choose) ||
      any(!is.finite(K_choose)) || any(K_choose != floor(K_choose)) ||
      anyDuplicated(K_choose) || any(K_choose < 2L) ||
      any(K_choose > min(n, 100L))) {
    stop("K_choose must contain distinct integers between 2 and min(n, 100)")
  }
  K_choose <- as.integer(K_choose)
  if (num_init < 1L || num_init > M) stop("num_init must be between 1 and M")
  if (M < 3L || n %% 3L != 0L) stop("CVA requires M >= 3 and n divisible by 3")
  
  set.seed(seed)
  cov_matrix <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
  X <- array(
    MASS::mvrnorm(n = M * n, mu = rep(0, p), Sigma = cov_matrix),
    dim = c(M, n, p)
  )
  Y <- matrix(0, nrow = M, ncol = n)
  ep <- matrix(rnorm(M * n), nrow = M, ncol = n)
  z0 <- array(0, dim = c(M, n, K_true))
  
  if (isTRUE(balance)) {
    if (n %% K_true != 0L) stop("Balanced design requires n divisible by K_true")
    counts <- matrix(n / K_true, nrow = M, ncol = K_true)
  } else {
    if (!is.matrix(balance) || !all(dim(balance) == c(M, K_true))) {
      stop("balance must be TRUE or an M by K_true matrix of proportions")
    }
    if (any(balance < 0) || any(abs(rowSums(balance) - 1) > 1e-8)) {
      stop("Each row of balance must contain nonnegative proportions summing to one")
    }
    raw_counts <- n * balance
    if (any(abs(raw_counts - round(raw_counts)) > 1e-8)) {
      stop("Every n * balance[m, k] must be an integer")
    }
    counts <- round(raw_counts)
  }
  
  for (m in seq_len(M)) {
    ends <- cumsum(counts[m, ])
    starts <- c(1L, head(ends, -1L) + 1L)
    for (k in seq_len(K_true)) {
      if (counts[m, k] == 0L) next
      ids <- starts[k]:ends[k]
      Y[m, ids] <- X[m, ids, ] %*% beta0[k, ] + ep[m, ids]
      z0[m, ids, k] <- 1
    }
  }
  
  scores <- list(
    KR = rep(NA_real_, length(K_choose)),
    E_FICR = rep(NA_real_, length(K_choose)),
    AA = rep(NA_real_, length(K_choose)),
    SMA = rep(NA_real_, length(K_choose))
  )
  beta_init_list <- vector("list", length(K_choose))
  names(beta_init_list) <- paste0("K", K_choose)
  or_na <- function(x) if (is.null(x)) NA_real_ else as.numeric(x)
  
  for (i_K in seq_along(K_choose)) {
    K_hat <- K_choose[i_K]
    beta_init <- do.call(rbind, lapply(seq_len(num_init), function(ini) {
      INIT_beta_1K(
        X, Y, K_hat = K_hat, seed = 1000L * ini + seed,
        lam = 0.001, select_m = ini
      )
    }))
    beta_init_list[[paste0("K", K_hat)]] <- beta_init
    
    scores$KR[i_K] <- or_na(cva_M(
      X, Y, beta_init, z0, beta0, K_hat,
      fit_alg = CICR_alg, seed = seed, shred = shred, iters = iters
    ))
    scores$E_FICR[i_K] <- or_na(cva_M(
      X, Y, beta_init, z0, beta0, K_hat,
      fit_alg = E_FICR_alg, seed = seed, shred = shred, iters = iters
    ))
    scores$AA[i_K] <- or_na(cva_M_forAA_alg(
      X, Y, beta_init, z0, beta0, K_hat,
      fit_alg = AA_alg, seed = seed, shred = shred, iters = iters
    ))
    scores$SMA[i_K] <- or_na(cva_n(
      X, Y, beta_init, z0, beta0, K_hat,
      fit_alg = SC_alg, seed = seed, shred = shred, iters = iters
    ))
  }
  
  score_names <- paste0("K", K_choose)
  cva_result <- lapply(
    scores,
    function(x) c(seed = seed, stats::setNames(x, score_names))
  )
  if (!isTRUE(evaluate_selected)) return(cva_result)
  
  # MSE2 and ARI remain valid when the selected K differs from the true K.
  evaluate_fit <- function(beta_hat, z_hat) {
    mse2 <- 0
    true_labels <- integer(M * n)
    fitted_labels <- integer(M * n)
    for (m in seq_len(M)) {
      mse2 <- mse2 + sum(
        (z0[m, , ] %*% beta0 - z_hat[m, , ] %*% beta_hat)^2
      )
      ids <- ((m - 1L) * n + 1L):(m * n)
      true_labels[ids] <- max.col(z0[m, , ], ties.method = "first")
      fitted_labels[ids] <- max.col(z_hat[m, , ], ties.method = "first")
    }
    c(
      mse2 = mse2 / (M * n),
      ARI = mclust::adjustedRandIndex(true_labels, fitted_labels),
      l = l_calculate(X, Y, beta_hat, z_hat)
    )
  }
  
  best_full_fit <- function(method, K_hat) {
    beta_init <- beta_init_list[[paste0("K", K_hat)]]
    if (method == "AA") {
      return(AA_alg(
        X, Y, beta_init, z0, beta0, iters, shred,
        num_init = num_init, seed = seed, save_estimate = TRUE
      ))
    }
    
    best <- NULL
    for (ini in seq_len(num_init)) {
      rows <- ((ini - 1L) * K_hat + 1L):(ini * K_hat)
      b0 <- beta_init[rows, , drop = FALSE]
      fit <- switch(
        method,
        CICR = CICR_alg(
          X, Y, b0, z0, beta0, iters, shred, save_estimate = TRUE
        ),
        E_FICR = E_FICR_alg(
          X, Y, b0, z0, beta0, iters, shred, save_estimate = TRUE
        ),
        SC = SC_alg(
          X, Y, b0, z0, beta0, iters, shred,
          single_m = 1L, save_estimate = TRUE
        )
      )
      if (!is.null(fit) && is.finite(fit$l) &&
          (is.null(best) || fit$l < best$l)) {
        best <- fit
      }
    }
    best
  }
  
  method_map <- c(CICR = "KR", E_FICR = "E_FICR", AA = "AA", SC = "SMA")
  selected_k <- vapply(method_map, function(score_name) {
    method_scores <- scores[[score_name]]
    # Exclude incomplete candidate sets from K selection.
    if (!all(is.finite(method_scores))) return(NA_integer_)
    K_choose[which.min(method_scores)]
  }, integer(1))
  
  selected_metrics <- do.call(rbind, lapply(names(method_map), function(method) {
    K_hat <- selected_k[[method]]
    fit <- if (is.na(K_hat)) NULL else best_full_fit(method, K_hat)
    metrics <- if (is.null(fit)) {
      c(mse2 = NA_real_, ARI = NA_real_, l = NA_real_)
    } else {
      evaluate_fit(fit$beta_hat, fit$z_hat)
    }
    data.frame(
      seed = seed, method = method, k_hat = K_hat,
      mse2 = metrics[["mse2"]], ARI = metrics[["ARI"]],
      l = metrics[["l"]], row.names = NULL
    )
  }))
  
  oracle_beta <- matrix(NA_real_, nrow = K_true, ncol = p)
  for (k in seq_len(K_true)) {
    sigma <- matrix(0, nrow = p, ncol = p)
    xy <- numeric(p)
    for (m in seq_len(M)) {
      w <- z0[m, , k]
      sigma <- sigma + crossprod(X[m, , ], X[m, , ] * w)
      xy <- xy + crossprod(X[m, , ], Y[m, ] * w)
    }
    oracle_beta[k, ] <- solve(sigma, xy)
  }
  oracle_z <- update_Z(X, Y, oracle_beta)
  oracle_metrics <- evaluate_fit(oracle_beta, oracle_z)
  selected_metrics <- rbind(
    data.frame(
      seed = seed, method = "Oracle", k_hat = K_true,
      mse2 = oracle_metrics[["mse2"]], ARI = oracle_metrics[["ARI"]],
      l = oracle_metrics[["l"]], row.names = NULL
    ),
    selected_metrics
  )
  
  list(cva = cva_result, selected_metrics = selected_metrics)
}
