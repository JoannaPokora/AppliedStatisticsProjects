
library(igraph)
library(cvCovEst)
library(grpSLOPE)
library(glasso)
library(PCGLASSO)
library(ggplot2)
library(tidyr)
library(dplyr)

set.seed(12)

# 1

f_grad <- function(vals){
  x <- vals[1]
  y <- vals[2]
  z <- vals[3]
  
  c(
    4*(x-1)^3 - 12*(x-1)^2 + 8*(x-1) + y + 1,
    4*(y+1)^3 - 12*(y+1)^2 + 8*(y+1) + x + 1 - z,
    4*(z-2)^3 - 12*(z-2)^2 + 8*(z-2) - y - 1
  )
}

f <- function(vals){
  x <- vals[1]
  y <- vals[2]
  z <- vals[3]
  
  (x-1)^4 - 4*(x-1)^3 + 4*(x-1)^2 + (y+1)^4 - 4*(y+1)^3 + 4*(y+1)^2 +
    (z-2)^4 - 4*(z-2)^3 + 4*(z-2)^2 + (x-1)*(y+1) - (y+1)*(z-2)
}

# Task 8

eta_gd <- 0.05
end_gd <- 0.01

f_star <- 0

iter_gd <- 0
prev_vals_gd <- c(1.2, -0.8, 2.1)

while(f(prev_vals_gd) - f_star >= end_gd || iter_gd == 0){
  iter_gd <- iter_gd + 1
  new_vals <- prev_vals_gd - eta_gd*f_grad(prev_vals_gd)
  prev_vals_gd <- new_vals
}

iter_gd

prev_vals_gd

f(prev_vals_gd)

# Task 10

eta_cgd <- 0.05
end_cgd <- 0.01

iter_cgd <- 0
vals_cgd <- c(1.2, -0.8, 2.1)

while(f(vals_cgd) - f_star >= end_cgd || iter_cgd == 0){
  iter_cgd <- iter_cgd + 1
  for(coord in 1:3){
    vals_cgd[coord] <- vals_cgd[coord] - eta_cgd*f_grad(vals_cgd)[coord] 
  }
}

iter_cgd

vals_cgd

f(vals_cgd)

# Task 11

f_sec_grad <- function(vals){
  x <- vals[1]
  y <- vals[2]
  z <- vals[3]
  
  matrix(
    c(
      12*x^2 - 48*x + 44, 1, 0,
      1, 12*y^2 - 4, -1,
      0, -1, 12*z^2 - 72*z + 104
    ),
    byrow = TRUE,
    nrow = 3
  )
}

direction <- function(vals){
  grad <- f_grad(vals)
  sec_grad <- f_sec_grad(vals)
  
  -solve(sec_grad) %*% grad
}
  
end_n <- 0.01

iter_n <- 0
prev_vals_n <- c(1.2, -0.8, 2.1)

while(f(prev_vals_n) - f_star >= end_n || iter_n == 0){
  iter_n <- iter_n + 1
  new_vals_n <- prev_vals_n + direction(prev_vals_n)
  prev_vals_n <- new_vals_n
}

iter_n

prev_vals_n

f(vals_cgd)

#################

no_sum_df <- data.frame(
  "Method" = c("gradient descent", "cyclic coordinate gradient descent",
               "Newton's"),
  "Number of iterations" = c(iter_gd, iter_cgd, iter_n),
  "Final iterate" = c(
    paste0("(", paste(round(prev_vals_gd, 3), collapse = ", "), ")"),
    paste0("(", paste(round(vals_cgd, 3), collapse = ", "), ")"),
    paste0("(", paste(round(prev_vals_n, 3), collapse = ", "), ")")
  ),
  "Final value" = round(c(f(prev_vals_gd), f(vals_cgd), f(prev_vals_n)),3),
  check.names = FALSE
)

# Problem 2

n <- 300
p <- 100

# creating block sigma
rho <- 0.8
block <- (1 - rho) * diag(nrow = 50) + rho * matrix(1, nrow = 50, ncol = 50)
block_sigma <- diag(nrow = p)
block_sigma[1:50, 1:50] <- block
block_sigma[51:100, 51:100] <- block

block_omega <- solve(block_sigma)

block_graph <- graph_from_adjacency_matrix(block_omega != 0,
                                           mode = "undirected",
                                           diag = FALSE)

# creating omega
hub_omega <- diag(nrow = 100)
hub_omega[1,1] <- 60
hub_omega[1, 2:100] <- -0.7
hub_omega[2:100, 1] <- -0.7

# graph
hub_graph <- igraph::graph_from_adjacency_matrix(hub_omega != 0,
                                                 mode = "undirected",
                                                 diag = FALSE)

hub_sigma <- solve(hub_omega)

Sigma <- list("setting 1" = block_sigma, "setting 2" = hub_sigma)

true_omega <- list("setting 1" = block_omega, "setting 2" = hub_omega)

# Functions

simulate_X <- function(Sigma, n, p){
  lapply(Sigma, function(sigma_matr){
    MASS::mvrnorm(n, rep(0, p), sigma_matr)
  })
}

get_optim_tuning_par <- function(model_fun, matr, true_omega, lam = NA,
                                 min_par = 0.1, max_par = 1, num_par = 10,
                                 include_plot = FALSE, show_plot = FALSE){
  par_seq <- seq(min_par, max_par, length.out = num_par)
  
  est_error <- sapply(par_seq, function(par){
    par_res <- if(!all(is.na(lam))){
      scaled_lam <- par * lam
      print(scaled_lam[1])
      do.call(model_fun, list(matr, scaled_lam))
    }else{
      do.call(model_fun, list(matr, par))
    }
    
    omega_est <- switch(
      model_fun,
      "glasso" = par_res$wi,
      "gslope_new" = par_res$precision_matrix,
      "pcglasso" = par_res
    )
      
    sum((omega_est - true_omega)^2) / (nrow(omega_est)^2)
  })
  
  optim_par <- par_seq[which.min(est_error)]
  
  print(paste0("Model: ", model_fun, ", optimal parameter: ", optim_par))
  
  if(include_plot){
    est_error_df <- data.frame("estimation error" = est_error,
                               "tuning parameter" = par_seq,
                               check.names = FALSE)
    
    plt <- ggplot(est_error_df, aes(x = `tuning parameter`, y = `estimation error`)) +
      geom_line() +
      geom_point() +
      theme_light()
    
    optim_par <- list(optim_par, plt)
    
    if(show_plot == TRUE)
      print(plt)
  }
  
  optim_par
}

get_prec_matr_est <- function(X, true_omega){
  mapply(function(X_set, true_omega_set, setting){
    cov_matr <- cov(X_set)
    
    lambda_glasso <- rep(1, ncol(X_set) * (ncol(X_set) - 1) / 2)
    glasso_optim_par <- get_optim_tuning_par(
      "gslope_new", cov_matr, true_omega_set, lambda_glasso,
      min_par = ifelse(setting == "setting 1", 0.001, 0),
      max_par = ifelse(setting == "setting 1", 0.1, 0.009),
      include_plot = TRUE
    )
    
    lambda_bh <- lambdaBH(ncol(X_set), nrow(X_set))
    gslope_optim_par <- get_optim_tuning_par(
      "gslope_new", cov_matr, true_omega_set, lambda_bh,
      min_par = ifelse(setting == "setting 1", 0.12, 0),
      max_par = ifelse(setting == "setting 1", 0.3, 0.009),
      include_plot = TRUE
    )
    
    pcglasso_optim_par <- get_optim_tuning_par(
      "pcglasso", cov_matr, true_omega_set, max_par = 0.5, num_par = 5,
      include_plot = TRUE
    )
    
    list(
      estimators = list(
        "MLE" = solve(cov_matr),
        "LW" = solve(linearShrinkLWEst(X_set)),
        "gLASSO" = gslope_new(cov_matr, lambda_glasso * glasso_optim_par[[1]])$precision_matrix,
        "gSLOPE" = gslope_new(cov_matr, lambda_bh * gslope_optim_par[[1]])$precision_matrix,
        "pcgLASSO" = pcglasso(cov_matr, pcglasso_optim_par[[1]]) 
      ),
      optim_par = list(
        "MLE" = NA,
        "LW" = NA,
        "gLASSO" = glasso_optim_par[[1]],
        "gSLOPE" = gslope_optim_par[[1]],
        "pcgLASSO" = pcglasso_optim_par[[1]]
      ),
      plot = list(
        "gLASSO" = glasso_optim_par[[2]],
        "gSLOPE" = gslope_optim_par[[2]],
        "pcgLASSO" = pcglasso_optim_par[[2]]
      )
    )
  }, X, true_omega, names(X), SIMPLIFY = FALSE)
}

get_optim_params_df <- function(omega_est){
  params_lst <- lapply(omega_est, function(omega_set){
    omega_set[[2]]
  })
  
  data.frame(
    `optimal parameter` = unlist(params_lst),
    estimator = names(omega_est[[1]][[1]]),
    setting = rep(names(omega_est), each = length(names(omega_est[[1]][[1]]))),
    check.names = FALSE
  ) %>%
    pivot_wider(names_from = estimator, values_from = `optimal parameter`)%>%
    select(where(~all(!is.na(.))))
}

get_mse_df <- function(omega_est, true_omega, long = FALSE){
  mse_lst <- mapply(function(omega_set, true_omega_mat){
    lapply(omega_set[[1]], function(omega){
      1/p^2 * sum((omega - true_omega_mat)^2)
    })
  }, omega_est, true_omega, SIMPLIFY = FALSE)
  
  df <- data.frame(
    MSE = unlist(mse_lst),
    estimator = names(mse_lst[[1]]),
    setting = rep(c("Block matrix", "Hub matrix"),
                  each = length(names(mse_lst[[1]])))
  )
  
  if(long){
    df
  }else{
    df %>%
      pivot_wider(names_from = estimator, values_from = MSE)
  }
  
}

get_graph_from_prec_matr <- function(omega_est){
  lapply(omega_est, function(omega_set){
    lapply(omega_set[[1]], function(omega){
      mode <- if(isSymmetric(omega))
        "undirected"
      else
        "upper"
      
      graph_from_adjacency_matrix(omega != 0,
                                  mode = mode,
                                  diag = FALSE)
    })
  })
}

get_error_count_df <- function(omega_est, true_omega, long = FALSE){
  error_count_lst <- mapply(function(omega_set, true_omega_mat){
    lapply(omega_set[[1]], function(omega){
      ut_omega_true <- true_omega_mat[upper.tri(true_omega_mat)]
      ut_omega <- omega[upper.tri(omega)]
      type1 <- ut_omega[ut_omega_true == 0] != 0
      type2 <- ut_omega[ut_omega_true != 0] == 0
      list(
        "fp" = ifelse(length(type1) == 0, 0, sum(type1)),
        "fn" = ifelse(length(type2) == 0, 0, sum(type2))
      )
    })
  }, omega_est, true_omega, SIMPLIFY = FALSE)
  
  df <- data.frame(
    count = unlist(error_count_lst),
    `error type` = c("false positive", "false negative"),
    estimator = rep(names(omega_est[[1]][[1]]), each = 2),
    setting = rep(c("Block matrix", "Hub matrix"),
                  each = length(names(omega_est[[1]][[1]])) * 2),
    check.names = FALSE
  )
  
  if(long){
    df
  }else{
    df %>%
      pivot_wider(names_from = estimator, values_from = count)
  }
}

####################

X <- simulate_X(Sigma, n, p)

omega_est <- get_prec_matr_est(X, true_omega)

optim_params_df <- get_optim_params_df(omega_est)

plots <- lapply(omega_est, function(omega_set){
    omega_set[[3]]
})

est_mse <- get_mse_df(omega_est, true_omega)

est_graph <- get_graph_from_prec_matr(omega_est)

plot(est_graph[[1]][[1]],
     vertex.size = 10, vertex.label.cex = 0.5,
     edge.curved = FALSE,
     layout = layout_nicely)

error_count_df <- get_error_count_df(omega_est, true_omega)

# Problem 2 replications

n_rep <- 25

X_rep <- lapply(1:n_rep, function(x){
  simulate_X(Sigma, n, p)
})

omega_est_rep <- lapply(X_rep, function(X){
  get_prec_matr_est(X, true_omega)
})

optim_params_df_rep_lst <- lapply(omega_est_rep, function(omega_est){
  get_optim_params_df(omega_est)
})
optim_params_df_rep <- bind_rows(optim_params_df_rep_lst)
optim_params_rep_long <- optim_params_df_rep %>%
  pivot_longer(cols = c(gLASSO, gSLOPE, pcgLASSO),
               names_to = "Estimator",
               values_to = "Parameter") %>%
  mutate(setting = case_when(setting == "setting 1" ~ "Block matrix",
                             .default = "Hub matrix"))

params_rep_plt <- ggplot(optim_params_rep_long, aes(x = Estimator, y = Parameter)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.3, stackratio = 1.1) +
  facet_wrap(~ setting, nrow = 2) +
  theme_light()

est_mse_rep <- lapply(omega_est_rep, function(omega_est){
  get_mse_df(omega_est, true_omega, long = TRUE)
})
est_mse_rep_df <- do.call(rbind, est_mse_rep)
est_mse_rep_df$estimator <- case_when(est_mse_rep_df$estimator == "LW" ~ "inv. LW",
                                      .default = est_mse_rep_df$estimator)
est_mse_rep_df$Estimator <- factor(
  est_mse_rep_df$estimator,
  levels = c("MLE", "inv. LW", "gLASSO", "gSLOPE", "pcgLASSO")
)
est_mse_rep_plt <- ggplot(est_mse_rep_df, aes(y = MSE, x = Estimator)) +
  geom_boxplot() +
  facet_wrap(~ setting) +
  theme_light()
est_mse_rep_df_wide <- est_mse_rep_df %>%
  group_by(Estimator, setting) %>%
  summarise(`mean MSE` = mean(MSE)) %>%
  pivot_wider(names_from = Estimator, values_from = `mean MSE`)

est_graph_rep <- lapply(omega_est_rep, function(omega_est){
  get_graph_from_prec_matr(omega_est)
})

plot(est_graph_rep[[1]][[1]][[3]],
     vertex.size = 10, vertex.label.cex = 0.5,
     edge.curved = FALSE,
     layout = layout_nicely)

error_count_rep <- lapply(omega_est_rep, function(omega_est){
  get_error_count_df(omega_est, true_omega, long = TRUE)
})
error_count_rep_df <- do.call(rbind, error_count_rep)
error_count_rep_df$estimator <- factor(
  error_count_rep_df$estimator,
  levels = c("MLE", "LW", "gLASSO", "gSLOPE", "pcgLASSO")
)
error_count_rep_plt <- ggplot(error_count_rep_df, aes(y = count, x = estimator)) +
  geom_boxplot() +
  facet_grid(setting ~ `error type`) +
  theme_light()

error_count_rep_df_wide <- error_count_rep_df %>%
  group_by(setting, estimator, `error type`) %>%
  summarise(`mean count` = mean(count)) %>%
  pivot_wider(names_from = estimator, values_from = `mean count`)



