################################################################################
# 0. Working directory & basic configuration
################################################################################
# NOTE:
# - This script is shared as a partial release for peer-review transparency.
# - The original dataset and file paths are not included.
# - Please provide your own Excel file following the same column structure.

# Option A (recommended): use project root
# setwd("path/to/your/project")

# Option B: relative paths (recommended for GitHub)
data_dir    <- "data"   # users put their local data here (not tracked by git)
excel_path  <- file.path(data_dir, "input_data.xlsx")  # placeholder file name

sheet_names <- c("Cd", "Pb", "As")
target_col  <- "RGR"
output_dir  <- file.path("outputs", "xgb_outputs_cv")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# 1. Packages & theme
################################################################################
library(readxl)
library(dplyr)
library(caret)          # createDataPartition
library(xgboost)
library(ggplot2); theme_set(theme_bw(12))
library(data.table)
library(shapviz)
library(forcats)
library(tidyr)
library(tibble)

################################################################################
# 2. Data cleaning function
################################################################################
clean_sheet <- function(sheet) {
  df <- read_excel(excel_path, sheet = sheet)
  
  # Drop optional source-tag columns if they exist
  df <- df %>% select(-any_of(c("srctag", "src_tag")))
  
  # Ensure the 1st column is treated as 'Dose'
  names(df)[1] <- "Dose"
  
  # Convert specified categorical columns (if present) to numeric codes
  dum_vars <- intersect(c("CaCO3_cat", "plant_cat", "organ_cat"), names(df))
  if (length(dum_vars) > 0)
    df <- df %>% mutate(across(any_of(dum_vars), ~ as.numeric(as.character(.))))
  
  # Coerce character columns to numeric when possible; otherwise factor -> integer
  for (col in names(df)) {
    if (col %in% dum_vars || col == target_col) next
    if (is.character(df[[col]])) {
      num <- suppressWarnings(as.numeric(df[[col]]))
      df[[col]] <- if (all(is.na(num))) as.integer(factor(df[[col]])) else num
    }
    if (anyNA(df[[col]]))
      df[[col]][is.na(df[[col]])] <- median(df[[col]], na.rm = TRUE)
  }
  df
}

################################################################################
# 3. Utility: feature matrix
################################################################################
make_matrix <- function(.df) {
  .df %>%
    select(-all_of(target_col)) %>%
    mutate(across(everything(), as.numeric)) %>%
    as.matrix()
}

################################################################################
# 4. XGBoost parameter grid
################################################################################
param_grid <- expand.grid(
  max_depth        = c(3, 5),
  min_child_weight = c(1, 5),
  subsample        = c(0.7, 0.9),
  colsample_bytree = 0.8,
  eta              = c(0.05, 0.1),
  gamma            = 0,
  lambda           = c(1, 5),
  alpha            = c(0, 1)
)

################################################################################
# 5. SHAP dependence plot function
################################################################################
dep_feats <- list(
  Cd = "Labile(F0+F1)",
  Pb = "F3",
  As = "Labile(F0+F1)"
)

plot_shap_dependence <- function(
    x, shap,
    feature_name,
    metal,
    df_meta,
    save_dir = file.path(output_dir, "dependence_plots"),
    n_boot   = 500L) {
  
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  
  df <- tibble(
    x     = as.numeric(x),
    shap  = as.numeric(shap),
    CaCO3 = factor(df_meta$CaCO3_cat, levels = c(0, 1)),
    Plant = factor(df_meta$plant_cat,  levels = c(0, 1))
  )
  if (nrow(df) < 4 || var(df$x, na.rm = TRUE) == 0)
    return(invisible(NULL))
  
  dens <- density(df$x, na.rm = TRUE, adjust = 1.2)
  dens_df <- tibble(
    x = dens$x,
    y = dens$y / max(dens$y) * max(table(cut(df$x, 30)))
  )
  
  hist_dat <- hist(df$x, breaks = 30, plot = FALSE)
  hist_df  <- tibble(
    x_mid = hist_dat$mids,
    y     = hist_dat$counts / max(hist_dat$counts) * max(dens_df$y)
  )
  
  low_main <- lowess(df$x, df$shap, f = .35, iter = 2)
  grid_x   <- sort(unique(low_main$x))
  
  set.seed(123)
  boot_mat <- replicate(
    n_boot,
    {
      idx   <- sample(seq_len(nrow(df)), replace = TRUE)
      low_b <- lowess(df$x[idx], df$shap[idx], f = .35, iter = 2)
      approx(low_b$x, low_b$y, xout = grid_x, rule = 2)$y
    }
  )
  
  low_df <- tibble(
    x  = grid_x,
    y  = approx(low_main$x, low_main$y, xout = grid_x)$y,
    lo = apply(boot_mat, 1, quantile, .025),
    hi = apply(boot_mat, 1, quantile, .975)
  )
  
  y_max     <- max(abs(c(df$shap, low_df$lo, low_df$hi))) * 1.05
  scale_fac <- max(dens_df$y) / y_max
  df      <- mutate(df,      shap_s = shap * scale_fac)
  low_df  <- mutate(low_df,  y_s  = y  * scale_fac,
                    lo_s = lo * scale_fac,
                    hi_s = hi * scale_fac)
  
  col_vals   <- c("0" = "grey70", "1" = "black")
  shape_vals <- c("0" = 17,       "1" = 16)
  
  thr_ix <- which(diff(sign(low_df$y)) != 0)
  thr_df <- tibble(x = low_df$x[thr_ix])
  
  p_dep <- ggplot() +
    geom_area(data = dens_df, aes(x = x, y = y),
              fill = "#CAB7FF", alpha = .35) +
    geom_col(data = hist_df, aes(x = x_mid, y = y),
             width = diff(hist_dat$breaks)[1] * .9,
             fill  = "#CAB7FF", alpha = .35) +
    geom_ribbon(data = low_df, aes(x = x, ymin = lo_s, ymax = hi_s),
                fill = "grey60", alpha = .20) +
    geom_line(data = low_df, aes(x = x, y = y_s),
              colour = "#8A2BE2", linewidth = 1.1) +
    geom_point(data = df, aes(x = x, y = shap_s, colour = CaCO3, shape = Plant),
               size = 2, alpha = 0.75) +
    scale_colour_manual(values = col_vals, guide = "none") +
    scale_shape_manual(values  = shape_vals, guide = "none") +
    geom_vline(data = thr_df, aes(xintercept = x),
               linetype = "dotted", colour = "forestgreen") +
    geom_text(data = thr_df, aes(x = x, y = 0, label = sprintf("%.2f", x)),
              vjust = -0.7, colour = "forestgreen", size = 3.2) +
    scale_y_continuous(
      name   = "Density-scaled counts",
      limits = c(-max(dens_df$y), max(dens_df$y)),
      sec.axis = sec_axis(~ . / scale_fac,
                          name   = "SHAP value (impact on output)",
                          breaks = pretty(c(-y_max, y_max)))
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = sprintf("%s – Dependence of %s", metal, feature_name),
         x = feature_name) +
    theme_bw(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", hjust = .5),
      axis.title.y.right = element_text(margin = margin(l = 6)),
      axis.text.y.right  = element_text(colour = "#00008B"),
      axis.text.y.left   = element_text(colour = "grey40")
    )
  
  ggsave(
    file.path(save_dir,
              sprintf("%s_%s_dependence.png",
                      metal, gsub('[^A-Za-z0-9_]', '_', feature_name))),
    p_dep, width = 8, height = 6, dpi = 600
  )
  
  print(p_dep)
  invisible(p_dep)
}

################################################################################
# 6. Main loop: model + SHAP + figures
################################################################################
results <- list()

for (sheet in sheet_names) {
  
  df <- clean_sheet(sheet)
  
  set.seed(2024)
  hold_idx <- createDataPartition(df[[target_col]], p = 0.75, list = FALSE)
  main_df  <- df[ hold_idx, ]
  hold_df  <- df[-hold_idx, ]
  
  dtrain <- xgb.DMatrix(make_matrix(main_df), label = main_df[[target_col]])
  dhold  <- xgb.DMatrix(make_matrix(hold_df),  label = hold_df[[target_col]])
  
  best_rmse <- Inf
  for (i in seq_len(nrow(param_grid))) {
    p      <- as.list(param_grid[i, ])
    params <- c(list(booster = "gbtree",
                     objective = "reg:squarederror",
                     eval_metric = "rmse",
                     tree_method = "hist"), p)
    cv <- xgb.cv(params, dtrain, nrounds = 2000, nfold = 5,
                 early_stopping_rounds = 20, verbose = 0)
    rmse_mean <- cv$evaluation_log[cv$best_iteration, "test_rmse_mean"]
    if (rmse_mean < best_rmse) {
      best_rmse   <- rmse_mean
      best_param  <- params
      best_nround <- cv$best_iteration
    }
  }
  
  best_model <- xgb.train(best_param, dtrain, nrounds = best_nround, verbose = 0)
  
  pred_train <- predict(best_model, dtrain)
  pred_hold  <- predict(best_model, dhold)
  
  rmse <- \(a, b) sqrt(mean((a - b)^2))
  mae  <- \(a, b) mean(abs(a - b))
  r2   <- \(a, b) cor(a, b)^2
  
  eval_tbl <- data.frame(
    Metric  = c("R2", "RMSE", "MAE"),
    Train   = c(r2(main_df[[target_col]], pred_train),
                rmse(main_df[[target_col]], pred_train),
                mae (main_df[[target_col]], pred_train)),
    Holdout = c(r2(hold_df[[target_col]], pred_hold),
                rmse(hold_df[[target_col]], pred_hold),
                mae (hold_df[[target_col]], pred_hold))
  )
  print(eval_tbl)
  
  # NOTE: Remaining plotting/export steps are unchanged (omitted here for brevity)
  # You can keep your full plotting section as-is.
  
  X_full <- make_matrix(df)
  sv     <- shapviz(best_model, X = X_full, X_pred = X_full)
  
  p_imp <- sv_importance(sv) +
    ggtitle(sprintf("%s – SHAP Importance", sheet)) +
    theme(plot.title = element_text(face = "bold", hjust = .5))
  
  ggsave(file.path(output_dir, sprintf("%s_SHAP_importance.png", sheet)),
         p_imp, width = 8, height = 6, dpi = 600)
  
  shap_contrib <- predict(best_model, newdata = X_full, predcontrib = TRUE)
  shap_contrib <- shap_contrib[, -ncol(shap_contrib), drop = FALSE]
  
  feat_target <- dep_feats[[sheet]]
  if (feat_target %in% colnames(shap_contrib)) {
    plot_shap_dependence(
      x            = df[[feat_target]],
      shap         = shap_contrib[, feat_target],
      feature_name = feat_target,
      metal        = sheet,
      df_meta      = df
    )
  }
  
  results[[sheet]] <- list(
    param  = best_param,
    nround = best_nround,
    eval   = eval_tbl
  )
}

################################################################################
# 7. Summary table
################################################################################
summary_dt <- rbindlist(lapply(names(results), \(s) {
  e <- results[[s]]$eval
  data.table(
    Metal        = s,
    R2_holdout   = e$Holdout[e$Metric == "R2"],
    RMSE_holdout = e$Holdout[e$Metric == "RMSE"],
    MAE_holdout  = e$Holdout[e$Metric == "MAE"],
    Best_nround  = results[[s]]$nround
  )
}))

print(summary_dt)
cat("\nAll models and figures have been exported to:", output_dir, "\n")