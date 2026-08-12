library(tidyverse)
library(dplyr)
library(ggplot2)
library(tidyr)
library(car)
library(patchwork)
library(gam)

df_final_merged <- read.csv("./data/faster_tdt_merged.csv")
df_final_merged <- df_final_merged %>%
  mutate(
    site = as.numeric(str_extract(id, "(?<=S)\\d+")),
    
    population_category = case_when(
      site == 13 ~ "tableland",
      site %in% c(14, 12, 7) ~ "montane",
      site %in% c(5, 6, 9, 10) ~ "sub alpine",
      TRUE ~ NA_character_
    ),
    
    population_category = factor(
      population_category,
      levels = c("tableland", "montane", "sub alpine")
    )
  )



###############################  Fig. 4: Heat map of bivariate correlations ###############################
vars <- c("Amax", "breadth_90","T_opt","gsw_mean_above30")



df_sub <- df_final_merged %>% 
  dplyr::select(all_of(vars))

# correlations
cor_mat <- cor(df_sub, use = "pairwise.complete.obs")

# p-values
rc <- Hmisc::rcorr(as.matrix(df_sub))
p_mat <- rc$P

# long format
cor_long <- as.data.frame(as.table(cor_mat))
p_long   <- as.data.frame(as.table(p_mat))

colnames(cor_long) <- c("Var1","Var2","value")
colnames(p_long)   <- c("Var1","Var2","pvalue")

# combine
cor_long <- dplyr::left_join(
  cor_long,
  p_long,
  by = c("Var1","Var2")
)

# keep lower triangle
cor_long <- cor_long %>%
  mutate(
    Var1 = factor(Var1, levels = vars),
    Var2 = factor(Var2, levels = vars),
    
    row_num = as.numeric(Var1),
    col_num = as.numeric(Var2),
    
    keep = row_num >= col_num,
    
    sig = round(pvalue, 3) <= 0.050,
    
    R2 = value^2,
    
    fill_group = case_when(
      !keep ~ "blank",
      !sig ~ "ns",
      sig & value > 0 ~ "pos",
      sig & value < 0 ~ "neg",
      TRUE ~ "blank"
    )
  ) %>%
  filter(keep)
cor_long

pretty_labels <- c(
  Amax = expression(italic(A)[max]),
  breadth_90 = expression(italic(T)[breadth]),
  T_opt = expression(italic(T)[opt]),
  E_D = expression(italic(E)[d*","~A]),
  T50_prime = expression(italic(CT)["max,1min"]),
  gsw_mean_above30 = expression(italic(g)[s])
)

p <- ggplot(cor_long, aes(x = Var2, y = Var1)) +
  
  geom_tile(
    aes(fill = fill_group),
    color = "white"
  ) +
  
  geom_text(
    aes(label = ifelse(Var1 == Var2, "", sprintf("%.2f", value))),
    size = 4,
    color = "black"
  )+
  
  scale_fill_manual(
    values = c(
      neg = "#00798c",
      pos = "#edae49",
      ns = "grey85",
      blank = "white"
    ),
    guide = "none"
  ) +
  
  scale_x_discrete(
    labels = pretty_labels,
    drop = FALSE
  ) +
  
  scale_y_discrete(
    labels = pretty_labels,
    drop = FALSE,
    limits = rev(vars)   # flip y-axis only
  ) +
  
  labs(x = NULL, y = NULL) +
  
  coord_fixed() +
  
  theme_minimal() +
  
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14)
  )

p
ggsave("./figures_updated/Fig4.pdf", plot = p, width = 6, height = 5)


###########################################  Fig. 5: gsw by Tleaf #########################################
gsw_by_tleaf <- read.csv("./data/at_clean_trimmed_with_site.csv")


#define temperature bins
curve_coverage <- gsw_by_tleaf %>%
  mutate(Tbin = round(Tleaf, 1)) %>% 
  group_by(Tbin) %>%
  summarise(
    n_curves = n_distinct(curveID),
    .groups = "drop"
  ) %>%
  mutate(
    prop_curves = n_curves / n_distinct(df_final_merged$curveID)
  )


majority_range <- curve_coverage %>%
  filter(prop_curves >= 0.7) %>%
  summarise(
    Tmin = min(Tbin),
    Tmax = max(Tbin)
  )
majority_range$Tmin
majority_range$Tmax

# temperature grid restricted to majority-supported range
T_grid <- seq(
  majority_range$Tmin,
  majority_range$Tmax,
  length.out = 200
)

library(mgcv)
curve_preds <- gsw_by_tleaf %>%
  group_by(curveID) %>%
  group_split() %>%
  map_dfr(function(dat) {
    
    fit <- mgcv::gam(
      gsw ~ s(Tleaf, k = 4),
      data = dat,
      method = "REML"
    )
    
    data.frame(
      curveID = unique(dat$curveID),
      Tleaf = T_grid,
      gsw_pred = predict(
        fit,
        newdata = data.frame(Tleaf = T_grid)
      )
    )
  })

# average across all curves
all_curve_average <- curve_preds %>%
  group_by(Tleaf) %>%
  summarise(
    gsw_mean = mean(gsw_pred, na.rm = TRUE),
    gsw_median = median(gsw_pred, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(gsw_by_tleaf, aes(x = Tleaf, y = gsw)) +
  
  # individual curves underneath
  geom_smooth(
    aes(group = curveID,
        colour = population_category),
    method = "gam",
    formula = y ~ s(x, k = 4),
    se = FALSE,
    linewidth = 0.25
  ) +
  
  # pooled population-average curve
  geom_line(
    data = all_curve_average,
    aes(
      x = Tleaf,
      y = gsw_mean
    ),
    colour = "black",
    linewidth = 1.5
  )+
  scale_colour_manual(
    values = c(
      "tableland"   = alpha("#800080", 0.4),
      "montane"     = alpha("#FFAE03", 0.4),
      "sub alpine"  = alpha("#8AAA79", 0.4)
    )
  ) +
  
  theme_classic() +
  theme(legend.position = "right") +
  labs(
    x = "Leaf temperature (°C)",
    y = expression(g[s])
  )

ggsave("./figures_updated/Fig5.pdf", plot = last_plot(), width = 15, height = 10, units = "cm")

# Linear mixed-effects model: gs ~ Tleaf within the majority-
# supported temperature window, with curve identity as a
# random intercept (matches Methods text).

library(lme4)
library(lmerTest)

gsw_mixed_data <- gsw_by_tleaf %>%
  filter(
    Tleaf >= majority_range$Tmin,
    Tleaf <= majority_range$Tmax
  )

mod_gsw_mixed <- lmer(
  gsw ~ Tleaf + (1 | curveID),
  data = gsw_mixed_data
)

summary(mod_gsw_mixed)          # fixed effect estimate + Satterthwaite p-value (lmerTest)
confint(mod_gsw_mixed, method = "Wald")   # optional: CI on the Tleaf slope

######################################  Fig. 6: Partial correlations  ######################################
library(car)

df_partial <- df_final_merged

# optional: nicer plotting layout
par(mfrow = c(2,2))

###### Tbreadth #######
mod_breadth <- lm(
  Amax ~ breadth_90 + gsw_mean_above30,
  data = df_partial
)

summary(mod_breadth)

###### Topt #####
mod_Topt <- lm(
  Amax ~ T_opt + gsw_mean_above30,
  data = df_partial
)

summary(mod_Topt)


var_labels <- c(
  T_opt = expression(italic(T)[opt]~"(°C)"),
  breadth_90 = "Thermal breadth (°C)"
)

plot_partial_gsw <- function(var, data) {
  
  df <- data %>%
    select(
      Amax,
      all_of(var),
      gsw_mean_above30,
      population_category
    ) %>%
    na.omit()
  
  # residualize Amax against gsw
  res_A <- resid(
    lm(Amax ~ gsw_mean_above30, data = df)
  )
  
  # residualize focal predictor against gsw
  res_X <- resid(
    lm(
      reformulate(
        "gsw_mean_above30",
        response = var
      ),
      data = df
    )
  )
  
  df_plot <- data.frame(
    res_X = res_X,
    res_A = res_A,
    population_category = df$population_category
  )
  
  # statistics
  summary_model <- summary(lm(res_A ~ res_X))
  
  R2   <- round(summary_model$r.squared, 3)
  pval <- signif(summary_model$coefficients[2,4], 3)
 # cor(res_X, res_Y)
  
  
  label_text <- sprintf(
    "atop(R^2 == %.3f, p == %.3g)",
    R2,
    pval
  )
  
  ggplot(df_plot, aes(x = res_X, y = res_A)) +
    
    geom_point(
      aes(color = population_category),
      alpha = 0.6,
      size = 2
    ) +
    
    geom_smooth(
      method = "lm",
      color = "black",
      fill = "black",
      alpha = 0.2,
      se = FALSE,
      linewidth = 2
    ) +
    
    scale_color_manual(values = c(
      "tableland" = "#800080",
      "montane" = "#FFAE03",
      "sub alpine" = "#8AAA79"
    )) +
    
    xlab(var_labels[[var]]) +
    
    ylab(expression(
      italic(A)[max]~
        "(" * mu*"mol"~m^{-2}~s^{-1}*")"
    )) +
    
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = label_text,
      hjust = 1.3,
      vjust = 1.2,
      size = 5,
      parse = TRUE
    ) +
    
    theme_classic() +
    
    theme(
      legend.position = "none",
      axis.title = element_text(size = 14),
      axis.text  = element_text(size = 12)
    )
}

########Check VIF
vif(mod_Topt)
vif(mod_breadth)

########## Make plots
p_Topt <- plot_partial_gsw("T_opt", df_partial)
p_breadth <- plot_partial_gsw("breadth_90", df_partial)

########## Display 
p_Topt
p_breadth

all<-(p_breadth | p_Topt)
all

ggsave(
  "./figures_updated/Figure5.pdf",
  all,
  height = 8,
  width = 10,
  units = "in"
)

########## Save 
ggsave(
  "./figures_updated/Amax_versus_Topt_partial_gsw.pdf",
  p_Topt,
  height = 12,
  width = 12,
  units = "cm"
)

ggsave(
  "./figures_updated/Amax_versus_breadth_partial_gsw.pdf",
  p_breadth,
  height = 12,
  width = 12,
  units = "cm"
)

########################### Fig. 7 Topt versus Ed,A & T50 ###########################
plot_partial_model <- function(response, focal_predictor,
                               predictors, data,
                               xlab_expr, ylab_expr,
                               xlim = NULL, ylim = NULL) {
  
  # Variables required for this model
  df <- data %>%
    select(
      all_of(c(response, predictors)),
      population_category
    ) %>%
    na.omit()
  
  # Full model
  full_mod <- lm(
    reformulate(predictors, response = response),
    data = df
  )
  
  full_s <- summary(full_mod)
  
  full_R2 <- full_s$r.squared
  
  # Control variables
  control_vars <- predictors[predictors != focal_predictor]
  
  # Partial residuals
  res_Y <- resid(
    lm(
      reformulate(control_vars, response = response),
      data = df
    )
  )
  
  res_X <- resid(
    lm(
      reformulate(control_vars, response = focal_predictor),
      data = df
    )
  )
  
  df_plot <- data.frame(
    res_X = res_X,
    res_Y = res_Y,
    population_category = df$population_category
  )
  
  # Partial model
  partial_mod <- lm(res_Y ~ res_X)
  partial_s <- summary(partial_mod)
  
  partial_r <- cor(res_X, res_Y)
  partial_R2 <- partial_s$r.squared
  pval <- coef(partial_s)[2, 4]
  
  label_text <- sprintf(
    "r = %.3f\npartial R² = %.3f\np = %.3g\nfull R² = %.3f",
    partial_r,
    partial_R2,
    pval,
    full_R2
  )
  
  # Plot
  p <- ggplot(df_plot, aes(res_X, res_Y)) +
    
    geom_point(
      aes(colour = population_category),
      alpha = 0.6,
      size = 2
    ) +
    
    geom_smooth(
      method = "lm",
      colour = "black",
      fill = "black",
      alpha = 0.2,
      linewidth = 2,
      se = FALSE
    ) +
    
    scale_color_manual(values = c(
      "tableland" = "#800080",
      "montane" = "#FFAE03",
      "sub alpine" = "#8AAA79"
    )) +
    
    xlab(xlab_expr) +
    ylab(ylab_expr) +
    
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = label_text,
      hjust = 1.3,
      vjust = 1.2,
      size = 4.5
    ) +
    
    theme_classic() +
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.position = "none"
    )
  
  # Apply limits WITHOUT removing observations
  if (!is.null(xlim) || !is.null(ylim)) {
    p <- p + coord_cartesian(
      xlim = xlim,
      ylim = ylim
    )
  }
  
  return(p)
}

get_partial_x <- function(focal_predictor, data) {
  
  predictors <- c("T_opt", "T50", "gsw_mean_above30")
  
  df <- data %>%
    select(
      all_of(predictors),
      population_category
    ) %>%
    na.omit()
  
  control_vars <- predictors[predictors != focal_predictor]
  
  resid(
    lm(
      reformulate(control_vars, response = focal_predictor),
      data = df
    )
  )
}


x_Topt <- get_partial_x("T_opt", df_final_merged)
x_T50  <- get_partial_x("T50", df_final_merged)

x_range <- range(
  c(x_Topt, x_T50),
  na.rm = TRUE
)

x_range
x_range <- x_range + c(-1, 1)

get_partial_y <- function(response, data) {
  
  predictors <- c("T_opt", "T50", "gsw_mean_above30")
  
  df <- data %>%
    select(
      all_of(c(response, predictors)),
      population_category
    ) %>%
    na.omit()
  
  control_vars <- predictors[predictors != response]
  
  resid(
    lm(
      reformulate(control_vars, response = response),
      data = df
    )
  )
}


y_Tbreadth <- get_partial_y("breadth_90", df_final_merged)
y_EdA      <- get_partial_y("E_D", df_final_merged)

y_Tbreadth_range <- range(y_Tbreadth, na.rm = TRUE)
y_EdA_range      <- range(y_EdA, na.rm = TRUE)

y_Tbreadth_range <- y_Tbreadth_range + c(-1, 1)
y_EdA_range      <- y_EdA_range + c(-0.1, 0.1)

p_Tbreadth <- list(
  
  plot_partial_model(
    response = "breadth_90",
    focal_predictor = "T_opt",
    predictors = c("T_opt", "T50", "gsw_mean_above30"),
    data = df_final_merged,
    xlab_expr = expression(italic(T)[opt] ~ "(°C)"),
    ylab_expr = expression(italic(T)[breadth] ~ "(°C)"),
    xlim = x_range,
    ylim = y_Tbreadth_range
  ),
  
  plot_partial_model(
    response = "breadth_90",
    focal_predictor = "T50",
    predictors = c("T_opt", "T50", "gsw_mean_above30"),
    data = df_final_merged,
    xlab_expr = expression(italic(T)[50] ~ "(°C)"),
    ylab_expr = expression(italic(T)[breadth] ~ "(°C)"),
    xlim = x_range,
    ylim = y_Tbreadth_range
  )
)


p_EdA <- list(
  
  plot_partial_model(
    response = "E_D",
    focal_predictor = "T_opt",
    predictors = c("T_opt", "T50", "gsw_mean_above30"),
    data = df_final_merged,
    xlab_expr = expression(italic(T)[opt] ~ "(°C)"),
    ylab_expr = expression(italic(E)[d*","*A] ~ "(eV)"),
    xlim = x_range,
    ylim = y_EdA_range
  ),
  
  plot_partial_model(
    response = "E_D",
    focal_predictor = "T50",
    predictors = c("T_opt", "T50", "gsw_mean_above30"),
    data = df_final_merged,
    xlab_expr = expression(italic(T)[50] ~ "(°C)"),
    ylab_expr = expression(italic(E)[d*","*A] ~ "(eV)"),
    xlim = x_range,
    ylim = y_EdA_range
  )
)
all <- (p_Tbreadth[[1]] | p_Tbreadth[[2]]) /
  (p_EdA[[1]]     | p_EdA[[2]])

all

########## Save 
ggsave(
  "./figures_updated/Tbreadth_Topt.pdf",
  p_Tbreadth[[1]],
  height = 12,
  width = 12,
  units = "cm"
)

ggsave(
  "./figures_updated/Tbreadth_T50.pdf",
  p_Tbreadth[[2]],
  height = 12,
  width = 12,
  units = "cm"
)
ggsave(
  "./figures_updated/Ed_Topt.pdf",
  p_EdA[[1]],
  height = 12,
  width = 12,
  units = "cm"
)

ggsave(
  "./figures_updated/Ed_T50.pdf",
  p_EdA[[2]],
  height = 12,
  width = 12,
  units = "cm"
)
############## Table 2 ##############
get_partial_stats <- function(var, data) {
  
  df <- data %>%
    select(
      Amax,
      all_of(var),
      gsw_mean_above30
    ) %>%
    na.omit()
  
  # residualize Amax against gsw
  res_A <- resid(
    lm(Amax ~ gsw_mean_above30, data = df)
  )
  
  # residualize predictor against gsw
  res_X <- resid(
    lm(
      reformulate(
        "gsw_mean_above30",
        response = var
      ),
      data = df
    )
  )
  
  # partial regression
  mod <- lm(res_A ~ res_X)
  mod_summary <- summary(mod)
  
  estimate <- coef(mod_summary)[2, "Estimate"]
  pvalue <- coef(mod_summary)[2, "Pr(>|t|)"]
  
  # partial R2 (matches plot annotation)
  partial_R2 <- mod_summary$r.squared
  
  data.frame(
    predictor = var,
    estimate = estimate,
    pvalue = pvalue,
    partial_R2 = partial_R2,
    n = nrow(df)
  )
}

partial_results <- bind_rows(
  get_partial_stats("E_D", df_partial),
  get_partial_stats("breadth_90", df_partial),
  get_partial_stats("T_opt", df_partial),
  get_partial_stats("T50_prime", df_partial),
  get_partial_stats("T50", df_partial)
)


partial_results

########## Check VIF — 3-predictor models (Fig. 7 / Table 2) ##########

# These reproduce the full_mod fit inside plot_partial_model() and
# get_partial_stats(), just so vif() has a model object to act on.

vif_breadth_full <- lm(
  breadth_90 ~ T_opt + T50 + gsw_mean_above30,
  data = df_final_merged
)
vif(vif_breadth_full)

vif_EdA_full <- lm(
  E_D ~ T_opt + T50 + gsw_mean_above30,
  data = df_final_merged
)
vif(vif_EdA_full)

# Topt-as-response models (used for the Topt_partial_results table)
vif_Topt_EdA <- lm(
  T_opt ~ E_D + gsw_mean_above30,
  data = df_final_merged
)
vif(vif_Topt_EdA)

vif_Topt_T50 <- lm(
  T_opt ~ T50 + gsw_mean_above30,
  data = df_final_merged
)
vif(vif_Topt_T50)

# Table 2's remaining 2-predictor Amax models (E_D, T50_prime, T50 —
# breadth_90 and T_opt are already covered by mod_breadth/mod_Topt above)
vif_Amax_EdA   <- lm(Amax ~ E_D       + gsw_mean_above30, data = df_final_merged)
vif_Amax_T50p  <- lm(Amax ~ T50_prime + gsw_mean_above30, data = df_final_merged)
vif_Amax_T50   <- lm(Amax ~ T50       + gsw_mean_above30, data = df_final_merged)

vif(vif_Amax_EdA)
vif(vif_Amax_T50p)
vif(vif_Amax_T50)

############################### Supporting Information ##############################
############## Table 1 ##############
site_sample_sizes <- df_final_merged %>%
  count(site, name = "n") %>%
  arrange(site)

site_sample_sizes

population_sample_sizes <- df_final_merged %>%
  count(population_category, name = "n") %>%
  arrange(population_category)

population_sample_sizes

vars <- c(
  "Amax", "breadth_90", "T_opt", "E_D", "T50",
  "Ea_kJmol", "z", "T50_prime", "lnA", "gsw_mean_above30"
)

sample_sizes <- df_final_merged %>%
  summarise(
    across(all_of(vars), ~ sum(!is.na(.)))
  ) %>%
  tidyr::pivot_longer(
    everything(),
    names_to = "Trait",
    values_to = "Sample_size"
  )

sample_sizes


sample_sizes_pop <- df_final_merged %>%
  group_by(population_category) %>%
  summarise(
    across(all_of(vars), ~ sum(!is.na(.))),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    -population_category,
    names_to = "Trait",
    values_to = "N"
  ) %>%
  tidyr::pivot_wider(
    names_from = population_category,
    values_from = N
  )

sample_sizes_pop



##############Table S2 ###########
df_sub <- df_final_merged %>%
  select(breadth_90, Amax, T_opt, E_D, T50,Ea_kJmol,lnA,T50,z,T50_prime, gsw_mean_above30, population_category) %>%
  drop_na()

trait_summary <- df_sub %>%
  select(-population_category) %>%
  
  summarise(
    across(
      everything(),
      list(
        min = ~round(min(.x, na.rm = TRUE), 2),
        max = ~round(max(.x, na.rm = TRUE), 2),
        mean = ~round(mean(.x, na.rm = TRUE), 2),
        median = ~round(median(.x, na.rm = TRUE), 2),
        SE = ~round(sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))), 2)
      )
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("trait", "statistic"),
    names_sep = "_(?=[^_]+$)",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = statistic,
    values_from = value
  )
trait_summary %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))
trait_order <- c(
  "Amax",
  "breadth_90",
  "T_opt",
  "E_D",
  "Ea_kJmol",
  "T50",
  "z",
  "T50_prime",
  "lnA",
  "gsw_mean_above30"
)

Table2 <- trait_summary %>%
  slice(match(trait_order, trait))
Table2

####Topt versus Ed,A and T50 ##############
get_partial_stats_Topt <- function(var, data) {
  
  df <- data %>%
    select(
      T_opt,
      all_of(var),
      gsw_mean_above30
    ) %>%
    na.omit()
  
  # residualize T_opt against gsw
  res_Topt <- resid(
    lm(T_opt ~ gsw_mean_above30, data = df)
  )
  
  # residualize predictor against gsw
  res_X <- resid(
    lm(
      reformulate(
        "gsw_mean_above30",
        response = var
      ),
      data = df
    )
  )
  
  # partial regression
  mod <- lm(res_Topt ~ res_X)
  mod_summary <- summary(mod)
  
  estimate <- coef(mod_summary)[2, "Estimate"]
  pvalue <- coef(mod_summary)[2, "Pr(>|t|)"]
  
  # partial R2 (matches plot annotation)
  partial_R2 <- mod_summary$r.squared
  
  data.frame(
    predictor = var,
    estimate = estimate,
    pvalue = pvalue,
    partial_R2 = partial_R2,
    n = nrow(df)
  )
}

Topt_partial_results <- bind_rows(
  get_partial_stats_Topt("E_D", df_final_merged),
  get_partial_stats_Topt("T50", df_final_merged)
)

Topt_partial_results

############# Fig S2 ###############

library(maps)
library(dplyr)

df_clean <- df_final_merged %>%
  mutate(population_category = case_when(
    site == 13 ~ "tableland",
    site %in% c(14, 12, 7) ~ "montane",
    site %in% c(5, 6, 9, 10) ~ "sub alpine",
    TRUE ~ NA_character_  # in case there are other sites
  ))

# Custom colors for site groups
site_colors <- c(
  # Purple group 5,6,7 (slightly lighter/darker)
  "5" = "#8AAA79",  # original
  "6" =   "#8AAA79",# lighter purple
  "7" = "#FFAE03",  # darker purple
  
  # Red group 9,10
  "9" = "#8AAA79",  # original
  "10" = "#8AAA79", # lighter red
  
  # Yellow group 12,13,14
  "12" = "#FFAE03",  # original
  "13" =  "#800080", # lighter yellow
  "14" =   "#FFAE03" # darker yellow
)


df_clean <- df_clean %>% mutate(site = factor(site))

# Australia map
aus_map <- map_data("world", region = "Australia")

# Coordinates for major cities
cities <- data.frame(
  name = c("Sydney", "Canberra", "Melbourne"),
  lat  = c(-33.8688, -35.2809, -37.8136),
  long = c(151.2093, 149.1300, 144.9631)
)

site_labels <- df_clean %>%
  group_by(site) %>%
  summarise(
    long = first(long),
    lat = first(lat),
    .groups = "drop"
  ) %>%
  mutate(
    x_nudge = case_when(
      site %in% c("5", "9", "12") ~ 0.08,
      TRUE ~ -0.08
    ),
    y_nudge = case_when(
      site %in% c("5", "6")  ~ 0.08,
      site %in% c("9", "10") ~ -0.08,
      TRUE ~ 0.08
    )
  )
library(ggrepel)


# Plot
map<-ggplot() +
  geom_polygon(
    data = aus_map,
    aes(x = long, y = lat, group = group),
    fill = "white", color = "black"
  ) +
  geom_point(
    data = df_clean %>% filter(site %in% names(site_colors)),
    aes(x = long, y = lat, color = site),
    size = 3, alpha = 0.8
  ) +
  geom_point(
    data = cities,
    aes(x = long, y = lat),
    color = "black", size = 2
  ) +
  geom_text(
    data = cities,
    aes(x = long, y = lat, label = name),
    nudge_y = 0.2, size = 3, fontface = "bold"
  ) +
  geom_text(
    data = site_labels,
    aes(
      x = long + x_nudge,
      y = lat + y_nudge,
      label = site
    ),
    size = 3,
    fontface = "bold"
  )+
  scale_color_manual(values = site_colors, name = "Site Group") +
  coord_fixed(
    ratio = 1.3,
    xlim = c(144, 154),  # zoom into NSW
    ylim = c(-39, -29.5)
  ) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.6
    ),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 9, colour = "black"),
    axis.ticks = element_line(colour = "black"),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )
map
#ggsave("figures_updated/FigS2_v2.pdf", map, height = 5, width = 7, units = "in")

summary_table <- df_final_merged %>%
  group_by(population_category,site) %>%
  summarise(
    n = n(),
    elev_mean=mean(elev_m, na.rm = T),
    elev_se=sd(elev_m, na.rm = TRUE) / sqrt(n),
    MAP_mean = mean(BIO12, na.rm = TRUE),
    MAP_se   = sd(BIO12, na.rm = TRUE) / sqrt(n),
    MTWQ_mean = mean(BIO10, na.rm = TRUE),
    MTWQ_se   = sd(BIO10, na.rm = TRUE) / sqrt(n)
  )

summary_table

population_colors<-c(
  "tableland"   = alpha("#800080", 1),
  "montane"     = alpha("#FFAE03", 1),
  "sub alpine"  = alpha("#8AAA79", 1))


plot_data <- summary_table %>%
  select(
    population_category,
    elev_mean, elev_se,
    MTWQ_mean, MTWQ_se,
    MAP_mean, MAP_se
  ) %>%
  pivot_longer(
    cols = -population_category,
    names_to = c("variable", "stat"),
    names_pattern = "(elev|MTWQ|MAP)_(mean|se)",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value
  ) %>%
  mutate(
    variable = factor(
      variable,
      levels = c("elev", "MTWQ", "MAP"),
      labels = c(
        "Elevation",
        "Mean temperature of\nwarmest quarter",
        "Annual precipitation"
      )
    )
  )
plot_data$population_category <- factor(plot_data$population_category, levels = c("tableland","montane","sub alpine"))

p_climate <- ggplot(
  plot_data,
  aes(
    x = population_category,
    y = mean,
    colour = population_category
  )
) +
  geom_errorbar(
    aes(
      ymin = mean - se,
      ymax = mean + se
    ),
    width = 0.12,
    linewidth = 0.7
  ) +
  geom_point(size = 3) +
  facet_wrap(
    ~ variable,
    nrow = 1,
    scales = "free_y"
  ) +
  scale_colour_manual(values = population_colors) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(
      size = 14
    ))

p_climate

#ggsave("figures_updated/FigS2_climate_by_site.pdf", p_climate, height = 4, width = 10, units = "in")


############# Fig S3 ###############
df_sub <- df_final_merged

get_lm_label <- function(df, yvar) {
  fit <- lm(reformulate("elev_m", yvar), data = df)
  p <- summary(fit)$coefficients[2, 4]
  
  p_lab <- ifelse(p < 0.001, "p < 0.001",
                  paste0("p = ", signif(p, 2)))
  
  return(p_lab)
}

add_p_label <- function(p, label) {
  p + annotate(
    "text",
    x = -Inf, y = Inf,
    label = label,
    hjust = -0.1, vjust = 1.1,
    size = 4
  )
}

#df_sub$population_category <- as.factor(df_sub$population_category, levels = c("tableland","montane","sub alpine"))
Amax_by_elev<-ggplot(df_sub,aes(x=elev_m,y=Amax,color=population_category))+
  geom_point()+
  scale_colour_manual(
    values = c(
      "tableland"   = alpha("#800080", 0.4),
      "montane"     = alpha("#FFAE03", 0.4),
      "sub alpine"  = alpha("#8AAA79", 0.4)
    )
  ) +
  scale_x_continuous(
    breaks = c(800, 1200, 1600),
    limits = c(800, 1600)
  )+
  ylab(expression(italic(A)[max]~
                    "(" * mu*"mol"~m^{-2}~s^{-1}*")"))+
  #xlab(expression(italic(E)[a]~"(kJ mol"^{-1}*")"))+
  xlab("Elevation (m)")+
  theme_classic()+
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none")
Amax_by_elev
summary(lm(Amax~elev_m,data=df_sub))

Tbreadth_by_elev<-ggplot(df_sub,aes(x=elev_m,y=breadth_90,color=population_category))+
  geom_point()+
  scale_colour_manual(
    values = c(
      "tableland"   = alpha("#800080", 0.4),
      "montane"     = alpha("#FFAE03", 0.4),
      "sub alpine"  = alpha("#8AAA79", 0.4)
    )
  ) +
  scale_x_continuous(
    breaks = c(800, 1200, 1600),
    limits = c(800, 1600)
  )+
  scale_y_continuous(
    breaks = c(6,10, 14,18),
    limits = c(5, 18)
  )+
  ylab(expression(italic(T)[breadth]~"(°C)"))+
  #xlab(expression(italic(E)[a]~"(kJ mol"^{-1}*")"))+
  xlab("Elevation (m)")+
  theme_classic()+
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none")
Tbreadth_by_elev

Topt_by_elev<-ggplot(df_sub,aes(x=elev_m,y=T_opt,color=population_category))+
  geom_point()+
  scale_colour_manual(
    values = c(
      "tableland"   = alpha("#800080", 0.4),
      "montane"     = alpha("#FFAE03", 0.4),
      "sub alpine"  = alpha("#8AAA79", 0.4)
    )
  ) +
  scale_x_continuous(
    breaks = c(800, 1200, 1600),
    limits = c(800, 1600)
  )+
  ylab(expression(italic(T)[opt]~"(°C)"))+
  #xlab(expression(italic(E)[a]~"(kJ mol"^{-1}*")"))+
  xlab("Elevation (m)")+
  theme_classic()+
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none")
Topt_by_elev

Ed_by_elev<-ggplot(df_sub,aes(x=elev_m,y=E_D,color=population_category))+
  geom_point()+
  scale_colour_manual(
      values = c(
        "tableland"   = alpha("#800080", 0.4),
        "montane"     = alpha("#FFAE03", 0.4),
        "sub alpine"  = alpha("#8AAA79", 0.4)
      )
  ) +
  scale_x_continuous(
    breaks = c(800, 1200, 1600),
    limits = c(800, 1600)
  )+
  ylab(expression(italic(E)[d]~"(eV)")) +
  #xlab(expression(italic(E)[a]~"(kJ mol"^{-1}*")"))+
  xlab("Elevation (m)")+
  theme_classic()+
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none")
Ed_by_elev

T50_by_elev<-ggplot(df_sub,aes(x=elev_m,y=T50,color=population_category))+
  geom_point()+
  scale_colour_manual(
    values = c(
      "tableland"   = alpha("#800080", 0.4),
      "montane"     = alpha("#FFAE03", 0.4),
      "sub alpine"  = alpha("#8AAA79", 0.4)
    )
  ) +
  scale_x_continuous(
    breaks = c(800, 1200, 1600),
    limits = c(800, 1600)
  )+
  ylab(expression(italic(T)[50]~"(°C)")) +
  #xlab(expression(italic(E)[a]~"(kJ mol"^{-1}*")"))+
  xlab("Elevation (m)")+
  theme_classic()+
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none")
T50_by_elev

Ea_by_elev<-ggplot(df_sub,aes(x=elev_m,y=Ea_kJmol,color=population_category))+
  geom_point()+
  scale_colour_manual(
    values = c(
      "tableland"   = alpha("#800080", 0.4),
      "montane"     = alpha("#FFAE03", 0.4),
      "sub alpine"  = alpha("#8AAA79", 0.4)
    )
  ) +
  scale_x_continuous(
    breaks = c(800, 1200, 1600),
    limits = c(800, 1600)
  )+
  ylab(expression(italic(E)[a]~"(kJ mol"^{-1}*")"))+
  #xlab(expression(italic(E)[a]~"(kJ mol"^{-1}*")"))+
  xlab("Elevation (m)")+
  theme_classic()+
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none")
Ea_by_elev

z_by_elev<-ggplot(df_sub,aes(x=elev_m,y=z,color=population_category))+
  geom_point()+
  scale_colour_manual(
    values = c(
      "tableland"   = alpha("#800080", 0.4),
      "montane"     = alpha("#FFAE03", 0.4),
      "sub alpine"  = alpha("#8AAA79", 0.4)
    )
  ) +
  scale_x_continuous(
    breaks = c(800, 1200, 1600),
    limits = c(800, 1600)
  )+
  ylab(expression(italic(z)~"dimensionless"))+
  #xlab(expression(italic(E)[a]~"(kJ mol"^{-1}*")"))+
  xlab("Elevation (m)")+
  theme_classic()+
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none")
z_by_elev

T50_prime_by_elev <- ggplot(
  df_sub,
  aes(x = elev_m, y = T50_prime, color = population_category)
) +
  geom_point() +
  scale_colour_manual(
    values = c(
      "tableland"  = alpha("#800080", 0.4),
      "montane"    = alpha("#FFAE03", 0.4),
      "sub alpine" = alpha("#8AAA79", 0.4)
    )
  ) +
  scale_x_continuous(
    breaks = c(800, 1200, 1600),
    limits = c(800, 1600)
  ) +
  ylab(expression(italic(CT)[paste("max, 1 min")] ~ "(°C)")) +
  xlab("Elevation (m)") +
  theme_classic() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

T50_prime_by_elev

gsw_by_elev<-ggplot(df_sub,aes(x=elev_m,y=gsw_mean_above30,color=population_category))+
  geom_point()+
  scale_colour_manual(
    values = c(
      "tableland"  = alpha("#800080", 0.4),
      "montane"    = alpha("#FFAE03", 0.4),
      "sub alpine" = alpha("#8AAA79", 0.4)
    )
  ) +
  scale_x_continuous(
    breaks = c(800, 1200, 1600),
    limits = c(800, 1600)
  )+
  #ylab(expression(g[s]))+
  ylab(expression(italic(g)[s]~
                    "(mol "~m^{-2}~s^{-1}*")")) +
  xlab("Elevation (m)")+
  theme_classic()+
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none")
gsw_by_elev

library(patchwork)
Amax_by_elev <- add_p_label(Amax_by_elev, get_lm_label(df_sub, "Amax"))
Tbreadth_by_elev <- add_p_label(Tbreadth_by_elev, get_lm_label(df_sub, "breadth_90"))
Topt_by_elev     <- add_p_label(Topt_by_elev, get_lm_label(df_sub, "T_opt"))
Ed_by_elev       <- add_p_label(Ed_by_elev, get_lm_label(df_sub, "E_D"))
T50_by_elev      <- add_p_label(T50_by_elev, get_lm_label(df_sub, "T50"))
Ea_by_elev       <- add_p_label(Ea_by_elev, get_lm_label(df_sub, "Ea_kJmol"))
z_by_elev        <- add_p_label(z_by_elev, get_lm_label(df_sub, "z"))
T50_prime_by_elev<- add_p_label(T50_prime_by_elev, get_lm_label(df_sub, "T50_prime"))
gsw_by_elev      <- add_p_label(gsw_by_elev, get_lm_label(df_sub, "gsw_mean_above30"))

combined_plot <- (Amax_by_elev | Tbreadth_by_elev |Topt_by_elev) /
  (Ed_by_elev | T50_by_elev |Ea_by_elev) /
  (z_by_elev | T50_prime_by_elev|gsw_by_elev)
combined_plot

# Save as PDF
ggsave("./figures_updated/FigS3.pdf", combined_plot,
       width = 8, height = 8, units = "in")

#### Fig S4 ###############
df_final_merged
vars <- c("Amax", "breadth_90", "T_opt", "E_D", "T50",
          "Ea_kJmol", "z", "T50_prime","lnA", "gsw_mean_above30")
df_sub <- df_final_merged %>%
  dplyr::select(all_of(vars)) 

# correlations
cor_mat <- cor(df_sub, use = "pairwise.complete.obs")

# p-values
rc <- Hmisc::rcorr(as.matrix(df_sub))
p_mat <- rc$P

# long format
cor_long <- as.data.frame(as.table(cor_mat))
p_long   <- as.data.frame(as.table(p_mat))

colnames(cor_long) <- c("Var1", "Var2", "value")
colnames(p_long)   <- c("Var1", "Var2", "pvalue")

# combine
cor_long <- left_join(
  cor_long,
  p_long,
  by = c("Var1", "Var2")
)


# formatting
cor_long <- cor_long %>%
  mutate(
    Var1 = factor(Var1, levels = vars),
    Var2 = factor(Var2, levels = vars),
    
    row_num = match(as.character(Var1), vars),
    col_num = match(as.character(Var2), vars),
    
    # keep lower triangle only
    keep = row_num > col_num,
    
    sig = pvalue < 0.05,
    
    # add R2
    R2 = value^2,
    
    fill_group = case_when(
      !keep ~ "blank",
      !sig ~ "ns",
      value > 0 ~ "pos",
      value < 0 ~ "neg",
      TRUE ~ "blank"
    )
  ) %>%
  filter(keep)
#write.csv(cor_long, "figures_updated/cor_long_all.csv", row.names = FALSE)
cor_long

pretty_labels <- c(
  Amax = expression(italic(A)[max]),
  breadth_90 = expression(italic(T)[breadth]),
  T_opt = expression(italic(T)[opt]),
  E_D = expression(italic(E)[d*","~A]),
  T50 = expression(italic(T)[50]),
  Ea_kJmol = expression(italic(E)[a*","~PSII]),
  z = expression(italic(z)),
  T50_prime = expression(italic(CT)["max,1min"]),
  gsw_mean_above30 = expression(italic(g)[s])
)


p <- ggplot(cor_long, aes(x = Var2, y = Var1)) +
  
  geom_tile(
    aes(fill = fill_group),
    color = "white"
  ) +
  
  geom_text(
    aes(label = sprintf("%.2f", value)),
    size = 4,
    color = "black"
  ) +
  
  scale_fill_manual(
    values = c(
      pos = "#edae49",
      neg = "#00798c",
      ns = "grey85",
      blank = "white"
    ),
    guide = "none"
  ) +
  
  scale_x_discrete(
    labels = pretty_labels,
    drop = FALSE
  ) +
  
  scale_y_discrete(
    limits = rev(vars),
    labels = pretty_labels,
    drop = FALSE
  ) +
  
  labs(x = NULL, y = NULL) +
  
  coord_fixed() +
  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14)
  )

p
write.csv(cor_long,"pairwise_all_vars.csv")
ggsave(
  "./figures_updated/FigS4.pdf",
  plot = p,
  width = 8,
  height = 7
)

#Correlation between Topt and Ed,A within curves
cor.test(
  df_final_merged$T_opt,
  df_final_merged$E_D,
  use = "complete.obs",
  method = "pearson"
)


####################Fig S5  ###########
library(dplyr)
library(ggplot2)

# Exposure durations
exposure_times <- 1:100000

# Calculate correlation at each exposure duration
cor_results <- lapply(exposure_times, function(t) {
  
  CTmax <- df_final_merged$T50_prime -
    df_final_merged$z * log10(t)
  
  test <- cor.test(
    df_final_merged$z,
    CTmax,
    method = "pearson"
  )
  
  data.frame(
    time = t,
    r = unname(test$estimate),
    p = test$p.value
  )
}) %>%
  bind_rows() %>%
  mutate(
    significance = case_when(
      p < 0.05 & r > 0 ~ "Positive significant",
      p < 0.05 & r < 0 ~ "Negative significant",
      TRUE ~ "Non-significant"
    )
  )

ggplot(cor_results, aes(x = time, y = r)) +
  
  # Measured exposure-duration window
  annotate(
    "rect",
    xmin = 5,
    xmax = 120,
    ymin = -Inf,
    ymax = Inf,
    fill = "grey70",
    alpha = 0.25
  ) +
  
  # r = 0
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey40"
  ) +
  
  # Correlation curve
  #  geom_line(
  #    colour = "black",
  #    linewidth = 0.7
  #  ) +
  
  # Points coloured by significance/direction
  geom_point(
    aes(colour = significance),
    size = 2
  ) +
  
  scale_colour_manual(
    values = c(
      "Positive significant" = "red",
      "Non-significant" = "grey60",
      "Negative significant" = "blue"
    )
  ) +
  
  # Log10-scaled x axis
  scale_x_log10(
    breaks = c(1, 5, 10, 30, 60, 120, 1000,10000,100000)
  ) +
  
  scale_y_continuous(
    limits = c(-1, 1)
  ) +
  
  theme_classic() +
  
  labs(
    x = "Exposure duration (min)",
    y = expression(italic(r) ~ "(z vs. CT"[max] * ")"),
    colour = NULL
  ) +
  
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  )
ggsave("./figures_updated/FigS5.pdf", last_plot(),
       width = 8, height = 6, units = "in")

#Crossing temp
cor_results %>%
  mutate(
    sign = sign(r),
    previous_sign = lag(sign)
  ) %>%
  filter(sign != previous_sign) %>%
  select(time, r, previous_sign, sign)

################# Fig S6 ##################
library(broom)
# restrict to majority-supported temperature range
df_sub <- gsw_by_tleaf %>%
  filter(
    Tleaf >= majority_range$Tmin,
    Tleaf <= majority_range$Tmax
  )

# calculate slope within this temperature range for each curve
slopes <- df_sub %>%
  group_by(curveID, Population, population_category) %>%
  do({
    fit <- lm(gsw ~ Tleaf, data = .)
    tidy(fit)
  }) %>%
  filter(term == "Tleaf") %>%
  transmute(
    curveID,
    Population,
    population_category,
    slope = estimate,
    slope_se = std.error,
    pvalue = p.value
  )

# Create site-level climate table
climate_site <- df_final_merged %>%
  mutate(
    site = paste0("S", sprintf("%02d", site))
  ) %>%
  group_by(site) %>%
  summarise(
    elev_m = first(elev_m),
    BIO10  = first(BIO10),
    BIO12  = first(BIO12),
    .groups = "drop"
  )

# Add site identifier to slopes
slopes_climate <- slopes %>%
  mutate(site = Population) %>%
  left_join(climate_site, by = "site")
summary(lm(slope~elev_m,slopes_climate))


gsw_slope_elev <- ggplot(slopes_climate,aes(x=elev_m,y=slope,color=population_category))+
  geom_point()+
  scale_colour_manual(
    values = c(
      "tableland"   = alpha("#800080", 0.4),
      "montane"     = alpha("#FFAE03", 0.4),
      "sub alpine"  = alpha("#8AAA79", 0.4)
    )
  ) +
  ylab("Slope of stomatal conductance\nby leaf temperature")+
  #xlab(expression(italic(E)[a]~"(kJ mol"^{-1}*")"))+
  xlab("Elevation (m)")+
  theme_classic()+
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none")
gsw_slope_elev
summary(lm(slope~elev_m,slopes_climate))


gsw_slope_mtwq <-ggplot(slopes_climate,aes(x=BIO10,y=slope,color=population_category))+
  geom_point()+
  scale_colour_manual(
    values = c(
      "tableland"   = alpha("#800080", 0.4),
      "montane"     = alpha("#FFAE03", 0.4),
      "sub alpine"  = alpha("#8AAA79", 0.4)
    )
  ) +
  ylab("Slope of stomatal conductance\nby leaf temperature")+
  #xlab(expression(italic(E)[a]~"(kJ mol"^{-1}*")"))+
  xlab("Mean temperature warmest quarter (C)")+
  theme_classic()+
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none")
gsw_slope_mtwq
gsw_slope_precip <-ggplot(slopes_climate,aes(x=BIO12,y=slope,color=population_category))+
  geom_point()+
  scale_colour_manual(
    values = c(
      "tableland"   = alpha("#800080", 0.4),
      "montane"     = alpha("#FFAE03", 0.4),
      "sub alpine"  = alpha("#8AAA79", 0.4)
    )
  ) +
  ylab("Slope of stomatal conductance\nby leaf temperature")+
  #xlab(expression(italic(E)[a]~"(kJ mol"^{-1}*")"))+
  xlab("Mean annual precipitation (mm)")+
  theme_classic()+
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none")
gsw_slope_precip
library(patchwork)

combined_plot <- (gsw_slope_elev | gsw_slope_mtwq) /
  (gsw_slope_precip | plot_spacer())
combined_plot

# Save as PDF
ggsave("./figures_updated/FigS6.pdf", combined_plot,
       width = 8, height = 7, units = "in")


##################Fig S7 PCA ##################
df_pca <- df_final_merged %>%
  select(breadth_90, Amax, T_opt, E_D, T50,Ea_kJmol,lnA,T50,z,T50_prime, gsw_mean_above30, population_category) %>%
  drop_na()

# Run PCA (scale = TRUE is important!)# Rungsw_mean_above30 PCA (scale = TRUE is important!)
pca_res <- prcomp(
  df_pca %>% select(-population_category),
  scale = TRUE
)

# View PCA summary
summary(pca_res)

# Extract scores for plotting
scores <- as.data.frame(pca_res$x)
scores$population_category <- df_pca$population_category

# Variance explained for axis labels
var_exp <- round((pca_res$sdev^2) / sum(pca_res$sdev^2) * 100, 1)

library(ggplot2)

# Variable loadings (arrows)
loadings <- as.data.frame(pca_res$rotation)
loadings$var <- rownames(loadings)

scores$PC1 <- -scores$PC1
loadings$PC1 <- -loadings$PC1

#scores$PC2 <- -scores$PC2
#loadings$PC2 <- -loadings$PC2

# Scale arrows to visible length
arrow_multiplier <- 3

pca_plot <- ggplot(scores, aes(PC1, PC2, colour = population_category)) +
  geom_point(size = 3, alpha = 0.8) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  
  geom_segment(data = loadings,
               aes(x = 0, y = 0,
                   xend = PC1 * arrow_multiplier,
                   yend = PC2 * arrow_multiplier),
               arrow = arrow(length = unit(0.25, "cm")),
               linewidth = 0.8,
               inherit.aes = FALSE) +
  
  geom_text(data = loadings,
            aes(x = PC1 * arrow_multiplier * 1.1,
                y = PC2 * arrow_multiplier * 1.1,
                label = var),
            size = 5,
            inherit.aes = FALSE) +
  
  scale_colour_manual(values = c(
    "tableland" = "#800080",
    "montane" = "#FFAE03",
    "sub alpine" = "#8AAA79"
  )) +
  
  labs(
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)")
  ) +
  
  theme_classic(base_size = 14)+
  theme(legend.position = "none")
pca_plot


pca_res$rotation #Table 3
summary(pca_res)

ggsave("./figures_updated/FigS6_v2.pdf", pca_plot,
       width = 8, height = 8, units = "in")

