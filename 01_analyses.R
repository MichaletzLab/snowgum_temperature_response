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
vars <- c("Amax", "Tbreadth","Topt","gsw_mean_above30")



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
  Tbreadth = expression(italic(T)[breadth]),
  Topt = expression(italic(T)[opt]),
  Ed_A = expression(italic(E)[d*","~A]),
  CTmax_1min = expression(italic(CT)["max,1min"]),
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
  Amax ~ Tbreadth + gsw_mean_above30,
  data = df_partial
)

summary(mod_breadth)

###### Topt #####
mod_Topt <- lm(
  Amax ~ Topt + gsw_mean_above30,
  data = df_partial
)

summary(mod_Topt)


var_labels <- c(
  Topt = expression(italic(T)[opt]~"(°C)"),
  Tbreadth = "Thermal breadth (°C)"
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
p_Topt <- plot_partial_gsw("Topt", df_partial)
p_breadth <- plot_partial_gsw("Tbreadth", df_partial)

########## Display 
p_Topt
p_breadth

fig6 <- (p_breadth | p_Topt)
fig6

ggsave(
  "./figures_updated/Fig6.pdf",
  fig6,
  height = 8,
  width = 10,
  units = "in"
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
  
  predictors <- c("Topt", "T50", "gsw_mean_above30")
  
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


x_Topt <- get_partial_x("Topt", df_final_merged)
x_T50  <- get_partial_x("T50", df_final_merged)

x_range <- range(
  c(x_Topt, x_T50),
  na.rm = TRUE
)

x_range
x_range <- x_range + c(-1, 1)

get_partial_y <- function(response, data) {
  
  predictors <- c("Topt", "T50", "gsw_mean_above30")
  
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


y_Tbreadth <- get_partial_y("Tbreadth", df_final_merged)
y_EdA      <- get_partial_y("Ed_A", df_final_merged)

y_Tbreadth_range <- range(y_Tbreadth, na.rm = TRUE)
y_EdA_range      <- range(y_EdA, na.rm = TRUE)

y_Tbreadth_range <- y_Tbreadth_range + c(-1, 1)
y_EdA_range      <- y_EdA_range + c(-0.1, 0.1)

p_Tbreadth <- list(
  
  plot_partial_model(
    response = "Tbreadth",
    focal_predictor = "Topt",
    predictors = c("Topt", "T50", "gsw_mean_above30"),
    data = df_final_merged,
    xlab_expr = expression(italic(T)[opt] ~ "(°C)"),
    ylab_expr = expression(italic(T)[breadth] ~ "(°C)"),
    xlim = x_range,
    ylim = y_Tbreadth_range
  ),
  
  plot_partial_model(
    response = "Tbreadth",
    focal_predictor = "T50",
    predictors = c("Topt", "T50", "gsw_mean_above30"),
    data = df_final_merged,
    xlab_expr = expression(italic(T)[50] ~ "(°C)"),
    ylab_expr = expression(italic(T)[breadth] ~ "(°C)"),
    xlim = x_range,
    ylim = y_Tbreadth_range
  )
)


p_EdA <- list(
  
  plot_partial_model(
    response = "Ed_A",
    focal_predictor = "Topt",
    predictors = c("Topt", "T50", "gsw_mean_above30"),
    data = df_final_merged,
    xlab_expr = expression(italic(T)[opt] ~ "(°C)"),
    ylab_expr = expression(italic(E)[d*","*A] ~ "(eV)"),
    xlim = x_range,
    ylim = y_EdA_range
  ),
  
  plot_partial_model(
    response = "Ed_A",
    focal_predictor = "T50",
    predictors = c("Topt", "T50", "gsw_mean_above30"),
    data = df_final_merged,
    xlab_expr = expression(italic(T)[50] ~ "(°C)"),
    ylab_expr = expression(italic(E)[d*","*A] ~ "(eV)"),
    xlim = x_range,
    ylim = y_EdA_range
  )
)
fig7 <- (p_Tbreadth[[1]] | p_Tbreadth[[2]]) /
  (p_EdA[[1]]     | p_EdA[[2]])

fig7

ggsave(
  "./figures_updated/Fig7.pdf",
  fig7,
  height = 10,
  width = 10,
  units = "in"
)

############## Partial regression stats reported in Fig. 6 / Fig. 7 captions ##############
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
  get_partial_stats("Ed_A", df_partial),
  get_partial_stats("Tbreadth", df_partial),
  get_partial_stats("Topt", df_partial),
  get_partial_stats("CTmax_1min", df_partial),
  get_partial_stats("T50", df_partial)
)


partial_results

########## Check VIF — 3-predictor models (Fig. 7 / Table 2) ##########

# These reproduce the full_mod fit inside plot_partial_model() and
# get_partial_stats(), just so vif() has a model object to act on.

vif_breadth_full <- lm(
  Tbreadth ~ Topt + T50 + gsw_mean_above30,
  data = df_final_merged
)
vif(vif_breadth_full)

vif_EdA_full <- lm(
  Ed_A ~ Topt + T50 + gsw_mean_above30,
  data = df_final_merged
)
vif(vif_EdA_full)

# Topt-as-response models (used for the Topt_partial_results table)
vif_Topt_EdA <- lm(
  Topt ~ Ed_A + gsw_mean_above30,
  data = df_final_merged
)
vif(vif_Topt_EdA)

vif_Topt_T50 <- lm(
  Topt ~ T50 + gsw_mean_above30,
  data = df_final_merged
)
vif(vif_Topt_T50)

# Table 2's remaining 2-predictor Amax models (Ed_A, CTmax_1min, T50 —
# Tbreadth and Topt are already covered by mod_breadth/mod_Topt above)
vif_Amax_EdA   <- lm(Amax ~ Ed_A       + gsw_mean_above30, data = df_final_merged)
vif_Amax_T50p  <- lm(Amax ~ CTmax_1min + gsw_mean_above30, data = df_final_merged)
vif_Amax_T50   <- lm(Amax ~ T50       + gsw_mean_above30, data = df_final_merged)

vif(vif_Amax_EdA)
vif(vif_Amax_T50p)
vif(vif_Amax_T50)

############################### Supporting Information ##############################
############## Table S1: Site and population sampling description ##############
# Taxon labels: S10T18P03 (site 10) is recorded as "niphophila" in the
# source metadata, while every other site-10 individual is "hedraia" --
# per review, grouped with hedraia here (site 10 total n = 10) rather
# than broken out as its own row.
taxon_labels <- c(
  "pauciflora" = "Eucalyptus pauciflora subsp. pauciflora",
  "acerina"    = "Eucalyptus pauciflora subsp. acerina",
  "hedraia"    = "Eucalyptus pauciflora subsp. hedraia",
  "niphophila" = "Eucalyptus pauciflora subsp. hedraia",
  "lacrimans"  = "Eucalyptus lacrimans"
)

elev_class_labels <- c(
  "tableland"  = "Low (tableland)",
  "montane"    = "Mid (montane)",
  "sub alpine" = "High (subalpine)"
)

site_elev <- df_final_merged %>%
  group_by(site) %>%
  summarise(elev_m = round(mean(elev_m, na.rm = TRUE)), .groups = "drop")

table_s1 <- df_final_merged %>%
  mutate(Taxon = taxon_labels[species]) %>%
  count(population_category, Taxon, site, name = "n") %>%
  group_by(population_category) %>%
  mutate(N = sum(n)) %>%
  ungroup() %>%
  left_join(site_elev, by = "site") %>%
  mutate(`Source elevation class` = elev_class_labels[as.character(population_category)]) %>%
  arrange(population_category, site) %>%
  select(
    `Source elevation class`, N, Taxon, n,
    `Site identifier` = site, `Elevation (m)` = elev_m
  )

table_s1
write.csv(table_s1, "./figures_updated/TableS1.csv", row.names = FALSE)

############## Table S2: Sample sizes by parameter and source population ##############
vars <- c(
  "Amax", "Tbreadth", "Topt", "Ed_A", "Ea_PSII",
  "T50", "z", "CTmax_1min", "lnA_PSII", "gsw_mean_above30"
)

n_by_class <- df_final_merged %>%
  group_by(population_category) %>%
  summarise(across(all_of(vars), ~ sum(!is.na(.))), .groups = "drop") %>%
  tidyr::pivot_longer(-population_category, names_to = "trait", values_to = "n") %>%
  tidyr::pivot_wider(names_from = population_category, values_from = n)

trait_summary <- df_final_merged %>%
  summarise(
    across(
      all_of(vars),
      list(
        N      = ~ sum(!is.na(.x)),
        Min    = ~ round(min(.x, na.rm = TRUE), 2),
        Max    = ~ round(max(.x, na.rm = TRUE), 2),
        Mean   = ~ round(mean(.x, na.rm = TRUE), 2),
        Median = ~ round(median(.x, na.rm = TRUE), 2),
        SE     = ~ round(sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))), 2)
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

table_s2 <- trait_summary %>%
  left_join(n_by_class, by = "trait") %>%
  slice(match(vars, trait)) %>%
  select(
    Trait = trait, N,
    Low = tableland, Mid = montane, High = `sub alpine`,
    Min, Max, Mean, Median, SE
  )

table_s2
write.csv(table_s2, "./figures_updated/TableS2.csv", row.names = FALSE)

####Topt versus Ed,A and T50 ##############
get_partial_stats_Topt <- function(var, data) {
  
  df <- data %>%
    select(
      Topt,
      all_of(var),
      gsw_mean_above30
    ) %>%
    na.omit()
  
  # residualize Topt against gsw
  res_Topt <- resid(
    lm(Topt ~ gsw_mean_above30, data = df)
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
  get_partial_stats_Topt("Ed_A", df_final_merged),
  get_partial_stats_Topt("T50", df_final_merged)
)

Topt_partial_results

############# Fig S1 ###############

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

# One row per elevation class (not per site) -- the panels show
# class-level means/SEs across all individuals, per the Fig. S1 caption.
summary_table <- df_final_merged %>%
  group_by(population_category) %>%
  summarise(
    n = n(),
    elev_mean = mean(elev_m, na.rm = TRUE),
    elev_se   = sd(elev_m, na.rm = TRUE) / sqrt(n),
    MAP_mean  = mean(BIO12, na.rm = TRUE),
    MAP_se    = sd(BIO12, na.rm = TRUE) / sqrt(n),
    MTWQ_mean = mean(BIO10, na.rm = TRUE),
    MTWQ_se   = sd(BIO10, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
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

fig_s1 <- map | p_climate
fig_s1

ggsave("figures_updated/FigS1.pdf", fig_s1, height = 5, width = 12, units = "in")


############# Fig S2 ###############
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

Tbreadth_by_elev<-ggplot(df_sub,aes(x=elev_m,y=Tbreadth,color=population_category))+
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

Topt_by_elev<-ggplot(df_sub,aes(x=elev_m,y=Topt,color=population_category))+
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

Ed_by_elev<-ggplot(df_sub,aes(x=elev_m,y=Ed_A,color=population_category))+
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

Ea_by_elev<-ggplot(df_sub,aes(x=elev_m,y=Ea_PSII,color=population_category))+
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

CTmax_1min_by_elev <- ggplot(
  df_sub,
  aes(x = elev_m, y = CTmax_1min, color = population_category)
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

CTmax_1min_by_elev

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
Tbreadth_by_elev <- add_p_label(Tbreadth_by_elev, get_lm_label(df_sub, "Tbreadth"))
Topt_by_elev     <- add_p_label(Topt_by_elev, get_lm_label(df_sub, "Topt"))
Ed_by_elev       <- add_p_label(Ed_by_elev, get_lm_label(df_sub, "Ed_A"))
T50_by_elev      <- add_p_label(T50_by_elev, get_lm_label(df_sub, "T50"))
Ea_by_elev       <- add_p_label(Ea_by_elev, get_lm_label(df_sub, "Ea_PSII"))
z_by_elev        <- add_p_label(z_by_elev, get_lm_label(df_sub, "z"))
CTmax_1min_by_elev<- add_p_label(CTmax_1min_by_elev, get_lm_label(df_sub, "CTmax_1min"))
gsw_by_elev      <- add_p_label(gsw_by_elev, get_lm_label(df_sub, "gsw_mean_above30"))

combined_plot <- (Amax_by_elev | Tbreadth_by_elev |Topt_by_elev) /
  (Ed_by_elev | T50_by_elev |Ea_by_elev) /
  (z_by_elev | CTmax_1min_by_elev|gsw_by_elev)
combined_plot

# Save as PDF
ggsave("./figures_updated/FigS2.pdf", combined_plot,
       width = 8, height = 8, units = "in")

#### Fig S3 ###############
df_final_merged
vars <- c("Amax", "Tbreadth", "Topt", "Ed_A", "T50",
          "Ea_PSII", "z", "CTmax_1min","lnA_PSII", "gsw_mean_above30")
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
cor_long

pretty_labels <- c(
  Amax = expression(italic(A)[max]),
  Tbreadth = expression(italic(T)[breadth]),
  Topt = expression(italic(T)[opt]),
  Ed_A = expression(italic(E)[d*","~A]),
  T50 = expression(italic(T)[50]),
  Ea_PSII = expression(italic(E)[a*","~PSII]),
  z = expression(italic(z)),
  CTmax_1min = expression(italic(CT)["max,1min"]),
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
ggsave(
  "./figures_updated/FigS3.pdf",
  plot = p,
  width = 8,
  height = 7
)

#Correlation between Topt and Ed,A within curves
cor.test(
  df_final_merged$Topt,
  df_final_merged$Ed_A,
  use = "complete.obs",
  method = "pearson"
)

# Measurement-error correction: Topt and Ed_A are fit simultaneously from
# the same curve, so their observed correlation could partly reflect the
# fitting procedure rather than real among-individual variation. Check
# this by comparing each parameter's mean SE (fitting uncertainty) to its
# among-individual SD, then disattenuate the correlation for that
# measurement error (Spearman 1904): r_true = r_obs / sqrt(rel_x * rel_y),
# where reliability = 1 - mean(SE)^2 / var(estimate).
topt_edA <- df_final_merged %>%
  filter(!is.na(Topt), !is.na(Ed_A), !is.na(Topt_SE), !is.na(Ed_A_SE))

r_obs <- cor(topt_edA$Topt, topt_edA$Ed_A)

pct_SE_Topt <- mean(topt_edA$Topt_SE) / sd(topt_edA$Topt) * 100
pct_SE_ED   <- mean(topt_edA$Ed_A_SE)   / sd(topt_edA$Ed_A)   * 100

rel_Topt <- 1 - mean(topt_edA$Topt_SE)^2 / var(topt_edA$Topt)
rel_ED   <- 1 - mean(topt_edA$Ed_A_SE)^2   / var(topt_edA$Ed_A)

r_corrected <- r_obs / sqrt(rel_Topt * rel_ED)

cat(sprintf(
  "n = %d\nmean SE as %% of among-individual SD: Topt = %.1f%%, Ed_A = %.1f%%\nuncorrected r2 = %.3f\ncorrected r2   = %.3f\n",
  nrow(topt_edA), pct_SE_Topt, pct_SE_ED, r_obs^2, r_corrected^2
))


####################Fig S4  ###########
library(dplyr)
library(ggplot2)

# Exposure durations
exposure_times <- 1:100000

# Calculate correlation at each exposure duration
cor_results <- lapply(exposure_times, function(t) {
  
  CTmax <- df_final_merged$CTmax_1min -
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
ggsave("./figures_updated/FigS4.pdf", last_plot(),
       width = 8, height = 6, units = "in")

#Crossing temp
cor_results %>%
  mutate(
    sign = sign(r),
    previous_sign = lag(sign)
  ) %>%
  filter(sign != previous_sign) %>%
  select(time, r, previous_sign, sign)

################# Fig S5 ##################
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
ggsave("./figures_updated/FigS5.pdf", combined_plot,
       width = 8, height = 7, units = "in")


##################Fig S6 PCA ##################
df_pca <- df_final_merged %>%
  select(Tbreadth, Amax, Topt, Ed_A, T50, Ea_PSII, lnA_PSII, z, CTmax_1min, gsw_mean_above30, population_category) %>%
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


table_s3 <- as.data.frame(pca_res$rotation) #Table S3
table_s3
write.csv(table_s3, "./figures_updated/TableS3.csv")

summary(pca_res)

ggsave("./figures_updated/FigS6.pdf", pca_plot,
       width = 8, height = 8, units = "in")

############################### Excel stat summaries: Fig 4-7 ##############################
# One workbook per figure, containing the exact statistics annotated on
# the plot (and/or quoted in the Results text) for that figure. Each is
# recomputed here from df_final_merged directly, rather than reusing
# objects from earlier sections, so this block stays correct even if
# variable names above (e.g. "vars", "cor_mat") get reassigned later in
# the script.
library(openxlsx)
library(broom)

vif_to_df <- function(model, label) {
  v <- vif(model)
  data.frame(model = label, term = names(v), VIF = round(as.numeric(v), 3))
}

pairwise_corr_table <- function(data, vars) {
  sub <- data %>% dplyr::select(all_of(vars))
  rc <- Hmisc::rcorr(as.matrix(sub))
  pairs <- combn(vars, 2, simplify = FALSE)
  purrr::map_dfr(pairs, function(pr) {
    a <- pr[1]; b <- pr[2]
    data.frame(
      Variable_1 = a,
      Variable_2 = b,
      n = rc$n[a, b],
      r = round(rc$r[a, b], 4),
      r2 = round(rc$r[a, b]^2, 4),
      p = rc$P[a, b],
      significant = rc$P[a, b] < 0.05
    )
  })
}

# ---- Fig 4: pairwise correlations among Amax, Tbreadth, Topt, gs ----
fig4_corr <- pairwise_corr_table(
  df_final_merged, c("Amax", "Tbreadth", "Topt", "gsw_mean_above30")
)

wb4 <- createWorkbook()
addWorksheet(wb4, "Pairwise correlations")
writeData(wb4, "Pairwise correlations", fig4_corr)
saveWorkbook(wb4, "./figures_updated/Fig4_stats.xlsx", overwrite = TRUE)

# ---- Fig 5: gs ~ Tleaf (majority-supported range + mixed model) ----
fig5_range <- data.frame(
  Description = "Leaf temperature range where >=70% of curves have observations",
  Tmin_C = majority_range$Tmin,
  Tmax_C = majority_range$Tmax
)

mm_coefs <- as.data.frame(summary(mod_gsw_mixed)$coefficients)
mm_coefs <- data.frame(term = rownames(mm_coefs), mm_coefs, row.names = NULL)

wb5 <- createWorkbook()
addWorksheet(wb5, "Temperature range")
writeData(wb5, "Temperature range", fig5_range)
addWorksheet(wb5, "Mixed model gs~Tleaf")
writeData(wb5, "Mixed model gs~Tleaf", mm_coefs)
saveWorkbook(wb5, "./figures_updated/Fig5_stats.xlsx", overwrite = TRUE)

# ---- Fig 6: Amax ~ predictor + gs (partial regression) ----
fig6_partial <- bind_rows(
  get_partial_stats("Tbreadth", df_partial) %>% mutate(panel = "a: Tbreadth"),
  get_partial_stats("Topt", df_partial)      %>% mutate(panel = "b: Topt")
) %>% select(panel, predictor, n, pvalue, partial_R2)

fig6_full_coefs <- bind_rows(
  broom::tidy(mod_breadth) %>% mutate(model = "Amax ~ Tbreadth + gs"),
  broom::tidy(mod_Topt)    %>% mutate(model = "Amax ~ Topt + gs")
) %>% select(model, everything())

fig6_full_fit <- bind_rows(
  broom::glance(mod_breadth) %>% mutate(model = "Amax ~ Tbreadth + gs"),
  broom::glance(mod_Topt)    %>% mutate(model = "Amax ~ Topt + gs")
) %>% select(model, everything())

fig6_vif <- bind_rows(
  vif_to_df(mod_breadth, "Amax ~ Tbreadth + gs"),
  vif_to_df(mod_Topt, "Amax ~ Topt + gs")
)

wb6 <- createWorkbook()
addWorksheet(wb6, "Partial regression stats")
writeData(wb6, "Partial regression stats", fig6_partial)
addWorksheet(wb6, "Full model coefficients")
writeData(wb6, "Full model coefficients", fig6_full_coefs)
addWorksheet(wb6, "Full model fit")
writeData(wb6, "Full model fit", fig6_full_fit)
addWorksheet(wb6, "VIF")
writeData(wb6, "VIF", fig6_vif)
saveWorkbook(wb6, "./figures_updated/Fig6_stats.xlsx", overwrite = TRUE)

# ---- Fig 7: Tbreadth/Ed,A ~ Topt + T50 + gs (partial regression) ----
get_partial_stats_fig7 <- function(response, focal_predictor, predictors, data) {
  df <- data %>% dplyr::select(all_of(c(response, predictors))) %>% na.omit()
  full_mod <- lm(reformulate(predictors, response = response), data = df)
  full_R2 <- summary(full_mod)$r.squared
  control_vars <- predictors[predictors != focal_predictor]
  res_Y <- resid(lm(reformulate(control_vars, response = response), data = df))
  res_X <- resid(lm(reformulate(control_vars, response = focal_predictor), data = df))
  partial_s <- summary(lm(res_Y ~ res_X))
  data.frame(
    response = response, focal = focal_predictor,
    n = nrow(df),
    pvalue = coef(partial_s)[2, 4],
    partial_R2 = partial_s$r.squared,
    full_R2 = full_R2
  )
}

preds7 <- c("Topt", "T50", "gsw_mean_above30")
fig7_partial <- bind_rows(
  get_partial_stats_fig7("Tbreadth", "Topt", preds7, df_final_merged) %>% mutate(panel = "a: Tbreadth~Topt"),
  get_partial_stats_fig7("Tbreadth", "T50",   preds7, df_final_merged) %>% mutate(panel = "b: Tbreadth~T50"),
  get_partial_stats_fig7("Ed_A", "Topt", preds7, df_final_merged)        %>% mutate(panel = "c: Ed,A~Topt"),
  get_partial_stats_fig7("Ed_A", "T50",   preds7, df_final_merged)        %>% mutate(panel = "d: Ed,A~T50")
) %>% select(panel, response, focal, n, pvalue, partial_R2, full_R2)

fig7_full_coefs <- bind_rows(
  broom::tidy(vif_breadth_full) %>% mutate(model = "Tbreadth ~ Topt + T50 + gs"),
  broom::tidy(vif_EdA_full)     %>% mutate(model = "Ed,A ~ Topt + T50 + gs")
) %>% select(model, everything())

fig7_full_fit <- bind_rows(
  broom::glance(vif_breadth_full) %>% mutate(model = "Tbreadth ~ Topt + T50 + gs"),
  broom::glance(vif_EdA_full)     %>% mutate(model = "Ed,A ~ Topt + T50 + gs")
) %>% select(model, everything())

fig7_vif <- bind_rows(
  vif_to_df(vif_breadth_full, "Tbreadth ~ Topt + T50 + gs"),
  vif_to_df(vif_EdA_full, "Ed,A ~ Topt + T50 + gs")
)

wb7 <- createWorkbook()
addWorksheet(wb7, "Partial regression stats")
writeData(wb7, "Partial regression stats", fig7_partial)
addWorksheet(wb7, "Full model coefficients")
writeData(wb7, "Full model coefficients", fig7_full_coefs)
addWorksheet(wb7, "Full model fit")
writeData(wb7, "Full model fit", fig7_full_fit)
addWorksheet(wb7, "VIF")
writeData(wb7, "VIF", fig7_vif)
saveWorkbook(wb7, "./figures_updated/Fig7_stats.xlsx", overwrite = TRUE)

