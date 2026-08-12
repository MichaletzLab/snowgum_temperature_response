schoolfield <- function(temp, J_ref, E, E_D, T_opt, T_ref = 25) {
  k <- 8.62e-5  # Boltzmann's constant
  temp_K <- temp + 273.15
  T_ref_K <- T_ref + 273.15
  T_opt_K <- T_opt + 273.15
  
  # force numeric values (avoid named scalars causing recycling issues)
  J_ref <- as.numeric(J_ref)
  E     <- as.numeric(E)
  E_D   <- as.numeric(E_D)
  T_opt <- as.numeric(T_opt)
  T_ref <- as.numeric(T_ref)
  
  numerator <- J_ref * exp(E * (1/(k*T_ref_K) - 1/(k*temp_K)))
  denominator <- 1 + (E / (E_D - E)) * exp((E_D / k) * (1/T_opt_K - 1/temp_K))
  
  return(numerator / denominator)
}


fit_mod_schoolfield2 <- function(Data, curve_id_col = "curveID", 
                                 x = "Tleaf", y = "A", T_ref = 25, fitted_models = TRUE) {
  models <- c("schoolfield")
  aic_values <- list()
  fits <- list()  # Store model fits
  Data$T_ref = T_ref
  
  # Create output directory if it doesn't exist
  if (!dir.exists("SS_fits")) dir.create("SS_fits")
  
  # Standardize names of x and y
  if (x != "Tleaf") {
    Data$Tleaf = Data[[x]] 
    Data[[x]] = NULL      
  }
  if (y != "A") {
    Data$A = Data[[y]] 
    Data[[y]] = NULL      
  }
  
  Data = subset(Data, A > 0)
  Data$curve_type = "Schoolfield"
  Data <- Data[order(Data[[curve_id_col]]), ]
  
  # Check for missing values
  if (anyNA(Data$Tleaf) || anyNA(Data$A)) {
    stop("Missing values detected in the input data.")
  }
  
  schoolfield_results <- list()
  
  ######################
  # Curve fitting loop #
  ######################
  for (curve_id in unique(Data$curveID)) {
    curve_data <- Data %>%
      filter(curveID == curve_id) %>%
      select(Tleaf, A)
    
    print(paste("Curve ID:", curve_id, ", Number of observations:", nrow(curve_data)))
    
    plot <- ggplot(curve_data, aes(x = Tleaf, y = A)) +
      geom_point() +
      theme_classic() +
      ggtitle(paste("Curve ID:", curve_id))
    
    for (model in models) {
      if (model == "schoolfield") {
        fit <- tryCatch({
          nls.multstart::nls_multstart(
            A ~ schoolfield(temp = Tleaf, J_ref, E, E_D, T_opt, T_ref = 25),
            data = curve_data,
            iter = c(2,2,2,2),
            start_lower = c(J_ref = 0, E = 0, E_D = 0.2, T_opt = 0),
            start_upper = c(J_ref = 20, E = 2, E_D = 5, T_opt = 50),
            supp_errors = 'Y',
            na.action = na.omit,
            lower = c(J_ref = -10, E = 0, E_D = 0, T_opt = 0)
          )
        }, error = function(e) {
          print(paste("Fitting failed for Curve ID:", curve_id))
          print(e)
          return(NULL)
        })
        
        predicted_schoolfield <- fitted(fit)
        
        if (!is.null(fit)) {
          r_sq = 1 - sum(resid(fit)^2) / sum((curve_data$A - mean(curve_data$A))^2)
        } else {
          r_sq <- NA
        }
      }
      
      if (!is.null(fit)) {
        J_ref <- coef(fit)["J_ref"]
        J_ref_SE <- summary(fit)$coefficients[,'Std. Error']['J_ref']
        
        E <- coef(fit)["E"]
        E_SE <- summary(fit)$coefficients[,'Std. Error']['E']
        
        E_D <- coef(fit)["E_D"]
        E_D_SE <- summary(fit)$coefficients[,'Std. Error']['E_D']
        
        T_opt <- coef(fit)["T_opt"]
        
        # Test extraction of within-curve parameter variances
        test_vcov <- vcov(fit)
        
        # Within-curve variances
        J_ref_var <- test_vcov["J_ref", "J_ref"]
        E_var     <- test_vcov["E", "E"]
        E_D_var   <- test_vcov["E_D", "E_D"]
        T_opt_var <- test_vcov["T_opt", "T_opt"]
        
        # Within-curve covariance
        Topt_ED_cov <- test_vcov["T_opt", "E_D"]
        
        # Within-curve correlation
        param_cor <- cov2cor(test_vcov)
        Topt_ED_cor <- param_cor["T_opt", "E_D"]
        
        # Existing Amax calculation
        Amax <- predict(fit, newdata = data.frame(Tleaf = T_opt))
        
        T_opt_SE <- summary(fit)$coefficients[,'Std. Error']['T_opt']
        AIC <- AIC(fit)
        
        T_extrap <- seq(min(curve_data$Tleaf) - 5, max(curve_data$Tleaf) + 5, by = 0.1)
        pred_extrap <- schoolfield(
          temp = T_extrap,
          J_ref = coef(fit)["J_ref"],
          E     = coef(fit)["E"],
          E_D   = coef(fit)["E_D"],
          T_opt = coef(fit)["T_opt"]
        )
        
        # --- 80% breadth with fixed 15°C baseline ---
        
        # --- 90% of Amax threshold (no baseline) ---
        
        threshold_90 <- 0.95 * Amax
        
        # All predicted temperatures above threshold
        T_above_90 <- T_extrap[pred_extrap >= threshold_90]
        
        if (length(T_above_90) > 1) {
          Tbreadth_90 <- range(T_above_90)
          Breadth_90_value <- diff(Tbreadth_90)
        } else {
          Tbreadth_90 <- c(NA_real_, NA_real_)
          Breadth_90_value <- NA_real_
        }
        
        pred_df <- data.frame(
          Tleaf = T_extrap,
          A = pred_extrap
        )
        
        param_text <- paste0(
          "Topt = ", round(T_opt,1), "°C\n",
          "Amax = ", round(Amax,2), "\n",
          "J_ref = ", round(J_ref,2), "\n",
          "E = ", round(E,2), "\n",
          "E_D = ", round(E_D,2), "\n",
          "Tbreadth_95 = ", round(Tbreadth_90[1],1), "–", round(Tbreadth_90[2],1), "°C"
        )
        
        plot <- ggplot(curve_data, aes(x = Tleaf, y = A)) +
          geom_point(size = 1, color = "gray70") +
          geom_line(data = pred_df, aes(x = Tleaf, y = A), size = 1.2, color = "black", linetype = "dashed") +
          # Horizontal dashed line across full Tbreadth
          geom_segment(aes(
            x = Tbreadth_90[1],
            xend = Tbreadth_90[2],
            y = threshold_90,
            yend = threshold_90
          ),
          color = "#FFAD03", linetype = "dashed", size = 1)+
          geom_point(data = data.frame(Tleaf = T_opt, A = Amax),
                     aes(x = Tleaf, y = A),
                     color = "#990066", size = 3)+
          annotate(
            "text",
            x = mean(T_extrap),
            y = min(pred_extrap) + 0.05*(max(pred_extrap) - min(pred_extrap)),
            label = param_text,
            hjust = 0.5,
            vjust = 0,
            color = "black") +
          theme_classic() +
          ggtitle(paste("Curve ID:", curve_id))
      }
      
      clean1 <- function(x) { unname(x[1]) }
      
      schoolfield_results[[as.character(curve_id)]] <- data.frame(
        curveID           = curve_id,
        
        J_ref             = clean1(J_ref),
        J_ref_SE          = clean1(J_ref_SE),
        J_ref_var         = clean1(J_ref_var),
        
        E                 = clean1(E),
        E_SE              = clean1(E_SE),
        E_var             = clean1(E_var),
        
        E_D               = clean1(E_D),
        E_D_SE            = clean1(E_D_SE),
        E_D_var           = clean1(E_D_var),
        
        T_opt             = clean1(T_opt),
        T_opt_SE          = clean1(T_opt_SE),
        T_opt_var         = clean1(T_opt_var),
        
        Topt_ED_cov       = clean1(Topt_ED_cov),
        Topt_ED_cor       = clean1(Topt_ED_cor),
        
        Amax              = clean1(Amax),
        AIC               = clean1(AIC),
        r_sq              = clean1(r_sq),
        breadth_90        = clean1(diff(Tbreadth_90))
      )
      
      
      if (!is.null(fit)) {
        fits[[model]] <- fit
      }
    }
    
    print(plot)
    ggsave(
      filename = file.path("SS_fits", paste0("curve_", curve_id, ".png")),
      plot = plot,
      width = 6,
      height = 5,
      dpi = 300
    )
  }
  
  # Combine results
  schoolfield_results <- schoolfield_results[!sapply(schoolfield_results, is.null)]
  schoolfield_results <- do.call(rbind, schoolfield_results)
  rownames(schoolfield_results) <- NULL
  
  return(schoolfield_results)
}

