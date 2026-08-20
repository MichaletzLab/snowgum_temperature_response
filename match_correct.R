match_correct = function(data) {

  return_data = c()

  for (i in unique(data$curveID)) {

    # Grab each curve in our dataset and figure out how long it is
    a = subset(data, curveID == i)
    len = dim(a)[1]

    # We need to do some kind of check here to make sure that there exists one
    # match test point at the end
    if(a$co2_adj[1] == a$co2_adj[len]) {
      warning(paste("Warning: curveID", i, "appears to have no match test point."))
    }

    # Interpolate match offset corrections from original match value and test point
    # assuming match accumulates linearly with time (check this assumption)
    correction_co2 = (a$time-a$time[1])*(a$co2_adj[len]-a$co2_adj[len-1])/(a$time[len]-a$time[1])
    correction_h2o = (a$time-a$time[1])*(a$h2o_adj[len]-a$h2o_adj[len-1])/(a$time[len]-a$time[1])

    # Add to data frame
    a$correction_co2 = correction_co2
    a$correction_h2o = correction_h2o

    # Change sample CO2 and H2O values using new correction
    b=a
    b$CO2_s = b$CO2_s + b$correction_co2
    b$H2O_s = b$H2O_s + b$correction_h2o

    # Recalculate gas exchange variables and drop test point
    b = b %>% calc_licor6800()
    b = b[1:len-1,]

    # Bind together
    return_data = bind_rows(return_data, b)

    #ggplot(data = a, aes(x = Tleaf, y = A)) + geom_point()
    #
    # irga_data =bind_rows(
    #   a %>% select(Tleaf, A) %>% mutate(Condition = "Uncorrected"),
    #   b %>% select(Tleaf, A) %>% mutate(Condition = "Corrected"))
    #
    # irga_wide =
    #   merge(
    #     a %>% select(Tleaf, A) %>% rename(A_corr = A),
    #     b %>% select(Tleaf, A) %>% rename(A_uncorr = A))
    #
    #
    #
    # ggplot(irga_data, aes(x = Tleaf, y =A, color=Condition)) +
    #   geom_point(size = 1) +
    #   #my_theme +
    #   xlab("Leaf Temperature (ºC)") +
    #   ylab("Assimilation rate (µmol/m²s)")
    #
    # #p2 = ggplot(irga_wide, aes(x = A_uncorr, A_corr)) +
    #   geom_point(size = 1) +
    #   geom_abline(slope = 1, lty = 2) +
    #   #my_theme +
    #   xlab("Assimilation rate, uncorrected (µmol/m²s)") +
    #   ylab("Assimilation rate, corrected (µmol/m²s)")
    #
    # grid.arrange(p1, p2, ncol = 1)
  }

  return(return_data)


}
