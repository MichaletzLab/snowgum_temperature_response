read_6800_with_BLC <- function(filename) {
  
  raw_input = readLines(filename)                # Read in raw datafile
  
  # Look for stability criteria and chop them out because they can mess things up
  #raw_input = raw_input[-grep("Stability Definition", raw_input)]
  
  # So here we need to grab the BLC parameters and then later attach them
  blc_data = data.frame(raw = raw_input[grep("blc_a", raw_input):grep("blc_Po", raw_input)]) %>% 
    separate(raw, into = c(NA, "blc_params"), sep = ":") %>% 
    separate(blc_params, into = c("name", "value"), sep = "\t") 
  blc_data$value = as.numeric(blc_data$value)
  blc_data = blc_data %>% spread(name, value)
  
  eb_data = data.frame(raw = raw_input[grep("LTConst", raw_input)]) %>% 
    separate(raw, into = c(NA, "EBparams"), sep = ":") %>% 
    separate(EBparams, into = c("name", "value"), sep = "\t")
  eb_data$value = as.numeric(eb_data$value)
  eb_data = eb_data %>% spread(name, value)
  
  data_start = grep("\\[Data\\]", raw_input) + 2 # Find where data starts (at "[Data]" line)
  data_end = grep("\\[Header\\]", raw_input) - 1 # Find where data ends (at "[Header]" line)
  data_end = c(data_end, length(raw_input))[-1]  # Add data end at end of file
  
  data_end = data_end[!(data_end < data_start[1])]
  
  compiled_data = c() # Initiate blank holder
  
  # Loop over each data set in the file
  for (i in 1:length(data_start)) {
    trimmed = raw_input[data_start[i]:data_end[i]]      # Grab data
    
    if(!is_empty(grep("Stability Definition", trimmed))) {
      trimmed = trimmed[-grep("Stability Definition", trimmed)]# Remove stability definitions
    }
    if(!is_empty(grep("^[^\t]*\t?[^\t]*$", trimmed))) {
      trimmed = trimmed[-grep("^[^\t]*\t?[^\t]*$", trimmed)] # Remove comments, hopefully
    }
    
    trimmed = trimmed[-2]                               # Remove units line
    current_data = read.csv(text = trimmed, sep = "\t") # Convert to dataframe
    compiled_data = rbind(compiled_data, current_data)  # Merge together into one 
  }
  compiled_data$filename = filename
  
  compiled_data = bind_cols(compiled_data, blc_data, eb_data)
  
  return(compiled_data)
}


recalc_blc = function(data) {
  
  a = data$blc_a
  b = data$blc_b
  c = data$blc_c
  d = data$blc_d
  e = data$blc_e
  x = data$Fan_speed*data$Pa / (1000*data$blc_Po)
  y = pmax(pmin(data$S, data$blc_maxS) , data$blc_minS)

  
  # If broadleaf geometry is specified...
  BLC = ifelse(data$Geometry == "0: Broadleaf",
         # Then apply broadleaf formula.
         a + b*x + c*x*y*y + d*x*y + e*x*x,
         # Otherwise, check if needle geometry is specified...
         ifelse(data$Geometry == "1: Needle",
                # Then apply needle value.
                3,
                # Otherwise, apply custom BLC value
                data$Custom
         )
  )

  data$gbw = BLC
}


# After BLC correction, we need to correct everything else which is downstream of BLC
# gsw -> gtc -> Ci, 
# TleafEB -> TleafCnd (if Ebal is on) -> everything which depends on Tleaf

calc_licor6800 = function(licor) {
  
  x = licor$Fan_speed*licor$Pa / (1000*licor$blc_Po)
  y = pmax(pmin(licor$S, licor$blc_maxS) , licor$blc_minS)
  
  licor %>% 
    mutate(
      gbw = ifelse(Geometry == "0: Broadleaf",
                   # Then apply broadleaf formula.
                   blc_a + blc_b*x + blc_c*x*y*y + blc_d*x*y + blc_e*x*x,
                   # Otherwise, check if needle geometry is specified...
                   ifelse(Geometry == "1: Needle",
                          # Then apply needle value.
                          3,
                          # Otherwise, apply custom BLC value
                          Custom
                   )
            ),
      E = Flow * CorrFact * (H2O_s-H2O_r)/(100*S*(1000-CorrFact*H2O_s)),
      A = Flow * CorrFact * (CO2_r-CO2_s*(1000-CorrFact*H2O_r)/(1000-CorrFact*H2O_s))/(100*S),
      Ca = CO2_s - ifelse(CorrFact>1, A*S*100 , 0),
      # Error again....
      Rabs = Qin * convert,
      VPcham = H2O_s * (Pa+`ΔPcham`)/1000,
      SVPcham = 0.61365 * exp(17.502*Tair/(240.97+Tair)),
      RHcham = VPcham/SVPcham*100,
      TleafEB = (Tair+(Rabs+2*0.95*0.0000000567*(((Tair+deltaTw)+273)^4 - (Tair+273)^4)-44100*E)/(1.84*29.3*gbw+8*0.95*0.0000000567*(Tair+273)^3)),
      TleafCnd = fT1*Tleaf + fT2*Tleaf2 + fTeb*TleafEB,
      SVPleaf = 0.61365*exp(17.502*TleafCnd/(240.97+TleafCnd)),
      VPDleaf = (SVPleaf-H2O_s*(Pa+`ΔPcham`)/1000),
      LatHFlux = -E*44100,
      SenHFlux = 2*29.3*gbw*0.92*(Tair-TleafCnd),
      NetTherm = 2*0.95*0.0000000567*(((Tair+deltaTw)+273)^4-(TleafCnd+273)^4),
      EBSum = Rabs+NetTherm+LatHFlux+SenHFlux,
      gtw = E*(1000-(1000*0.61365*exp(17.502*TleafCnd/(240.97+TleafCnd))/(Pa+`ΔPcham`)+H2O_s)/2)/(1000*0.61365*exp(17.502*TleafCnd/(240.97+TleafCnd))/(Pa+`ΔPcham`)-H2O_s),
      gsw = 2 / ((1/gtw - 1/gbw) + sign(gtw)*sqrt((1/gtw - 1/gbw)^2 + 4*K/((K+1)^2)*(2/(gtw*gbw) - 1/(gbw^2)))),
      gtc = 1 / ((K+1)/(gsw/1.6)+1/(gbw/1.37)) + K/((K+1)/(gsw/1.6) + K/(gbw/1.37)),
      Ci = ((gtc-E/2)*Ca-A)/(gtc+E/2),
      Pci = Ci*(Pa+`ΔPcham`)/1000,
      Pca = (CO2_s - ifelse(CorrFact>1, A*S*100/(Fan*Fan_speed), 0)) * (Pa+`ΔPcham`)/1000
    )
  
}




