#########################################################################################################
#                                Use NCD Function to Scan Genome Data         
######################################################################################################### 
## Date: 12/08/2019
## Involved: Miles
## Task: Use Bitarello's NCD1 function to scan population sepcific genome data.
## Note: The updates to this, when compared to the original are in the readme.

library(dplyr)
library(data.table)

## The Bitarello 2019 NCD1 function (no outgroup needed). 
## NOTE: The W here is the window size, and S is the step size.
NCD1 <- function(X, W = 7500, S = 1500) {
  ## With my alterations to start at position 1, rather than the first variant site. 
  windows_dt <-
    data.table(POS = seq(1, X[nrow(X), POS], S))[
      , POS2 := POS + W][
        -length(POS)]

    print (paste0('Finished setting up coordinates for chr ', unique(X$CHR)))

    ## ordering the windows dataframe by POS1
    setkey(windows_dt, POS, POS2) #
    ## Creating new column called POS2 that is the same as POS1 since it's a single site.
    X[, POS2 := POS]
    ## Same as above, and a prerequisite for using foverlaps
    setkey(X, POS, POS2)
  
    ## Calculates the number sites in each of the windows in windows_dt
   X_windows <-
    foverlaps(X, windows_dt, type = "within", nomatch = 0L)[ #this is not ideal but for now it's fine
      , window := .GRP, by = .(POS, POS2)][
        order(window, i.POS)][
          , .(Win.ID = paste(CHR, POS, POS2, sep = "|"), MAF)][,N_Raw:=.N,by = Win.ID]

  print (paste0('Finished selecting SNPs per window for chr ', unique(X$CHR)))

  ## Calculates the deviation from the balanced target frequency
  X_NCD <-
    X_windows[MAF != 0 & MAF != 1][
      , .(N_Raw= N_Raw,
          N_SNPs_cor = .N,
          NCD1_tf0.5 = sqrt(sum((MAF-0.5)^2)/.N),  
          NCD1_tf0.4 = sqrt(sum((MAF-0.4)^2)/.N),
          NCD1_tf0.3=sqrt(sum((MAF-0.3)^2)/.N)),
      by = Win.ID]

  print (paste0('NCD1 calculations done for chr ', unique(X$CHR)))

  ## Gets the number of informative SNPs
  unique_windows <- X_windows[, .(N_SNPs_cor = sum(MAF != 0 & MAF != 1)), by = Win.ID]
  unique_windows$Win.ID %>% gsub(pattern = "\\|.*", replacement = "") %>% as.numeric() %>% cbind(unique_windows, Chr = .) -> unique_windows
  unique_windows$Win.ID %>% gsub(pattern = "^[0-9]{1,2}\\|", replacement = "") %>%  gsub(pattern = "\\|.*", replacement = "")  %>% as.numeric() %>% cbind(unique_windows, Start = .) -> unique_windows
  unique_windows$Win.ID %>% gsub(pattern = "^[0-9]{1,2}\\|.*\\|", replacement = "") %>% as.numeric() %>% cbind(unique_windows, End = .) -> unique_windows

  setkey(unique_windows, Win.ID, N_SNPs_cor)
  setkey(X_NCD, Win.ID, N_SNPs_cor)

  unique(X_NCD[unique_windows])-> X_NCD2
  setorder(X_NCD2, Start)
	return(X_NCD2)
}
