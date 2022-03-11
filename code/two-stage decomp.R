require(Hmisc)

Theil <- function(y) {
  
  y_ave <- mean(y)
  N <- length(y)
  
  Ttotal = (1 / N) * sum( (y / y_ave) * log(y / y_ave) )
  
  return(Ttotal)
}

Theil_w <- function(y, pop) {
  
  y_ave_w <- weighted.mean(y, pop)
  N <- sum(pop)
  
  Ttotal =  sum( (pop / N) * (y / y_ave_w) * log(y / y_ave_w) )
  
  return(Ttotal)
}

CV <- function(y) {
  
  y_ave <- mean(y)
  
  var <- var(y)
  
  CVtotal = var / y_ave^2
  
  return(CVtotal)
  
}

CV_w <- function(y, pop) {
  
  require(Hmisc)
  
  y_ave_w <- weighted.mean(y, pop)
  
  var_w <- wtd.var(y, pop)
  
  CVtotal = var_w / y_ave_w^2
  
  return(CVtotal)
  
}

sd_w <- function(y, pop) {
  sd_w=sqrt(wtd.var(y, pop))
  return(sd_w)
}

onestage.Theil.decomp <- function(df, macroregion, ple) {
  
  y_ave <- mean(df[[ple]], na.rm = T)
  
  N <- length(df[[ple]])

  Ttotal = Theil(df[[ple]])
  
  by_macroregion <- split(df, df[[macroregion]])
  
  Ni <- sapply(by_macroregion, function(x){length(x[[ple]])})
  yi_ave <- sapply(by_macroregion, function(x){mean(x[[ple]], na.rm = T)})

  Tmi <- sapply(by_macroregion, function(x){
    
    Theil(x[[ple]])
  })
  
  Twithin <-  sum( (Ni/N) * (yi_ave/y_ave) * Tmi)
  Tbetween <-  sum( (Ni/N) * (yi_ave/y_ave) * log( yi_ave / y_ave ) )
  
  result <- data.frame(Tbetween = Tbetween, Twithin = Twithin, Ttotal = Ttotal, Additive = Twithin + Tbetween)
  
  return(result)
}

onestage.Theil.decomp.w <- function(df, macroregion, ple, pop) {
  
  y_ave <- weighted.mean(df[[ple]], df[[pop]])
  
  N <- sum(df[[pop]])
  
  Ttotal = Theil_w(df[[ple]], df[[pop]])
  
  by_macroregion <- split(df, df[[macroregion]])
  
  Ni <- sapply(by_macroregion, function(x){sum(x[[pop]])})
  
  yi_ave <- sapply(by_macroregion, function(x){weighted.mean(x[[ple]], x[[pop]])})
  
  Tmi <- sapply(by_macroregion, function(x){
    
    Theil_w(x[[ple]], x[[pop]])
  })
  
  Twithin <-  sum( (Ni/N) * (yi_ave/y_ave) * Tmi)
  Tbetween <-  sum( (Ni/N) * (yi_ave/y_ave) * log( yi_ave / y_ave ) )
  
  result <- data.frame(Tbetween = Tbetween, Twithin = Twithin, Ttotal = Ttotal, Additive = Twithin + Tbetween)
  
  return(result)
}

onestage.CV.decomp <- function(df, macroregion, ple) {
  
  y_ave <- mean(df[[ple]])
  N <- length(df[[ple]])
  
  CVtotal = CV(df[[ple]])
  
  by_macroregion <- split(df, df[[macroregion]])
  
  yi_ave <- sapply(by_macroregion, function(x){mean(x[[ple]])})
  Ni <- sapply(by_macroregion, function(x){length(x[[ple]])})
  
  CVi <- sapply(by_macroregion, function(x){
    CV(x[[ple]])
  })
  
  print(CVi)
  print(Ni/N)
  print(yi_ave/y_ave)
  
  CVwithin <-  sum( (Ni/N) * (yi_ave/y_ave)^2 * CVi)
  CVbetween <-  sum( (Ni/N) * (yi_ave/y_ave)^2 * (1- (y_ave / yi_ave) ) )
  
  result <- data.frame(CVbetween = CVbetween, CVwithin = CVwithin, CVtotal = CVtotal, Additive = CVwithin + CVbetween)
  
  return(result)
}

onestage.var.decomp <- function(df, macroregion, ple) {
  
  y_ave <- mean(df[[ple]])
  N <- length(df[[ple]])
  
  Vartotal = var(df[[ple]])
  
  by_macroregion <- split(df, df[[macroregion]])
  
  vari <- sapply(by_macroregion, function(x){var(x[[ple]])})
  Ni <- sapply(by_macroregion, function(x){length(x[[ple]])})
  yi_ave <- sapply(by_macroregion, function(x){mean(x[[ple]])})
  

  byNMSw = (Ni/N) * vari
  
  NMSw = byNMSw[1]
  OMSw = byNMSw[2]
  
  byNMSb = (Ni/N) * (yi_ave-y_ave)^2
  
  NMSb = byNMSb[1]
  OMSb = byNMSb[2]
  
  Varwithin <-  sum( (Ni/N) * vari)
  Varbetween <-  sum( (Ni/N) * (yi_ave-y_ave)^2 )
  
  result <- data.frame(NMSw, OMSw, Varwithin = Varwithin, NMSb, OMSb, Varbetween = Varbetween,  Vartotal = Vartotal, Additive = Varwithin + Varbetween)
  
  return(result)
}

onestage.var.decomp.w <- function(df, macroregion, ple, pop) {
  
  y_ave <- weighted.mean(df[[ple]], df[[pop]])
  N <- sum(df[[pop]])
  
  Vartotal = wtd.var(df[[ple]], df[[pop]])
  
  by_macroregion <- split(df, df[[macroregion]])
  
  vari <- sapply(by_macroregion, function(x){wtd.var(x[[ple]], x[[pop]])})
  Ni <- sapply(by_macroregion, function(x){sum(x[[pop]])})
  yi_ave <- sapply(by_macroregion, function(x){weighted.mean(x[[ple]], x[[pop]])})
  
  byNMSw = (Ni/N) * vari
  
  NMSw = byNMSw[1]
  OMSw = byNMSw[2]
  
  byNMSb = (Ni/N) * (yi_ave-y_ave)^2
  
  NMSb = byNMSb[1]
  OMSb = byNMSb[2]
  
  print(NMSb)
  
  Varwithin <-  sum( (Ni/N) * vari)
  Varbetween <-  sum( (Ni/N) * (yi_ave-y_ave)^2 )
  
  result <- data.frame(NMSw, OMSw, Varwithin = Varwithin, NMSb, OMSb, Varbetween = Varbetween,  Vartotal = Vartotal, Additive = Varwithin + Varbetween)
  
  return(result)
}

twostage.Theil.decomp <- function(df, country, macroregion, ple) {
  
  y_ave <- mean(df[[ple]], na.rm = T)
  N <- length(df[[ple]])
  
  Ttotal = Theil(df, ple)[[1]]
  
  by_macroregion <- split(df, df[[macroregion]])
  
  yi_ave <- sapply(by_macroregion, function(x){mean(x[[ple]], na.rm = T)})
  Ni <- sapply(by_macroregion, function(x){length(x[[ple]])})
  
  Tbmr <- sum( (yi_ave/y_ave) * (Ni/N) * log(yi_ave/y_ave) )
  
  by_country <- split(df, df[,country])
  
  yij_ave <- sapply(by_country, function(x){mean(x[[ple]], na.rm = T)})
  Nij <- sapply(by_country, function(x){length(x[[ple]])})
  
  Tij <- sapply(by_country, function(x){
    Theil(x, ple)[[1]]
  })
  
  Tpi <- sapply(by_macroregion, function(x){
    by_country2 <- split(x, x[,country])
    Yij <- sapply(by_country2, function(x){mean(x[[ple]], na.rm = T)})
    Nij <- sapply(by_country2, function(x){length(x[[ple]])})
    sum( (Yij / mean(x[[ple]], na.rm = T)) * (Nij/length(x[[ple]])) * log(Yij/mean(x[[ple]], na.rm = T)) )
  })
  
  Tbc <- sum( (yi_ave/y_ave) * (Ni/N) * Tpi, na.rm = T)
  Twc <- sum( (yij_ave/y_ave) * (Nij/N) * Tij, na.rm = T)
  
  result <- data.frame(Tbmr = Tbmr, Tbc = Tbc, Twc = Twc, Ttotal = Ttotal, Tsum = (Tbmr+Tbc+Twc), N = N)
  
  return(result)
}

twostage.CV.decomp <- function(df, country, macroregion, ple) {
  
  y_ave <- mean(df[[ple]], na.rm = T)
  N <- length(df[[ple]])
  
  CVtotal = CV(df, ple)[[1]]
  
  by_macroregion <- split(df, df[[macroregion]])
  
  yi_ave <- sapply(by_macroregion, function(x){mean(x[[ple]], na.rm = T)})
  Ni <- sapply(by_macroregion, function(x){length(x[[ple]])})
  
  CVbmr <- sum( (Ni/N) * (yi_ave/y_ave)^2 * (1- (y_ave / yi_ave) ) )
  
  by_country <- split(df, df[,country])
  
  yij_ave <- sapply(by_country, function(x){mean(x[[ple]], na.rm = T)})
  Nij <- sapply(by_country, function(x){length(x[[ple]])})
  
  CVij <- sapply(by_country, function(x){
    CV(x, ple)[[1]]
  })
  
  CVpi <- sapply(by_macroregion, function(x){
    by_country2 <- split(x, x[,country])
    Yij <- sapply(by_country2, function(x){mean(x[[ple]], na.rm = T)})
    Nij <- sapply(by_country2, function(x){length(x[[ple]])})
    sum( (Nij/length(x[[ple]])) * (Yij / mean(x[[ple]], na.rm = T))^2 * (1-mean(x[[ple]], na.rm = T)/Yij) )
  })
  
  CVbc <- sum( (Ni/N) * (yi_ave/y_ave)^2 * CVpi, na.rm = T)
  CVwc <- sum( (Nij/N) * (yij_ave/y_ave)^2 * CVij, na.rm = T)
  
  result <- data.frame(CVbmr = CVbmr, CVbc = CVbc, CVwc = CVwc, CVtotal = CVtotal, CVsum = (CVbmr+CVbc+CVwc), N = N)
  
  return(result)
}