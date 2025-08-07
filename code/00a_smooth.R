rm(list = ls())

require(tidyverse)
library(MortalitySmooth) 
source("code/mort_smooth.R")


CZE = readRDS('temp/CZE_for_smooth.RDS')%>%
  rename(Dx = Deaths, Ex = Pop, Age = AgeGroup)

EST_all = readRDS('temp/EST_NUTS3_all_for_smooth.RDS')%>%
  rename(Dx = Deaths, Ex = Pop, Age = AgeGroup)

LTU = readRDS('temp/LTU_for_smooth.RDS')%>%
  rename(Dx = Deaths, Ex = Pop, Age = AgeGroup)

POL = readRDS('temp/POL_for_smooth.RDS')%>%
  rename(Dx = Deaths, Ex = Pop, Age = AgeGroup) %>%
  filter(Year >= 2006) ## The smoothing algorithm requires the same number of age groups over time

ROU = readRDS('temp/ROU_for_smooth.RDS')%>%
  rename(Dx = Deaths, Ex = Pop, Age = AgeGroup)

SVK = readRDS('temp/SVK_for_smooth.RDS')%>%
  rename(Dx = Deaths, Ex = Pop, Age = AgeGroup)


CZE_smooth <- CZE %>%
  group_by(Country, RegionCode, Sex, CauseCategory) %>%
  group_modify(~mort2Dsmooth.fun(.))

EST_all_smooth <- EST_all %>%
  group_by(Country, RegionCode, Sex, CauseCategory) %>%
  group_modify(~mort2Dsmooth.fun(.))

LTU_smooth <- LTU %>%
  group_by(Country, RegionCode, Sex, CauseCategory) %>%
  group_modify(~mort2Dsmooth.fun(.))

POL_all_smooth <- POL %>%
  filter(CauseCategory == 'All') %>%
  group_by(Country, RegionCode, Sex, CauseCategory) %>%
  group_modify(~mort2Dsmooth.fun(.))

POL_cod_smooth <- POL %>%
  filter(CauseCategory != 'All') %>%
  group_by(Country, RegionCode, Sex, CauseCategory) %>%
  group_modify(~mort2Dsmooth.fun(., infant = FALSE))

POL_smooth = rbind(POL_all_smooth, POL_cod_smooth)

ROU_smooth <- ROU %>%
  group_by(Country, RegionCode, Sex, CauseCategory) %>%
  group_modify(~mort2Dsmooth.fun(.))

SVK_smooth <- SVK %>%
  group_by(Country, RegionCode, Sex, CauseCategory) %>%
  group_modify(~mort2Dsmooth.fun(.))


smooth = rbind(CZE_smooth, EST_all_smooth, LTU_smooth, POL_smooth, ROU_smooth, SVK_smooth)

smooth = smooth %>%
  mutate(AgePattern = case_when(Country == 'CZE' & CauseCategory == 'All' ~ '21',
                                  Country == 'CZE' & CauseCategory != 'All' ~ '6',
                                  Country == 'EST' & CauseCategory == 'All' ~ '21',
                                  Country == 'LTU' ~ '6',
                                  Country == 'POL' & CauseCategory == 'All' ~ '18',
                                  Country == 'POL' & CauseCategory != 'All' ~ 'POL5',
                                  Country == 'ROU' & CauseCategory == 'All' ~ '18',
                                  Country == 'ROU' & CauseCategory != 'All' ~ 'ROU6',
                                  Country == 'SVK' ~ '16'))


sdr = smooth %>%
  group_by(Country, RegionCode, Year, Sex, CauseCategory, Iteration) %>%
  summarise(smoothed_sdr = ESP_sdr.fun(logmx = smoothed_logmx, AgePattern = AgePattern),
            raw_sdr = ESP_sdr.fun(logmx = log(Dx/Ex), AgePattern = AgePattern))


## Add Estonia's cod

PS1D_wrap_numderiv <- function(x, y, e, lambdas, WEI = NULL,
                               levels=c(95,50), ...) {
  if (is.null(WEI)) WEI <- rep(1, length(y))
  # fit the smoother
  fit1d <- Mort1Dsmooth(x      = x,
                        y      = y,
                        offset = log(e),
                        lambda = lambdas,
                        method = 3,
                        w      = WEI,
                        ...)
  eta   <- fit1d$logmortality
  n     <- length(eta)
  dx    <- diff(x)
  if (max(dx) - min(dx) > 1e-8) 
    stop("x must be equally spaced for numeric differentiation")
  h <- dx[1]
  # central differences in the interior, forward/back at ends
  ETA1t <- numeric(n)
  ETA1t[1]     <- (eta[2]   - eta[1])     / h
  ETA1t[n]     <- (eta[n]   - eta[n-1])   / h
  ETA1t[2:(n-1)] <- (eta[3:n] - eta[1:(n-2)])/(2*h)
  
  # now quantiles just like deltasFUN()
  p.up  <- (100 - (100-levels[2])/2)/100
  p.low <- ((100-levels[2])/2)/100
  delta.t.up  <- quantile(ETA1t, probs=p.up,  na.rm=TRUE)
  delta.t.low <- quantile(ETA1t, probs=p.low, na.rm=TRUE)
  
  # now build a PSinfant‐style object
  out <- list(
    Y       = matrix(y, 1), e       = matrix(e, 1),
    lambdas = fit1d$lambda,   WEI     = matrix(WEI, 1),
    res     = matrix(fit1d$residuals, 1),
    ALPHAS  = matrix(fit1d$coefficients, 1),
    ETA     = matrix(eta, 1),
    ETA1a   = matrix(0, nrow=1, ncol=n), 
    ETA1t   = matrix(ETA1t, nrow=1, ncol=n),
    deltas  = list(delta.a.up=0, delta.a.low=0,
                   delta.t.up=delta.t.up,
                   delta.t.low=delta.t.low),
    S       = matrix(1, 1, n),
    kappas  = c(1,1),
    BCa     = NULL,
    BCt     = NULL,
    st.alphas = NULL
  )
  class(out) <- "PSinfant"
  return(out)
}

mort1Dsmooth.fun <- function(x) {
  frame <- arrange(x, Year)
  y <- as.integer(frame$Year)
  n = length(y)
  e = as.numeric(frame$Ex)
  z = as.numeric(frame$Dx)
  
  model <- Mort1Dsmooth(x=y, y=z, offset=log(e))
  
  ## extract deviance residuals
  res <- model$res
  ## replace res == NA (where Y1==0) with random residuals
  ## anyway they will be taken randomly in the next step
  whi0 <- which(is.na(res))
  res1 <- c(res)[-whi0]
  res[whi0] <- sample(res1, size=length(whi0))
  ## small correction of actual deaths when equal to zero
  ## to avoid computational issues when log is taken and we numerically find the root of the function
  Y11 <- z
  Y11[whi0] <- 1e-3
  
  ## bootstrapping procedure 
  
  ## number of instances
  n.sim <- 100
  ## empty array for bootstrap deaths
  Y1s <- matrix(0, n, n.sim)
  ## function to invert the bootstrap-deaths
  InvDev <- function(X, Z, C){X - Z * log(X) - C}
  
  eps      <- 1e-6
  big      <- 1e8
  max_tries <- 10
  
  ## resample deaths
  for(k in 1:n.sim){
    ## resample (with replacement) the residuals
    resS <- sample(x=res, size=n, replace = TRUE)  
    ## computing C =  res^2/2 + Dth - Dth*ln(Dth)
    C <- 1/2*(resS^2) + Y11 - Y11*log(Y11)
    ## recomputing the deaths
    low  <- ifelse(resS > 0, eps,       Y11)
    up   <- ifelse(resS > 0, Y11,       big)
    up   <- ifelse(low >= up, low + 1e-8, up)
    
    
    Y1s[, k] <- mapply(function(Z, C_ij, LOW, UP) {
      fl <- InvDev(LOW, Z, C_ij)
      fu <- InvDev(UP,   Z, C_ij)
      
      ## if bracket fails, fall back to original
      if (!is.finite(fl) || !is.finite(fu) || fl * fu > 0) {
        warning(sprintf("Bracket failed at Z=%.3g → returning original", Z))
        return(Z)
      }
      
      uniroot(InvDev,
              lower = LOW,
              upper = UP,
              tol   = 1e-5,
              Z     = Z,
              C     = C_ij)$root
    },
    Z    = Y11,
    C_ij = C,
    LOW  = low,
    UP   = up)
    
    cat("Completed bootstrap", k, "\n")

  }
  
  DELTAs <- list()
  
  for (k in 1:n.sim) {
    Y1.k     <- Y1s[, k]          # now a vector of length n
    FITinf.k <- PS1D_wrap_numderiv(x = y, y=Y1.k, e=e, lambdas=model$lambda)
    DELTAs[[k]] <- deltasFUN(FITinf.k)
  }
  
  ## if you still need a “forecast” step in one dimension:
  ETA.hatCs <- matrix(0, n, n.sim,
                      dimnames = list(Year = y,
                                      Sim  = seq_len(n.sim)))
  
  for (k in 1:n.sim) {
    ## fit the 1-D spline to the k-th bootstrap deaths
    fit1d.k <- Mort1Dsmooth(
      x      = y,
      y      = Y1s[, k],
      offset = log(e),
      ndx    = floor(length(y)/5),
      deg    = 3,
      pord   = 2,
      lambda = model$lambda,   # your “time” smoothness
      method = 3
    )
    
    ## extract its fitted log-mortality
    ETA.hatCs[, k] <- fit1d.k$logmortality
    
    ## (optional) diagnostic print
    cat("Sim", k, "\n")
  }
  
  out <- as.data.frame.table(ETA.hatCs, responseName = "smoothed_logmx")
  names(out) <- c("Year", "Iteration", "smoothed_logmx")
  out$Year      <- as.integer(as.character(out$Year))
  out$Iteration <- as.integer(as.character(out$Iteration))
  
  out <- left_join(out, frame, by = "Year")
  
}

EST_cod = readRDS('temp/EST_NUTS3_cod_rates.rds') %>%
  group_by(Country, RegionCode, Year, Sex, CauseCategory) %>%
  mutate(Dx = (sdr/1000)*Pop) %>%
  rename(Ex = Pop) %>%
  group_by(Country, RegionCode, Sex, CauseCategory) %>%
  group_modify(~mort1Dsmooth.fun(.))


smooth = rbind(smooth, EST_cod)

saveRDS(smooth, 'temp/smooth.RDS')

EST_cod_sdr = EST_cod %>%
  mutate(smoothed_sdr = 1000*exp(smoothed_logmx),
            raw_sdr = sdr)

sdr = rbind(sdr, EST_cod_sdr)

saveRDS(sdr, 'temp/sdr.RDS')


# ##sanity checks of all smoothed sdrs
# 
# 
# ggplot(data = filter(sdr, Country == 'CZE'))+
#   geom_line(aes(x = Year, y = smoothed_sdr, color = RegionCode, linetype = Sex, group = interaction(RegionCode, Sex, Iteration, CauseCategory)))+
#   geom_point(aes(x = Year, y = raw_sdr, color = RegionCode, shape = Sex, group = interaction(RegionCode, Sex)))+
#   facet_wrap(CauseCategory ~ RegionCode, scales = 'free_y')
# 
# ggplot(data = filter(test1, Country == 'EST'))+
#   geom_line(aes(x = Year, y = smoothed_sdr, color = RegionCode, linetype = Sex, group = interaction(RegionCode, Sex, Iteration, CauseCategory)))+
#   geom_point(aes(x = Year, y = raw_sdr, color = RegionCode, shape = Sex, group = interaction(RegionCode, Sex)))+
#   facet_wrap(CauseCategory ~ RegionCode, scales = 'free_y')
# 
# ggplot(data = filter(sdr, Country == 'LTU'))+
#   geom_line(aes(x = Year, y = smoothed_sdr, color = RegionCode, linetype = Sex, group = interaction(RegionCode, Sex, Iteration, CauseCategory)))+
#   geom_point(aes(x = Year, y = raw_sdr, color = RegionCode, shape = Sex, group = interaction(RegionCode, Sex)))+
#   facet_wrap(CauseCategory ~ RegionCode, scales = 'free_y')
# 
# ggplot(data = filter(sdr, Country == 'POL' & RegionCode == 'PL414'))+
#   geom_line(aes(x = Year, y = smoothed_sdr, color = RegionCode, linetype = Sex, group = interaction(RegionCode, Sex, Iteration, CauseCategory)))+
#   geom_point(aes(x = Year, y = raw_sdr, color = RegionCode, shape = Sex, group = interaction(RegionCode, Sex)))+
#   facet_wrap(CauseCategory ~ Sex, scales = 'free_y')
# 
# ggplot(data = filter(sdr, Country == 'ROU'))+
#   geom_line(aes(x = Year, y = smoothed_sdr, color = RegionCode, linetype = Sex, group = interaction(RegionCode, Sex, Iteration, CauseCategory)))+
#   geom_point(aes(x = Year, y = raw_sdr, color = RegionCode, shape = Sex, group = interaction(RegionCode, Sex)))+
#   facet_wrap(CauseCategory ~ Sex, scales = 'free_y')
# 
# ggplot(data = filter(sdr, Country == 'SVK'))+
#   geom_line(aes(x = Year, y = smoothed_sdr, color = RegionCode, linetype = Sex, group = interaction(RegionCode, Sex, Iteration, CauseCategory)))+
#   geom_point(aes(x = Year, y = raw_sdr, color = RegionCode, shape = Sex, group = interaction(RegionCode, Sex)))+
#   facet_wrap(CauseCategory ~ RegionCode, scales = 'free_y')
# 



EST_LAU_all = readRDS('temp/EST_LAU_all_for_smooth.RDS')%>%
  rename(Dx = Deaths, Ex = Pop, Age = AgeGroup)
LTU_LAU = readRDS('temp/LTU_LAU_for_smooth.RDS')%>%
  rename(Dx = Deaths, Ex = Pop, Age = AgeGroup)
SVK_LAU = readRDS('temp/SVK_LAU_for_smooth.RDS')%>%
  rename(Dx = Deaths, Ex = Pop, Age = AgeGroup)


EST_all_smooth <- EST_LAU_all %>%
  group_by(Country, RegionCode, Sex, CauseCategory) %>%
  group_modify(~mort2Dsmooth.fun(.))

LTU_smooth <- LTU_LAU %>%
  group_by(Country, RegionCode, Sex, CauseCategory) %>%
  group_modify(~mort2Dsmooth.fun(.))

SVK_smooth <- SVK_LAU %>%
  group_by(Country, RegionCode, Sex, CauseCategory) %>%
  group_modify(~mort2Dsmooth.fun(., infant = FALSE))

smooth_LAU = rbind(EST_all_smooth, LTU_smooth, SVK_smooth)

smooth_LAU = smooth_LAU %>%
  mutate(AgePattern = case_when(Country == 'CZE' & CauseCategory == 'All' ~ '21',
                                Country == 'CZE' & CauseCategory != 'All' ~ '6',
                                Country == 'EST' & CauseCategory == 'All' ~ '21',
                                Country == 'LTU' ~ '6',
                                Country == 'POL' & CauseCategory == 'All' ~ '18',
                                Country == 'POL' & CauseCategory != 'All' ~ 'POL5',
                                Country == 'ROU' & CauseCategory == 'All' ~ '18',
                                Country == 'ROU' & CauseCategory != 'All' ~ 'ROU6',
                                Country == 'SVK' ~ '16'))

sdr_LAU = smooth_LAU %>%
  group_by(Country, RegionCode, Year, Sex, CauseCategory, Iteration) %>%
  summarise(smoothed_sdr = ESP_sdr.fun(logmx = smoothed_logmx, AgePattern = AgePattern),
            raw_sdr = ESP_sdr.fun(logmx = log(Dx/Ex), AgePattern = AgePattern))


EST_LAU_cod_smooth = readRDS('temp/EST_LAU_cod_rates.rds') %>%
  group_by(Country, RegionCode, Year, Sex, CauseCategory) %>%
  mutate(Dx = (sdr/1000)*Pop) %>%
  rename(Ex = Pop) %>%
  group_by(Country, RegionCode, Sex, CauseCategory) %>%
  group_modify(~mort1Dsmooth.fun(.))

smooth_LAU = rbind(smooth_LAU, EST_LAU_cod_smooth)

saveRDS(smooth_LAU, 'temp/smooth_LAU.RDS')

EST_LAU_cod_sdr = EST_LAU_cod_smooth %>%
  mutate(smoothed_sdr = 1000*exp(smoothed_logmx),
         raw_sdr = sdr)

sdr_LAU = rbind(sdr_LAU, EST_LAU_cod_sdr)

saveRDS(sdr_LAU, 'temp/sdr_LAU.RDS')

# ggplot(data = filter(sdr_LAU, Country == 'EST'))+
#   geom_line(aes(x = Year, y = smoothed_sdr, color = RegionCode, linetype = Sex, group = interaction(RegionCode, Sex, Iteration, CauseCategory)))+
#   geom_point(aes(x = Year, y = raw_sdr, color = RegionCode, shape = Sex, group = interaction(RegionCode, Sex)), size = 0.5)+
#   facet_wrap(CauseCategory ~ RegionCode, scales = 'free_y')
# 
# ggplot(data = filter(sdr_LAU, Country == 'LTU' & CauseCategory == 'Other' & Iteration == 10))+
#   geom_line(aes(x = Year, y = smoothed_sdr, color = RegionCode, linetype = Sex, group = interaction(RegionCode, Sex, Iteration, CauseCategory)))+
#   geom_point(aes(x = Year, y = raw_sdr, color = RegionCode, shape = Sex, group = interaction(RegionCode, Sex)), size = 0.5)+
#   facet_wrap(CauseCategory ~ RegionCode, scales = 'free_y')
# 
# ggplot(data = filter(sdr_LAU, Country == 'SVK'))+
#   geom_line(aes(x = Year, y = smoothed_sdr, color = RegionCode, linetype = Sex, group = interaction(RegionCode, Sex, Iteration, CauseCategory)))+
#   geom_point(aes(x = Year, y = raw_sdr, color = RegionCode, shape = Sex, group = interaction(RegionCode, Sex)), size = 0.5)+
#   facet_wrap(CauseCategory ~ RegionCode, scales = 'free_y')

