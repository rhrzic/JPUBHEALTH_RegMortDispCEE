## Code based on paper "Smooth Constrained Mortality Forecasting" by Carlo G. Camarda DOI: 10.4054/DemRes.2019.41.38
## Revised by Rok Hrzic (r.hrzic (at) maastrichtuniversity.nl) in 2025

## function to build up B-splines and associated bases for derivatives
BsplineGrad <- function(x, xl, xr, ndx=NULL, deg, knots=NULL){
  if(is.null(knots)){
    dx <- (xr - xl)/ndx
    knots <- seq(xl - deg * dx, xr + deg * dx, by=dx)
    knots <- round(knots, 8)
  }else{
    knots <- knots
    dx <- diff(knots)[1]
  }
  P <- outer(x, knots, MortSmooth_tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff=deg+1)/(gamma(deg+1)*dx^deg)
  B <- (-1)^(deg + 1) * P %*% t(D)
  ##
  knots1 <- knots[-c(1,length(knots))]
  P <- outer(x, knots1, MortSmooth_tpower, deg-1)
  n <- dim(P)[2]
  D <- diff(diag(n),diff=deg)/(gamma(deg)*dx^(deg-1))
  BB <- ((-1)^(deg) * P %*% t(D))/dx
  D <- diff(diag(ncol(BB) + 1))
  C <- BB %*% D
  ##
  out <- list(dx=dx, knots=knots, B=B, C=C)
}


## function for estimating two-dimensional P-spline 
## addressing infant mortality (when infant=TRUE)
## for a given set of smoothing parameters lambdas
PSinfant <- function(Y, E, lambdas, WEI, infant=TRUE, verbose=FALSE){
  ## dimensions
  m <- nrow(Y)
  n <- ncol(Y)
  a <- 1:m
  t <- 1:n
  ## w/o age 0
  a0 <- a[-1]
  m0 <- m-1
  ## original offset
  OFF <- log(E)
  ## B-splines basis
  ## with infant-specialized coeff
  if(infant){
    ## over ages w/o age 0
    a0min <- min(a0)
    a0max <- max(a0)
    nda0 <- floor(m0/5)
    dega <- 3
    BCa0 <- BsplineGrad(a0, a0min, a0max, nda0, dega)
    Ba0 <- BCa0$B
    nba0 <- ncol(Ba0)
    ## adding infant-specific basis
    Ba <- cbind(0, Ba0)
    Ba <- rbind(c(1, rep(0,nba0)), Ba)
    nba <- ncol(Ba)
    ## basis for the derivatives
    ## over ages w/o age 0
    Ca0 <- BCa0$C
    ## adding infant-specific basis
    Ca <- cbind(0, Ca0)
    Ca <- rbind(c(-1,Ba[2,1:dega+1],rep(0,nba0-dega)),
                Ca)
  }else{
    amin <- min(a)
    amax <- max(a)
    nda <- floor(m/5)
    dega <- 3
    BCa <- BsplineGrad(a, amin, amax, nda, dega)
    Ba <- BCa$B
    nba <- ncol(Ba)
    ## basis for the derivatives
    Ca <- BCa$C
  }
  ## over years
  tmin <- min(t)
  tmax <- max(t)
  ndt <- floor(n/5)
  degt <- 3
  BCt <- BsplineGrad(t, tmin, tmax, ndt, degt)
  Bt <- BCt$B
  nbt <- ncol(Bt)
  ## basis for the derivatives
  Ct <- BCt$C
  ## weights for exposures=0
  WEI[E==0] <- 0
  
  ## tensor product of B-splines for the GLAM
  Ba1 <- kronecker(matrix(1, ncol=nba, nrow=1), Ba)
  Ba2 <- kronecker(Ba, matrix(1, ncol=nba, nrow=1))
  RTBa <- Ba1 * Ba2
  Bt1 <- kronecker(matrix(1, ncol=nbt, nrow=1), Bt)
  Bt2 <- kronecker(Bt, matrix(1, ncol=nbt, nrow=1))
  RTBt <- Bt1 * Bt2
  
  ## penalty terms
  ## over ages
  Da <- diff(diag(nba), diff=2)
  ## no penalization over age for age 0,
  ## with infant=TRUE
  if(infant){
    Da[1,1] <- 0
  }
  tDDa <- t(Da) %*% Da
  ## over years
  Dt <- diff(diag(nbt), diff=2)
  tDDt <- t(Dt) %*% Dt
  ## kronecker product of difference matrices
  Pa <- kronecker(diag(nbt), tDDa)
  Pt <- kronecker(tDDt, diag(nba))
  ## smoothing parameters
  lambda.a <- lambdas[1]
  lambda.t <- lambdas[2]
  P <- lambda.a * Pa + lambda.t * Pt
  ## data in vector for the starting values
  y <- c(Y)
  e <- c(E)
  wei <- c(WEI)
  off <- log(e)
  off0 <- off
  off0[wei==0] <- 100
  ## "other" offset in matrices
  OFF0 <- matrix(off0, m, n)
  ## starting coeff
  aa <- rep(a, n)
  tt <- rep(t, each = m)
  fit0 <- glm(round(y) ~ aa + tt + offset(off0), 
              family = poisson, weights = wei)
  etaGLM <- matrix(log(fit0$fitted) - c(off0), m, n)
  eta0 <- log((Y + 1)) - OFF0
  eta0[WEI == 0] <- etaGLM[WEI == 0]
  BBa <- solve(t(Ba) %*% Ba + 1e-06*diag(nba), t(Ba))
  BBt <- solve(t(Bt) %*% Bt + 1e-06*diag(nbt), t(Bt))
  alphas <- MortSmooth_BcoefB(BBa, BBt, eta0)
  ## Poisson iteration within a GLAM framework
  for(it in 1:20){
    eta <- MortSmooth_BcoefB(Ba, Bt, alphas)
    mu <- exp(OFF0 + eta)
    W <- mu
    z <- eta + (1/mu) * (Y - mu)
    z[which(WEI == 0)] <- 0
    WW <- WEI * W
    BWB <- MortSmooth_BWB(RTBa, RTBt, nba, nbt, WW)
    BWBpP <- BWB + P
    BWz <- MortSmooth_BcoefB(t(Ba), t(Bt), (WW * z))
    alphas0 <- solve(BWBpP, c(BWz))
    alphas.old <- alphas
    alphas <- matrix(alphas0, nrow = nba)
    dalphas <- max(abs(alphas-alphas.old)/abs(alphas))
    if(verbose) cat(it, dalphas, "\n")
    if(dalphas<=10^-6 & it>=4) break
  }
  ## fitted values
  ALPHAS.hat <- alphas
  ETA.hat <- eta
  Y.hat <- exp(OFF0 + ETA.hat)
  ## fitted derivatives over ages
  ETA.hat1a <- MortSmooth_BcoefB(Ca, Bt, ALPHAS.hat)
  ## fitted derivatives over years
  ETA.hat1t <- MortSmooth_BcoefB(Ba, Ct, ALPHAS.hat)
  ## diagnostics
  H <- solve(BWBpP, BWB)
  h <- diag(H)
  ed <- sum(h)
  y0 <- y
  y0[y==0] <- 10^-8
  dev <- 2 * sum(wei * (y * log(y0/mu)))
  bic <- dev + log(sum(wei))*ed
  ## deviance residuals
  res0 <- sign(Y - Y.hat)
  res1 <- log(Y/Y.hat)
  res2 <- Y - Y.hat
  res <- res0 * sqrt(2 * (Y * res1 - res2))
  res[which(is.nan(res))] <- NA
  res <- matrix(res,m,n)
  
  ## return objects
  out <- list(
    ## original data
    Y=Y, E=E, lambdas=lambdas,
    ## diagnostic
    dev=dev, ed=ed, bic=bic, res=res,
    ## fitted values
    ALPHAS=ALPHAS.hat, ETA=ETA.hat, MU=Y.hat,
    ETA1a=ETA.hat1a, ETA1t=ETA.hat1t, 
    ## basis and associated objects
    BCa=ifelse(infant, BCa0, BCa), BCt=BCt,
    ## convergence objects
    dalphas=dalphas, it=it
  )
  return(out)
}


## function to extract deltas from a PSinfant object
## for a given confidence level, default 95% and 50% 
## over age and year as reccommened/used in the manuscript
deltasFUN <- function(object, levels=c(95,50)){
  ETA1a <- object$ETA1a
  ETA1t <- object$ETA1t
  ## compute levels and deltas over ages
  p.a.up <- (100 - (100-levels[1])/2)/100
  p.a.low <- ((100-levels[1])/2)/100
  delta.a.up <- apply(ETA1a, 1, quantile,
                      probs=p.a.up)
  delta.a.low <- apply(ETA1a, 1, quantile,
                       probs=p.a.low)
  ## compute levels and deltas over years
  p.t.up <- (100 - (100-levels[2])/2)/100
  p.t.low <- ((100-levels[2])/2)/100
  delta.t.up <- apply(ETA1t, 1, quantile,
                      probs=p.t.up)
  delta.t.low <- apply(ETA1t, 1, quantile,
                       probs=p.t.low)
  ## return objects
  out <- list(delta.a.up=delta.a.up,
              delta.a.low=delta.a.low,
              delta.t.up=delta.t.up,
              delta.t.low=delta.t.low)
  return(out)
}


## function for estimating CP-splines
## for a given set of smoothing parameters lambdas
## commonly taken from the estimated object, PSinfant()
CPSfunction <- function(Y, E, lambdas, WEI,
                        kappas=c(10^4, 10^4),
                        deltas, S,
                        infant=TRUE, verbose=FALSE,
                        st.alphas=NULL){
  ## dimensions
  m <- nrow(Y)
  n <- ncol(Y)
  a <- 1:m
  t <- 1:n
  ## w/o age 0
  a0 <- a[-1]
  m0 <- m-1
  ## original offset
  OFF <- log(E)
  ## B-splines basis
  ## with infant-specilized coeff
  if(infant){
    ## over ages w/o age 0
    a0min <- min(a0)
    a0max <- max(a0)
    nda0 <- floor(m0/5)
    dega <- 3
    BCa0 <- BsplineGrad(a0, a0min, a0max, nda0, dega)
    Ba0 <- BCa0$B
    nba0 <- ncol(Ba0)
    ## adding infant-specific basis
    Ba <- cbind(0, Ba0)
    Ba <- rbind(c(1, rep(0,nba0)), Ba)
    nba <- ncol(Ba)
    ## basis for the derivatives
    ## over ages w/o age 0
    Ca0 <- BCa0$C
    ## adding infant-specific basis
    Ca <- cbind(0, Ca0)
    Ca <- rbind(c(-1,Ba[2,1:dega+1],rep(0,nba0-dega)),
                Ca)
  }else{
    amin <- min(a)
    amax <- max(a)
    nda <- floor(m/5)
    dega <- 3
    BCa <- BsplineGrad(a, amin, amax, nda, dega)
    Ba <- BCa$B
    nba <- ncol(Ba)
    ## basis for the derivatives
    Ca <- BCa$C
  }
  ## over years
  tmin <- min(t)
  tmax <- max(t)
  ndt <- floor(n/5)
  degt <- 3
  BCt <- BsplineGrad(t, tmin, tmax, ndt, degt)
  Bt <- BCt$B
  nbt <- ncol(Bt)
  ## basis for the derivatives
  Ct <- BCt$C
  ## weights for exposures=0
  WEI[E==0] <- 0
  
  ## tensor product of B-splines for the GLAM
  Ba1 <- kronecker(matrix(1, ncol=nba, nrow=1), Ba)
  Ba2 <- kronecker(Ba, matrix(1, ncol=nba, nrow=1))
  RTBa <- Ba1 * Ba2
  Bt1 <- kronecker(matrix(1, ncol=nbt, nrow=1), Bt)
  Bt2 <- kronecker(Bt, matrix(1, ncol=nbt, nrow=1))
  RTBt <- Bt1 * Bt2
  ## tensor product for the derivatives
  Ca1 <- kronecker(matrix(1, ncol=nba, nrow=1), Ca)
  Ca2 <- kronecker(Ca, matrix(1, ncol=nba, nrow=1))
  RTCa <- Ca1 * Ca2
  Ct1 <- kronecker(matrix(1, ncol=nbt, nrow=1), Ct)
  Ct2 <- kronecker(Ct, matrix(1, ncol=nbt, nrow=1))
  RTCt <- Ct1 * Ct2
  
  
  ## penalty terms for smoothing
  ## over ages
  Da <- diff(diag(nba), diff=2)
  ## no penalization over age for age 0,
  ## with infant=TRUE
  if(infant){
    Da[1,1] <- 0
  }
  tDDa <- t(Da) %*% Da
  ## over years
  Dt <- diff(diag(nbt), diff=2)
  tDDt <- t(Dt) %*% Dt
  ## kronecker product of difference matrices
  Pa <- kronecker(diag(nbt), tDDa)
  Pt <- kronecker(tDDt, diag(nba))
  ## smoothing parameters
  lambda.a <- lambdas[1]
  lambda.t <- lambdas[2]
  P <- lambda.a * Pa + lambda.t * Pt
  
  ## penalty terms for constraints
  ## extract deltas
  delta.a.up <- deltas$delta.a.up
  delta.a.low <- deltas$delta.a.low
  delta.t.up <- deltas$delta.t.up
  delta.t.low <- deltas$delta.t.low
  ## construct G matrices
  ones12 <- matrix(1, n, 1)
  g.a.up <- kronecker(ones12, delta.a.up)
  G.a.up <- matrix(g.a.up, m, n)
  g.a.low <- kronecker(ones12, delta.a.low)
  G.a.low <- matrix(g.a.low, m, n)
  g.t.up <- kronecker(ones12, delta.t.up)
  G.t.up <- matrix(g.t.up, m, n)
  g.t.low <- kronecker(ones12, delta.t.low)
  G.t.low <- matrix(g.t.low, m, n)
  
  ## data in vector for the starting values
  y <- c(Y)
  e <- c(E)
  wei <- c(WEI)
  off <- log(e)
  off0 <- off
  off0[wei==0] <- 100
  ## "other" offset in matrices
  OFF0 <- matrix(off0, m, n)
  if(is.null(st.alphas)){
    ## starting coeff
    aa <- rep(a, n)
    tt <- rep(t, each = m)
    fit0 <- glm(round(y) ~ aa + tt + offset(off0), 
                family = poisson, weights = wei)
    etaGLM <- matrix(log(fit0$fitted) - c(off0), m, n)
    eta0 <- log((Y + 1)) - OFF0
    eta0[WEI == 0] <- etaGLM[WEI == 0]
    BBa <- solve(t(Ba) %*% Ba + 1e-06*diag(nba), t(Ba))
    BBt <- solve(t(Bt) %*% Bt + 1e-06*diag(nbt), t(Bt))
    alphas <- MortSmooth_BcoefB(BBa, BBt, eta0)
  }else{
    alphas=st.alphas
  }
  ## Poisson iteration with asymmetric penalty in a GLAM framework
  for(it in 1:100){
    ## over ages
    CaBt.alphas <- MortSmooth_BcoefB(Ca, Bt, alphas)
    ## up
    v.a.up <- CaBt.alphas > G.a.up
    v.a.up <- v.a.up * S
    P.a.up <- MortSmooth_BWB(RTCa, RTBt,
                             nba, nbt, v.a.up)
    P.a.up <- kappas[1] * P.a.up
    p.a.up <- MortSmooth_BcoefB(t(Ca), t(Bt),
                                v.a.up*G.a.up)
    p.a.up <- kappas[1] * p.a.up
    ## low
    v.a.low <- CaBt.alphas < G.a.low
    v.a.low <- v.a.low * S
    P.a.low <- MortSmooth_BWB(RTCa, RTBt,
                              nba, nbt, v.a.low)
    P.a.low <- kappas[1] * P.a.low
    p.a.low <- MortSmooth_BcoefB(t(Ca), t(Bt),
                                 v.a.low*G.a.low)
    p.a.low <- kappas[1] * p.a.low
    
    P.a <- P.a.up + P.a.low
    p.a <- p.a.up + p.a.low
    ## over years
    BaCt.alphas <- MortSmooth_BcoefB(Ba, Ct, alphas)
    ## up
    v.t.up <- BaCt.alphas > G.t.up
    v.t.up <- v.t.up * S  
    P.t.up <- MortSmooth_BWB(RTBa, RTCt, nba, nbt,
                             v.t.up)
    P.t.up <- kappas[2] * P.t.up
    p.t.up <- MortSmooth_BcoefB(t(Ba), t(Ct),
                                v.t.up*G.t.up)
    p.t.up <- kappas[2] * p.t.up
    ## low
    v.t.low <- BaCt.alphas < G.t.low
    v.t.low <- v.t.low * S  
    P.t.low <- MortSmooth_BWB(RTBa, RTCt, nba, nbt,
                              v.t.low)
    P.t.low <- kappas[2] * P.t.low
    p.t.low <- MortSmooth_BcoefB(t(Ba), t(Ct),
                                 v.t.low*G.t.low)
    p.t.low <- kappas[2] * p.t.low
    
    P.t <- P.t.up + P.t.low
    p.t <- p.t.up + p.t.low
    
    eta <- MortSmooth_BcoefB(Ba, Bt, alphas)
    
    # eta[eta > 20]  <- 20
    # eta[eta < -20] <- -20
    
    mu <- exp(OFF0 + eta)
    
    # mu[mu > 1e6] = 1e6
    
    # if (anyNA(mu) || any(!is.finite(mu))) {
    #   stop("Still got non‑finite mu even after fixes—check upstream!")
    # }
    
    W <- mu
    z <- eta + (1/mu) * (Y - mu)
    z[which(WEI == 0)] <- 0
    WW <- WEI * W
    
   # # Check eta
   #  if (anyNA(eta) || any(!is.finite(eta))) {
   #    stop("Found non‑finite eta:\n", summary(as.vector(eta)))
   #  }
   # 
   #  # Check mu and WW
   #  mu <- exp(OFF0 + eta)
   #  if (anyNA(mu) || any(!is.finite(mu))) {
   #    stop("Found non‑finite mu:\n", summary(as.vector(mu)))
   #  }
   #  WW <- WEI * mu
   #  if (anyNA(WW) || any(!is.finite(WW))) {
   #    stop("Found non‑finite WW:\n", summary(as.vector(WW)))
   #  }
    
    BWB <- MortSmooth_BWB(RTBa, RTBt, nba, nbt, WW)
    
    # message(" lambdas = ", paste(signif(lambdas,3), collapse=", "))
    # message(" kappas  = ", paste(signif(kappas,3), collapse=", "))
    # message(" anyNA(BWB)? ", anyNA(BWB),
    #         " any(!is.finite(BWB))? ", any(!is.finite(BWB)))
    # message(" anyNA(P)?   ", anyNA(P),
    #         " any(!is.finite(P))?   ", any(!is.finite(P)))
    
    BWBpP <- BWB + P + P.a + P.t
    
    ## Added adaptive ridge
    
    # make_system_pd <- function(M, rel_eps = 1e-6, max_iter = 10, 
    #                            cond_thresh = 1e12) {
    #   eps0 <- rel_eps * max(abs(diag(M)))
    #   M2   <- M
    #   for (i in seq_len(max_iter)) {
    #     cond <- tryCatch(kappa(M2), error = function(e) Inf)
    #     if (cond < cond_thresh) break
    #     eps0 <- eps0 * 10
    #     M2   <- M + eps0 * diag(nrow(M))
    #   }
    #   list(M_pd = M2, eps = eps0, cond = cond)
    # }
    # 
    # # in your loop
    # pd <- make_system_pd(BWBpP)
    # BWBpP <- pd$M_pd
    
    BWz <- MortSmooth_BcoefB(t(Ba), t(Bt), (WW * z))
    BWzpP <- BWz + p.a + p.t
    
    safe_solve <- function(mat, vec, rel_eps = 1e-8) {
      # 1) check finiteness
      if (anyNA(mat) || any(!is.finite(mat))) {
        warning("safe_solve: BWBpP contains non‑finite entries → returning NA")
        return(rep(NA_real_, length(vec)))
      }
      if (anyNA(vec) || any(!is.finite(vec))) {
        warning("safe_solve: BWzpP contains non‑finite entries → returning NA")
        return(rep(NA_real_, length(vec)))
      }

      # 2) add a ridge scaled to the problem size
      ridge <- rel_eps * max(abs(diag(mat)))
      mat   <- mat + ridge * diag(nrow(mat))

      # 3) solve & check output
      sol <- tryCatch(
        solve(mat, vec),
        error = function(e) e
      )
      if (inherits(sol, "error")) {
        warning("safe_solve: solve() failed: ", sol$message)
        return(rep(NA_real_, length(vec)))
      }
      if (anyNA(sol) || any(!is.finite(sol))) {
        warning("safe_solve: solution contains non‑finite values → returning NA")
        return(rep(NA_real_, length(vec)))
      }
      sol
    }
    
    
    # alphas0 <- solve(BWBpP, c(BWzpP))
    alphas0 <- safe_solve(BWBpP, c(BWzpP), rel_eps=1e-6)
    
    alphas.old <- alphas
    alphas <- matrix(alphas0, nrow = nba)
    dalphas <- max(abs(alphas-alphas.old)/abs(alphas))
    if(verbose) cat(it, dalphas, "\n")
    if(dalphas<=10^-4 & it>=4) break   
  }
  ## fitted values
  ALPHAS.hat <- alphas
  ETA.hat <- eta
  Y.hat <- exp(OFF0 + ETA.hat)
  ## fitted derivatives over ages
  ETA.hat1a <- MortSmooth_BcoefB(Ca, Bt, ALPHAS.hat)
  ## fitted derivatives over years
  ETA.hat1t <- MortSmooth_BcoefB(Ba, Ct, ALPHAS.hat)
  ## diagnostics
  H <- solve(BWBpP, BWB)
  h <- diag(H)
  ed <- sum(h)
  y0 <- y
  y0[y==0] <- 10^-8
  dev <- 2 * sum(wei * (y * log(y0/mu)))
  bic <- dev + log(sum(wei))*ed
  
  ## return objects
  out <- list(
    ## original data
    Y=Y, E=E, lambdas=lambdas,
    kappas=kappas, deltas=deltas, S=S,
    ## diagnostic
    dev=dev, ed=ed, bic=bic,
    ## fitted values
    ALPHAS=ALPHAS.hat, ETA=ETA.hat, MU=Y.hat,
    ETA1a=ETA.hat1a, ETA1t=ETA.hat1t, 
    ## basis and associated objects
    BCa=ifelse(infant, BCa0, BCa), BCt=BCt,
    ## convergence
    dalphas=dalphas, it=it
  )
  return(out)
}


mort2Dsmooth.fun <- function(x, infant=TRUE) {
  
  frame <- arrange(x, Age)
  
  ## ages
  ages <- as.integer(unique(frame$Age))
  m <- length(ages)
  ## observed years
  years <- as.integer(unique(frame$Year))
  n <- length(years)
  n1 = length(years[years<2020])
  
  ## selecting ages and years
  E0 <- select(frame, Age, Year, Ex)
  E0 <- pivot_wider(E0, id_cols = Age, names_from = Year, values_from = Ex)
  
  ## place in a mXn1 matrix
  E <- as.matrix(E0[-1], label = TRUE)
  rownames(E) <- ages
  
  mask <- is.na(E) & (col(E) > 1)
  
  E[mask] <- E[cbind(
    row(E)[mask],       # row index of each TRUE
    col(E)[mask] - 1    # column index minus one
  )]
  
  E1 = E[,0:n1]
  
  print(E1[1,1])
  
  ## selecting ages and years
  Y0 <- select(frame, Age, Year, Dx)
  Y0 <- pivot_wider(Y0, id_cols = Age, names_from = Year, values_from = Dx)
  
  ## place in a mXn1 matrix
  Y <- as.matrix(Y0[-1], label = TRUE)
  rownames(Y) <- ages
  
  mask <- is.na(Y) & (col(Y) > 1)
  
  Y[mask] <- Y[cbind(
    row(Y)[mask],       # row index of each TRUE
    col(Y)[mask] - 1    # column index minus one
  )]
  
  Y1 = Y[,0:n1]
  
  WEI1 <- matrix(1, m, n1)
  WEI <- cbind(WEI1, matrix(0, m, n-n1))
  
  
  # Y1[is.na(Y1)] <- 0.1 ##NAs are not allowed
  # E1[is.na(E1)] <- 0.1
  ## observed log-rates
  ETA1 <- log(Y1/E1)
  
  ## function for simply extract BIC, given deaths, exposures and weights
  BICinf <- function(par){
    FITinf <- PSinfant(Y=Y1, E=E1, lambdas=par, WEI=WEI1, infant = infant)
    FITinf$bic
  }
  
  # BICinf <- function(par) {
  #   tryCatch({
  #     FIT <- PSinfant(Y=Y1, E=E1, lambdas=par, WEI=WEI1)
  #     FIT$bic
  #   }, error = function(e) {
  #     warning(
  #       sprintf("PSinfant failed at λ=(%.3g,%.3g): %s",
  #               par[1], par[2], conditionMessage(e))
  #     )
  #     Inf   # tell cleversearch this is a terrible fit
  #   })
  # }
  
  ## optimizing lambdas using greedy grid search
  OPTinf <- cleversearch(BICinf, lower=c(-4, 1), upper=c(0, 5),
                         ngrid=5, logscale=TRUE, verbose=FALSE)
  
  ## estimating mortality with optimal lambdas
  FITinf <- PSinfant(Y=Y1, E=E1, lambdas=OPTinf$par, WEI=WEI1, infant = infant)
  ## extract estimated linear predictor, log-mortality
  ETA1.hatI <- FITinf$ETA
  
  ## extract deltas from PSinfant() fitted object 
  deltas <- deltasFUN(FITinf)
  
  ## where to apply the constraints
  S <- matrix(1, m, n)
  S[,1:n1] <- 0
  if(length(ages)<10){S[1,] <- 1}
  
  # ## modelling with CP-splines
  FITcon <- CPSfunction(Y=Y, E=E, lambdas=OPTinf$par, WEI=WEI, infant = infant, deltas=deltas, S=S, verbose=TRUE)
  # ## estimated and forecast linear predictor, log-mortality
  # ETA.hatC <- FITcon$ETA
  
  # ETA.hatC <- tryCatch(
  #   {
  #     FITcon <- CPSfunction(Y=Y1, E=E1, lambdas=OPTinf$par,
  #                           WEI=WEI1, deltas=deltas, S=S,
  #                           verbose=TRUE)
  #     FITcon$ETA
  #   },
  #   error = function(e) {
  #     warning("CPSfunction failed for this cohort: ", conditionMessage(e))
  #     # return an m × n matrix of NAs
  #     matrix(NA_real_, nrow = m, ncol = n)
  #   }
  # )
  
  ## extract deviance residuals
  res <- FITinf$res
  ## replace res == NA (where Y1==0) with random residuals
  ## anyway they will be taken randomly in the next step
  whi0 <- which(is.na(res))
  res1 <- c(res)[-whi0]
  res[whi0] <- sample(res1, size=length(whi0), replace=T)
  ## small correction of actual deaths when equal to zero
  ## to avoid computational issues when log is taken and we numerically find the root of the function
  Y11 <- Y1
  Y11[whi0] <- 1e-3
  
  ## bootstrapping procedure 
  
  ## number of instances
  n.sim <- 100
  ## empty array for bootstrap deaths
  Y1s <- array(0, dim=c(m,n1,n.sim))
  ## function to invert the bootstrap-deaths
  InvDev <- function(X, Z, C){X - Z * log(X) - C}
  
  ## resample deaths
  for(k in 1:n.sim){
    ## resample (with replacement) the residuals
    resS <- matrix(sample(x=res, size=m*n1, replace = TRUE), m, n1)  
    ## computing C =  res^2/2 + Dth - Dth*ln(Dth)
    C <- 1/2*(resS^2) + Y11 - Y11*log(Y11)
    ## recomputing the deaths
    LOW <- Y11
    LOW[which(resS>0)] <- 1e-6
    UP <- Y11
    UP[which(resS<0)] <- 1e6
    UP <- ifelse(LOW>=UP, UP+10^-8, UP)
    Y1s.i <- matrix(NA,m,n1)
    
    # for(i in 1:m){
    #   for(j in 1:n1){
    #     
    #     fl <- InvDev(LOW[i,j], Z=Y11[i,j], C=C[i,j])
    #     fu <- InvDev(UP[i,j],  Z=Y11[i,j], C=C[i,j])
    #     
    #     # Check for a proper sign‐change
    #     if (is.na(fl) || is.na(fu) || fl * fu > 0) {
    #       warning(sprintf("No bracket at [%d,%d]: f_low=%.3g, f_up=%.3g",
    #                       i, j, fl, fu))
    #       Y1s.i[i,j] <- NA_real_
    #     } else {
    #       Y1s.i[i,j] <- uniroot(InvDev,
    #                             lower = LOW[i,j],
    #                             upper = UP[i,j],
    #                             tol   = 1e-5,
    #                             Z     = Y11[i,j],
    #                             C     = C[i,j])$root
    #     }
    #     
    #     
    #     # Y1s.i[i,j]  <- uniroot(InvDev, lower=LOW[i,j], upper=UP[i,j],
    #     #                        tol = 0.00001, Z=Y11[i,j], C=C[i,j])$root
    #   }
    # }
    
    eps      <- 1e-6
    big      <- 1e8
    max_tries <- 10
    
    for (i in 1:m) {
      for (j in 1:n1) {
        root_found <- FALSE
        tries      <- 0
        
        while (!root_found && tries < max_tries) {
          tries <- tries + 1
          
          # draw one residual for this cell
          r_ij <- sample(res, 1, replace=TRUE)
          
          # recompute C, LOW and UP for this single cell
          C_ij   <- 0.5 * r_ij^2 + Y11[i,j] - Y11[i,j] * log(Y11[i,j])
          if (r_ij > 0) {
            low_ij <- eps
            up_ij  <- Y11[i,j]
          } else {
            low_ij <- Y11[i,j]
            up_ij  <- big
          }
          # make sure low < up
          if (low_ij >= up_ij) up_ij <- low_ij + 1e-8
          
          fl <- InvDev(low_ij, Z = Y11[i,j], C = C_ij)
          fu <- InvDev(up_ij,  Z = Y11[i,j], C = C_ij)
          
          if (is.finite(fl) && is.finite(fu) && fl * fu < 0) {
            # bracket OK, solve and exit loop
            Y1s.i[i,j] <- uniroot(InvDev,
                                  lower = low_ij,
                                  upper = up_ij,
                                  tol   = 1e-5,
                                  Z     = Y11[i,j],
                                  C     = C_ij)$root
            root_found <- TRUE
          }
          # else: loop and redraw r_ij
        }
        
        if (!root_found) {
          warning(sprintf(
            "Cell [%d,%d]: bracket failed after %d tries → using original Y11=%.3g",
            i, j, max_tries, Y11[i,j]
          ))
          # fallback: no resampling, just keep the original
          Y1s.i[i,j] <- Y11[i,j]
        }
      }
    }
    
    Y1s[,,k] <- Y1s.i
    cat(k, "\n")
  }
  
  ## fit each sampled/bootstrap matrix of death counts
  ## and from the derivatives, compute associated confidence intervals
  DELTAs <- list()
  for(k in 1:n.sim){
    Y1.k <- Y1s[,,k]
    FITinf.k <- PSinfant(Y=Y1.k, E=E1, lambdas=OPTinf$par, WEI=WEI1, infant = infant)
    deltas.k <- deltasFUN(FITinf.k)
    DELTAs[[k]] <- deltas.k
    cat(k, "\n")
  }
  
  ## fit and forecast each sampled/boostrap matrix of death counts
  ## with CP-splines
  ETA.hatCs <- array(0, dim=c(m, n, n.sim))
  
  for(k in 1:n.sim){
    Y.k <- Y
    Y.k[,1:n1] <- Y1s[,,k]
    FITcon.k <- CPSfunction(Y=Y.k, E=E, lambdas=OPTinf$par,
                            WEI=WEI, infant = infant, deltas=DELTAs[[k]], S=S,
                            st.alphas=FITcon$ALPHAS)
    ETA.hatCs[,,k] <- FITcon.k$ETA
    cat(k, FITcon.k$it, FITcon.k$dalphas, "\n")
  }
  
  dimnames(ETA.hatCs) <- list(
    Age       = ages,
    Year      = years,
    Iteration = seq_len(dim(ETA.hatCs)[3])
  )
  
  out <- as.data.frame.table(
    ETA.hatCs,
    responseName = "smoothed_logmx"      # name of the value‑column
  )
  
  out$Age       <- as.integer(as.character(out$Age))
  out$Year      <- as.integer(as.character(out$Year))
  out$Iteration <- as.integer(as.character(out$Iteration))
  
  out <- left_join(out, frame, by = c("Age", "Year"))
  
  return(out)
}


ESP_sdr.fun = function(logmx, AgePattern) {
  
  mx = exp(logmx)
  
  AgePattern = first(AgePattern)
  
  if(AgePattern=="21"){
    
    ESP <- read.csv(file = "data/european_standard_population.csv") %>%
      mutate(Age = as.numeric(str_sub(AgeGroup, start = 1, end = 2)),
             ESPw = EuropeanStandardPopulation/100000)
    
    
  }else if(AgePattern=="6"){
    
    ESP <- read.csv(file = "data/european_standard_population.csv") %>%
      mutate(AgeGroup = as.numeric(str_sub(AgeGroup, start = 1, end = 2)),
             Age = case_when(AgeGroup == 0 ~ 0,
                             AgeGroup %in% c(1, 5, 10, 15, 20) ~ 1,
                             AgeGroup %in% c(25, 30, 35, 40, 45) ~ 25,
                             AgeGroup %in% c(50, 55, 60, 65) ~ 50,
                             AgeGroup %in% c(70, 75, 80) ~ 70,
                             AgeGroup %in% c(85, 90, 95) ~ 85)) %>%
      group_by(Age) %>%
      summarise(ESPw = sum(EuropeanStandardPopulation)/100000)
    
  }else if(AgePattern=="POL5"){
    
    ESP <- read.csv(file = "data/european_standard_population.csv") %>%
      mutate(AgeGroup = as.numeric(str_sub(AgeGroup, start = 1, end = 2)),
             AgeGroup = case_when(AgeGroup %in% c(0, 1, 5, 10, 15, 20) ~ 0,
                                  AgeGroup %in% c(25, 30, 35, 40, 45) ~ 25,
                                  AgeGroup %in% c(50, 55, 60, 65) ~ 50,
                                  AgeGroup %in% c(70, 75, 80) ~ 70,
                                  AgeGroup %in% c(85, 90, 95) ~ 85)) %>%
      group_by(AgeGroup) %>%
      summarise(ESPw = sum(EuropeanStandardPopulation)/100000)
    
  }else if(AgePattern=="ROU6"){
    
    ESP <- read.csv(file = "data/european_standard_population.csv") %>%
      mutate(AgeGroup = as.numeric(str_sub(AgeGroup, start = 1, end = 2)),
             AgeGroup = case_when(AgeGroup %in% c(0,1) ~ 0,
                                  AgeGroup %in% c(5, 10, 15, 20) ~ 5,
                                  AgeGroup %in% c(25, 30, 35, 40, 45) ~ 25,
                                  AgeGroup %in% c(50, 55, 60, 65) ~ 50,
                                  AgeGroup %in% c(70, 75, 80) ~ 70,
                                  AgeGroup %in% c(85, 90, 95) ~ 85)) %>%
      group_by(AgeGroup) %>%
      summarise(ESPw = sum(EuropeanStandardPopulation)/100000)
    
  }else if(AgePattern=="16"){
  
  ESP <- read.csv(file = "data/european_standard_population.csv") %>%
    mutate(AgeGroup = as.numeric(str_sub(AgeGroup, start = 1, end = 2)),
           AgeGroup = case_when(AgeGroup %in% c(0, 1, 5) ~ 0,
                                AgeGroup %in% c(10, 15) ~ 10,
                                AgeGroup == 20 ~ 20,
                                AgeGroup == 25 ~ 25,
                                AgeGroup == 30 ~ 30,
                                AgeGroup == 35 ~ 35,
                                AgeGroup == 40 ~ 40,
                                AgeGroup == 45 ~ 45,
                                AgeGroup == 50 ~ 50,
                                AgeGroup == 55 ~ 55,
                                AgeGroup == 60 ~ 60,
                                AgeGroup == 65 ~ 65,
                                AgeGroup == 70 ~ 70,
                                AgeGroup == 75 ~ 75,
                                AgeGroup == 80 ~ 80,
                                AgeGroup %in% c(85, 90, 95) ~ 85)) %>%
    group_by(AgeGroup) %>%
    summarise(ESPw = sum(EuropeanStandardPopulation)/100000)
  
  }else if(AgePattern=="POL14"){
    
    
    ESP <- read.csv(file = "data/european_standard_population.csv") %>%
      mutate(AgeGroup = as.numeric(str_sub(AgeGroup, start = 1, end = 2)),
             AgeGroup = ifelse(AgeGroup >= 65, 65, AgeGroup),
             AgeGroup = ifelse(AgeGroup %in% 0:1, 0, AgeGroup)) %>%
      group_by(AgeGroup) %>%
      summarise(ESPw = sum(EuropeanStandardPopulation)/100000)
    
  }else if(AgePattern=="18"){
    
    ESP <- read.csv(file = "data/european_standard_population.csv") %>%
      mutate(AgeGroup = as.numeric(str_sub(AgeGroup, start = 1, end = 2)),
             AgeGroup = ifelse(AgeGroup >= 85, 85, AgeGroup),
             AgeGroup = ifelse(AgeGroup %in% 0:1, 0, AgeGroup)) %>%
      group_by(AgeGroup) %>%
      summarise(ESPw = sum(EuropeanStandardPopulation)/100000)
    
  }
  
  ESPw = pull(ESP, ESPw)

  sdr = 1000 * weighted.mean(mx, ESPw, na.rm = T)
  
  return(sdr)
}
