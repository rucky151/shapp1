

#' Generate an output matrix for one simulation
#'
#' This function simulates one parameter combination across nseasons once.
#'
#' Refer \code{\link{SeasModstoch}} for relative information.
#'

#' @inheritParams SeasModstoch
#' @importFrom stats median quantile rnorm var
#' @keywords seed health
#' @export
#' @examples
#' SeasModdet()

# Columns of output matrix
# col 1 - season timestep (initial time step is season 0)
# col 2 - HP healthy plant number
# col 3 - DP diseased plant number (after roguing)
# col 4 - HS healthy seed number
# col 5 - DS diseased seed number
# col 6 - pHS proportion healthy seed
# col 7 - pDS proportion diseased seed
# col 8 - Yld end of season yield
# col 9 - YL end of season yield loss
# col 10 - DPbr (diseased plants before roguing)


SeasModdet <- function(pHSinit=0.8, Kx = 100, betax = 0.02, wxtnormm = 0.8, hx = 1, mxtnormm = 1, axtnormm = 1, rx = 0.1, zxtnormm = 0.9, gx = 4, cx = 0.9, phix = 0,
                       maY = 100, miY = 0, thetax = 0.2, Ex = 0, nseasons = 10) {

  outm <- as.data.frame(matrix(data=-999, nrow=(nseasons+1), ncol=10, dimnames = list(1:(nseasons+1),c('season','HP', 'DP','HS','DS','pHS','pDS','Yld','YL', 'DPbr'))))

  outm[1,] <- NA# row one gives initial conditions
  outm$season <- 0:nseasons

  # Initial value for state variables
  outm$pHS[1] <- pHSinit # initial proportion healthy seed (for nseasons=0)
  outm$pDS[1] <- 1 - pHSinit

  # seasons 2 and higher

  for(si in 2:(nseasons+1)) {
    # Equation A1
    outm$HP[si] <- min(max(0,Kx * outm$pHS[si-1] - betax * wxtnormm * hx * mxtnormm * ((Kx * outm$pDS[si-1]) * (Kx * outm$pHS[si-1])+ Ex*(Kx * outm$pHS[si-1]))), 100)
    # max prevents <0 plants & min prevents >K plants
    # Equation A2
    outm$DPbr[si] <- min(max(0, Kx * outm$pDS[si-1] + betax * wxtnormm * hx * mxtnormm * ((Kx * outm$pDS[si-1]) * (Kx * outm$pHS[si-1]) + Ex*(Kx * outm$pHS[si-1]))), 100)
    # Equation A3
    outm$DP[si] <- axtnormm*outm$DPbr[si]
    # Equation B1
    outm$HS[si] <- max(0,gx * (outm$HP[si] + rx * outm$DP[si]))
    # Equation B2
    outm$DS[si] <- max(0,zxtnormm * cx * gx * (1 - rx) * outm$DP[si])
    # Equation C1
    outm$Yld[si]<- ((outm$HP[si]+outm$DP[si])/(outm$HP[si]+outm$DPbr[si]))*(miY+(maY-miY)*((1-(outm$DP[si]/(outm$DP[si]+outm$HP[si])))/((1-thetax)+ thetax*(1-(outm$DP[si]/(outm$DP[si]+outm$HP[si]))))^2))
    #Equation C2
    outm$YL[si]<- maY-outm$Yld[si]
    # Equation D1
    outm$pHS[si] <- phix + (1-phix) * outm$HS[si]/(outm$HS[si] + outm$DS[si])
    # Equation D2
    outm$pDS[si] <- 1 - outm$pHS[si]
  }

  #-----------------------------------------------
  # warning message

  if ( pHSinit < 0 | pHSinit > 1){
    warning(paste('pHSinit: your input value is', pHSinit,', it must be between 0 and 1'))
  } else if (betax < 0.001 | betax > 0.2) {
    warning(paste('betax: your input value is', betax,', it must be between 0.001 and 0.2'))
  } else if (wxtnormm < 0 | wxtnormm > 1) {
    warning(paste('wxtnormm: your input value is', wxtnormm,', it must be between 0 and 1'))
  } else if (hx < 0 | hx > 1) {
    warning(paste('hx: your input value is', hx,', it must be between 0 and 1'))
  } else if (mxtnormm < 0 | mxtnormm > 1) {
    warning(paste('mxtnormm: your input value is', mxtnormm,', it must be between 0 and 1'))
  } else if (axtnormm < 0 | axtnormm > 1) {
    warning(paste('axtnormm: your input value is', axtnormm,', it must be between 0 and 1'))
  } else if (rx < 0 | rx > 1) {
    warning(paste('rx: your input value is', rx,', it must be between 0 and 1'))
  } else if (zxtnormm < 0 | zxtnormm > 1) {
    warning(paste('zxtnormm: your input value is', zxtnormm,', it must be between 0 and 1'))
  } else if (gx < 0 | gx > 20) {
    warning(paste('gx: your input value is', gx,', it must be between 0 and 20'))
  } else if (cx < 0 | cx > 1) {
    warning(paste('cx: your input value is', cx,', it must be between 0 and 1'))
  } else if (phix < 0 | phix > 1) {
    warning(paste('phix: your input value is', phix,', it must be between 0 and 1'))
  } else if (thetax < -1 | thetax > 0.55) {
    warning(paste('thetax: your input value is', thetax,', it must be between -1 and 0.55'))
  } else if (Ex < 0 | Ex > 50) {
    warning(paste('Ex: your input value is', Ex,', it must be between 0 and 50'))
  }else {
    list(outm[(nseasons+1),])
  }
  #-----------------------------------------------

}
