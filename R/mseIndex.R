# mseIndex.R - DESC
# mse/R/mseIndex.R

# Copyright European Union, 2017
# Author: Iago Mosqueira (EC JRC) <iago.mosqueira@ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.


# mseIndex {{{

#' @title mseIndex
#' @description FUNCTION_DESCRIPTION
#' @param omp The population OM, an FLStock object extended to the projection years
#' @param sr The stock-recruit relationship, as an FLSR, predictModel or list ibjcts with 'model' and 'params'
#' @param cpue The observed coue series to be extended and used by the MP
#' @param cpuesel PARAM_DESCRIPTION
#' @param years Vector of years on which MP is to be evaluated
#' @param verbose Show output on screen?, Default: FALSE
#' @param hcr PARAM_DESCRIPTION, Default: ~tac * (1 + lambda * slope)
#' @param hcrparams PARAM_DESCRIPTION, Default: FLPar(lambda = 1.25, ny = 5, dltac = 0.15, dhtac = 0.15)
#' @param dlag PARAM_DESCRIPTION, Default: 1
#' @param mlag PARAM_DESCRIPTION, Default: 1
#' @param oemparams PARAM_DESCRIPTION, Default: FLPar(sd = 0, b = 0)
#' @param imparams PARAM_DESCRIPTION
#' @param tune PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' data(cod)
#' @seealso 
#'  \code{\link[utils]{txtProgressBar}},\code{\link[utils]{setTxtProgressBar}}
#' @rdname mseIndex


mseIndex <- function(
  # OM: FLStock + SR + RPs + cpue
  omp, sr, cpue, cpuesel,
  # years
  years, verbose=FALSE,
  # hcr
  hcr=~tac * (1 + lambda * slope),
  # hcrparams
  hcrparams=FLPar(lambda=1.25, ny=5, dltac=0.15, dhtac=0.15),
  # lags
  dlag=1, mlag=1, 
  # oem, imp
  oemparams=FLPar(sd=0, b=0), imparams, tune=FALSE) {

  # VARIABLES
  freq <- years[2] - years[1]

  # MESSAGES
  if(verbose)
    pb <- utils::txtProgressBar(min = years[1], max = years[length(years)],
      initial = 1, style=3)
  
  # TAC
  tac <- catch(omp)[, ac(seq(years[1] - dlag, years[length(years)] + freq))]

  # LOOP
  for (y in years) {

    # CATCH
    # TODO + E
    stk <- window(omp, end=c(y - dlag))

    # oem w/ selectivity[, y - dlag] in weight
    obs <- quantSums(cpue(stk[,ac(seq(y - dlag - freq, y - dlag))],
      sel.pattern=cpuesel, mass=TRUE))
    
    # EXTEND cpue from delta(obs)
    cpue[, ac(seq(y - dlag - freq + 1, y - dlag))] <- 
    # TODO ADD hyperstability
      # if obs[t] < obs[t-1], then obs[t-1] + E
    # TODO ADD bias
      # obs[t] * bias
      cpue[, ac(y - dlag - freq)] %*% obs[,-1] / obs[, -dim(obs)[2]] %*%
      # E: LN(0, 0.3) + b
      # TODO b from history, sd from OM or actual value
      rlnoise(1, FLQuant(0, dimnames=dimnames(obs[,-1])[-6]),
        sd=c(oemparams$sd), b=c(oemparams$b))
    # sd ~ sampling + E(stk~cpue)
    
    # INDICATOR
    dat <- data.table(as.data.frame(cpue[,
      ac(seq(y - dlag - hcrparams$ny - 1, y - dlag))], drop=FALSE))
    dat$iter <- as.numeric(dat$iter)

    # BUG: DT
    # setnames(dat, "data", "da")
    # CALCULATE slope
    # foo <- function(x)
    #  coef(lm(log(data)~year, data=x, na.action=na.exclude))[2]
    # slope <- dat[, list(V1=foo(x)), by = iter]$V1

    slope <- unlist(as.list(by(dat[, c("data", "year", "iter")], dat[, "iter"],
      function(x)
        coef(lm(log(data)~year, data=x, na.action=na.exclude))[2], simplify=TRUE)))

    # DECISION
    # TODO GENERALIZE based on formula (e.g. ~ssb)
    ytac <- eval(hcr[[2]], c(as(hcrparams, 'list'),
      list(tac=c(tac[,ac(y - dlag)]), slope=slope)))
    
    # CONSTRAINT in TAC change
    ptac <- c(tac[, ac(y-dlag)])
    ytac <- pmax(ptac * c((1 - hcrparams$dltac)),
      pmin(ptac * c((1 + hcrparams$dhtac)), ytac))

    # LOG tac
    tac[, ac(seq(y + mlag, length=freq))] <- rep(ytac, each=freq)

    # FWD w/IMP. ERROR + SR residuals
    # TODO ADD imp error
      # tac + N(tac, sigma)
    omp <- fwd(omp, sr=sr, residuals=sr$residuals,
      control=fwdControl(quant="catch", year=seq(y + mlag, length=freq),
      value=rep(ytac, freq)))

    # DONE
    if(verbose)
      utils::setTxtProgressBar(pb, y)
  }

  if(verbose)
    cat("\n")
  # END
  if(tune)
    return(window(omp, start=years[1] - dlag - 1, end=years[length(years)]))
  else
    return(list(om=window(omp, start=years[1] - dlag - 1, end=years[length(years)]),
      tac=tac, cpue=cpue))

} # }}}
