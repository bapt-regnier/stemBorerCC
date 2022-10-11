#' computeDT: get development times for different temperatures
#' @param nlsList list of nls fits for each life-stage
#' @param nInd number of individuals simulated
#' @param tempSeq vector of temperatures
computeDT <- function(nlsList, nInd = 5000, tempSeq = seq(0, 50, 0.01)){
  source(file = "scripts/function_getDT.R")
  ibmList <- lapply(seq_along(nlsList), function(i){
    models <- list(nlsList[[i]][[1]][[1]],
                   nlsList[[i]][[2]][[1]],
                   nlsList[[i]][[3]][[1]])
    eqName <- c(names(nlsList[[i]][[1]]),
                names(nlsList[[i]][[2]]),
                names(nlsList[[i]][[3]]))
    denList3 <- lapply(tempSeq, function(x){
      getDTx(Tx = x, eqName = eqName, models = models, nInd = nInd)
    })
    denList3 <- do.call(rbind, denList3)
    return(denList3)
  })
  return(ibmList)
}
#' getDTstats: get development times median and 95% prediction interval
#' @param ibmList simulated individual development times for each species (list)
getDTstats <- function(ibmList){
  statList <- lapply(seq_along(ibmList), function(i){
    med <- aggregate(x = ibmList[[i]]$days,
                     by = list(ibmList[[i]]$addTemp),
                     FUN = median, # median for each temperature
                     na.rm = TRUE)
    q <- aggregate(x = ibmList[[i]]$days,
                   by = list(ibmList[[i]]$addTemp),
                   FUN = quantile, # quantiles 95%
                   probs = c(0.025, 0.975),
                   na.rm = TRUE)
    return(data.frame(temp = med$Group.1,
                      med = med$x,
                      q1 = q$x[, 1],
                      q2 = q$x[, 2]))
  })
  return(statList)
}

#' plotDT: plot complete development time as a function of temperature
#' @param nlsList list of selected nls
#' @param statList list of median development time and 95% prediction interval
#' @param bioVal extracted biological parameters values of selected model fits
#' @param maxLim limit for median development time in days
plotDT <- function(nlsList, statList, bioVal, maxLim = 182){
  sp <- c("Chilo partellus",
          "Busseola fusca",
          "Ostrinia nubilalis",
          "Sesamia nonagrioides")
  colX <- c(rgb(223, 83, 107, 255/6, maxColorValue = 255),
            rgb(97, 208, 79, 255/6, maxColorValue = 255),
            rgb(34, 151, 230, 255/6, maxColorValue = 255),
            rgb(40, 226, 229, 255/6, maxColorValue = 255))
  colX2 <- c(rgb(223, 83, 107, 255, maxColorValue = 255),
             rgb(97, 208, 79, 255, maxColorValue = 255),
             rgb(34, 151, 230, 255, maxColorValue = 255),
             rgb(40, 226, 229, 255, maxColorValue = 255))
  medDevT <- statList[[1]]$med[statList[[1]]$med <= maxLim]
  addTemp <- statList[[1]]$temp[statList[[1]]$med  <= maxLim]
  maxXlim <- max(sapply(seq_along(statList), function(i){
    max(statList[[i]]$q2[statList[[i]]$med <= maxLim], na.rm = TRUE)
  }))
  pdf("./results/completeDev_allSp.pdf", width = 4, height = 4)
  par(mar = c(4, 4, 0.1, 0.1))
  plot(NULL,
       xlim = c(5, 40),
       ylim = c(0, maxXlim),
       ylab = "Development time (days)",
       xlab = "Temperature (°C)",
       mgp = c(2, 0.8, 0), bty = "l")
  grid()
  lapply(seq_along(statList), function(i){
    statList[[i]] <- statList[[i]][!is.na(statList[[i]]$med), ]
    medDevT <- statList[[i]]$med[statList[[i]]$med <= maxLim]
    addTemp <- statList[[i]]$temp[statList[[i]]$med <= maxLim]
    Tmin <- addTemp[1]
    Tmax <- addTemp[length(addTemp)]
    print(sp[i])
    print(paste0("Tmin = ", Tmin))
    print(paste0("Tmax = ", Tmax))
    q1 <- statList[[i]]$q1[statList[[i]]$med <= maxLim]
    q2 <- statList[[i]]$q2[statList[[i]]$med <= maxLim]

    lines(medDevT ~ addTemp, col = colX2[i], lwd = 2)
    polygon(
      c(min(addTemp),
        addTemp,
        max(addTemp),
        rev(addTemp), min(addTemp)),
      c(q1[addTemp == min(addTemp)],
        q1,
        q2[addTemp == max(addTemp)],
        rev(q2),
        q2[addTemp == min(addTemp)]),
      col = colX[i], lty = 0)
    abline(h = maxLim, lty = 2)
  })
  legend("bottomleft", legend = sp, col = colX2, lwd = 2, bty = "n")
  dev.off()
}

