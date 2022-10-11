#' fitNLS: fit 11 models to database using Levenberg-Marquardt algorithm
#' @param poolData pooled datasets
fitNLS <- function(poolData){
  lapply(unique(poolData$completeSp), function(speciesX){
    lapply(unique(poolData$STAGE), function(stageX){
      dfSub <- poolData[poolData$STAGE == stageX &
                          poolData$completeSp == speciesX,
                        c("TEMP", "DEVRATE")]
      AIO <- devRateModelAll(
        eqList = devRateEqList[c(
          "analytis_77",
          "ratkowsky_83",
          "hilbertLogan_83",
          "beta_95",
          "beta_16",
          "briere1_99",
          "briere2_99",
          "kontodimas_04",
          "shi_11",
          "perf2_11",
          "regniere_12"
        )],
        eqStartVal = devRateEqStartVal[c(
          "analytis_77",
          "ratkowsky_83",
          "hilbertLogan_83",
          "beta_95",
          "beta_16",
          "briere1_99",
          "briere2_99",
          "kontodimas_04",
          "shi_11",
          "perf2_11",
          "regniere_12"
        )],
        dfData = dfSub,
        algo = "LM",
        control = list(
          maxiter = 1024,
          maxfev = 10000
        )
      )
      return(AIO)
    })
  })
}

#' selectNLS: select the fits that corresponds to statistical and biological
#' criteria
#' @param AIOlist list of NLS fits and goodness of fit tables
#' @param poolData pooled dataset
selectNLS <- function(AIOlist, poolData){
  lapply(1:4, function(i){
    lapply(1:3, function(j){
      nlsList <- AIOlist[[i]][[j]][[1]]
      statDf <- AIOlist[[i]][[j]][[2]]
      statDf <- na.omit(statDf)
      selectedAIC <- statDf[statDf$AIC - min(statDf$AIC) <= 10,]
      bioX <- lapply(seq_along(nlsList), function(j){
        if(!is.null(nlsList[[j]])){
          if(j == 1){
            return(t(summary(nlsList[[j]])$parameters[4:5, 1]))
          }
          if(j == 2 | j == 4 | j == 5 | j == 6 | j == 7 | j == 8){
            return(t(summary(nlsList[[j]])$parameters[2:3, 1]))
          }
          if(j == 3 | j == 11){
            return(t(summary(nlsList[[j]])$parameters[3:4, 1]))
          }
          if(j == 9){
            return(t(summary(nlsList[[j]])$parameters[c(3, 5), 1]))
          }
          if(j == 10){
            return(t(summary(nlsList[[j]])$parameters[c(2, 4), 1]))
          }
        }else{
          return(NULL)
        }
      })
      condCTmin <- sapply(bioX, "[[", 1)
      condCTmin[sapply(condCTmin, is.null)] <- -999
      condCTmax <- sapply(bioX, "[[", 2)
      condCTmax[sapply(condCTmax, is.null)] <- 999
      thresholdBio <- (unlist(condCTmin) > 0 & unlist(condCTmin) < 50) & 
        (unlist(condCTmax) < 50 & unlist(condCTmax) > 0) &
        (unlist(condCTmax) - unlist(condCTmin) > 0)
      nlsList[!thresholdBio] <- NULL
      selected <- selectedAIC[selectedAIC$eqName %in% names(nlsList),]
      selected <- selectedAIC[!is.na(selectedAIC$nParam),]
      minAIC <- selected[selected$AIC == min(selected$AIC, na.rm = TRUE) &
                           !is.na(selected$nParam),]
      selectedNLS <- nlsList[names(nlsList) == minAIC$eqName]
      return(selectedNLS)
    })
  })
}

#' selectStats: write statistics and biological parameters used for model
#' selection in .csv
#' @param poolData pooled datasets
selectStats <- function(poolData){
  statList1 <- lapply(unique(poolData$completeSp), function(speciesX){
    statList2 <- lapply(unique(poolData$STAGE), function(stageX){
      dfSub <- poolData[poolData$STAGE == stageX &
                          poolData$completeSp == speciesX,
                        c("TEMP", "DEVRATE")]
      AIO <- devRateModelAll(
        eqList = devRateEqList[c(
          "analytis_77",
          "ratkowsky_83",
          "hilbertLogan_83",
          "beta_95",
          "beta_16",
          "briere1_99",
          "briere2_99",
          "kontodimas_04",
          "shi_11",
          "perf2_11",
          "regniere_12"
        )],
        eqStartVal = devRateEqStartVal[c(
          "analytis_77",
          "ratkowsky_83",
          "hilbertLogan_83",
          "beta_95",
          "beta_16",
          "briere1_99",
          "briere2_99",
          "kontodimas_04",
          "shi_11",
          "perf2_11",
          "regniere_12"
        )],
        dfData = dfSub,
        algo = "LM",
        control = list(
          maxiter = 1024,
          maxfev = 10000
        )
      )
      nlsList <- AIO[[1]]
      statDf <- data.frame(species = speciesX, stage = stageX, AIO[[2]])
      bioX <- lapply(seq_along(nlsList), function(j){
        if(!is.null(nlsList[[j]])){
          if(j == 1){
            return(t(summary(nlsList[[j]])$parameters[4:5, 1]))
          }
          if(j == 2 | j == 4 | j == 5 | j == 6 | j == 7 | j == 8){
            return(t(summary(nlsList[[j]])$parameters[2:3, 1]))
          }
          if(j == 3 | j == 11){
            return(t(summary(nlsList[[j]])$parameters[3:4, 1]))
          }
          if(j == 9){
            return(t(summary(nlsList[[j]])$parameters[c(3, 5), 1]))
          }
          if(j == 10){
            return(t(summary(nlsList[[j]])$parameters[c(2, 4), 1]))
          }
        }else{
          return(NA)
        }
      })
      bioX <- do.call(rbind, bioX)
      statDf$CTmin <- bioX[, 1]
      statDf$CTmax <- bioX[, 2]
      statDf[, c(7, 9, 10)] <- round(statDf[, c(7, 9, 10)], 3)
      return(statDf[, c(1:4, 7, 9, 10)])
    })
    statList2 <- do.call(rbind, statList2)
  })
  statList1 <- do.call(rbind, statList1)
  write.csv(x = statList1,
            file = "./results/statFits_allSp.csv")
  return(statList1)
}

#' plotStudyEffect : plot experimental data extracted from the literature as 
#' a function of articles and life-stages
#' @param data dataset
#' @param path path for pdf plot
plotStudyEffect <- function(data, path){
  pdf(path)
  par(mfrow = c(2, 2), mar = c(4, 4, 1, 1))
  lapply(unique(data$completeSp), function(spX){
    print(spX)
    sub1 <- data[data$completeSp == spX, ]
    nStudy <- length(unique(sub1$REF))
    plot(NULL, xlim = c(-0.1, 3), ylim = c(0, max(sapply(sub1$DEVRATE, max))),
         main = sub1$completeSp[1],
         axes = FALSE, xlab = NA, ylab = "taux de developpement",
         bty = "l")
    j <- 0
    for(stageX in c("egg", "larva", "pupa")){
      sub2 <- sub1[sub1$STAGE == stageX, ]
      print(stageX)
      j <- j + 1
      i <- 0
      dfX <- lapply(seq_along(sub2$DEVRATE), function(k){
        return(data.frame(rT = sub2$DEVRATE[[k]],
                          ref = as.factor(sub2$REF[k])))
      })
      dfX <- do.call(rbind, dfX)
      sumX <- summary(aov(dfX$rT ~ dfX$ref))
      for(refX in unique(sub1$REF)){
        i <- i + 1
        
        if(length(unique(sub1$REF)) == 2) colX <- c(2, 3)
        if(length(unique(sub1$REF)) == 3) colX <- c(2, 3, 4)
        if(length(unique(sub1$REF)) == 4) colX <- c(2, 3, 4, 5)
        
        if(j == 2){
          legend(x = 0.9,
                 y = max(sapply(sub1$DEVRATE, max)),
                 bty = "n",
                 legend = unique(sub1$REF),
                 col = colX,
                 pch = 16, 
                 xpd = TRUE)
        }
        sub3 <- sub2[sub2$REF == refX,]
        print(refX)
        points(y = unlist(sub3$DEVRATE),
               x = rep(j-i/4, length(unlist(sub3$DEVRATE))), cex = 0.5,
               col = colX[i],
               pch = 16)
        text(x = (j-i/4), y = 0,
             labels = length(unlist(sub3$DEVRATE)))
        
        if(i == 2){
          box(bty = "l")
          axis(side = 2)
          axis(side = 1, at = j-i/4, labels = stageX, las = 3)
          if(sumX[[1]]$`Pr(>F)`[1] <= 0.05){
            text(x = j-i/4, y = max(sapply(sub1$DEVRATE, max)), labels = "*") 
          } 
        }
      }
    }
  })
  dev.off()
}

#' plotResiduals : graphical diagnostic plots to verify nls assumptions for each
#' selected nls fit
#' @param nlsListl list of selected nls fits
#' @param pooledData pooled dataset
#' @param path path for pdf plot
#' @param ... additional arguments for pdf function
plotResiduals <- function(nlsList, pooledData , path, ...){
  pdf(path, ...)
  # res vs fit
  par(mfrow = c(4, 3), mar = c(4, 4, 1, 1))
  lapply(seq_along(nlsList), function(spX){
    lapply(seq_along(nlsList[[spX]]), function(stageX){
      nlsX <- nlsList[[spX]][[stageX]][[1]]
      resX <- residuals(nlsX)
      predX <- predict(nlsList[[spX]][[stageX]][[1]])
      plot(resX ~ predX, type = "n", bty="l",
           xlab = "fitted", ylab = "residuals")
      mtext(text = paste(unique(pooledData$completeSp)[spX],
                         unique(pooledData$STAGE)[stageX], sep = "/"),
            side = 3,
            cex = 0.8, adj = 0)
      grid()
      modloess <- loess(resX ~ predX, degree = 1)
      lines(predict(modloess,
                    newdata = seq(0, 1, 0.001)) ~ seq(0, 1, 0.001), col = 2)
      abline(h = 0, lty = 2, col = "grey")
      points(resX ~ predX, pch = 16)
    })
  })
  # qq plot
  par(mfrow = c(4, 3), mar = c(4, 4, 1, 1))
  lapply(seq_along(nlsList), function(spX){
    lapply(seq_along(nlsList[[spX]]), function(stageX){
      nlsX <- nlsList[[spX]][[stageX]][[1]]
      resX <- residuals(nlsX)
      qqnorm(resX, main = NA)
      mtext(text = paste(unique(pooledData$completeSp)[spX],
                         unique(pooledData$STAGE)[stageX], sep = "/"),
            side = 3,
            cex = 0.8, adj = 0)
      qqline(resX, col = 2)
    })
  })
  dev.off()
}

#' extractNfiltered: extract the number of fits filtered through AIC comparison
#' and biological criteria
#' @param statDf list of data.frame with statistics used for model selection

extractNfiltered <- function(statDf){
  dfSpecies <- lapply(unique(statDf$species), function(speciesX){
    dfStage <- lapply(unique(statDf$stage), function(stageX){
      dfSub <- statDf[statDf$species == speciesX &
                        statDf$stage == stageX,]
      dfSub <- na.omit(dfSub)
      nDelta10 <- nrow(dfSub[dfSub$AIC - min(dfSub$AIC) >= 10,])
      nDelta2 <- nrow(dfSub[dfSub$AIC - min(dfSub$AIC) <= 2,])
      nCTmin0 <- nrow(dfSub[dfSub$CTmin < 0,])
      nCTmax50 <- nrow(dfSub[dfSub$CTmax > 50, ])
      dfBio <- dfSub[dfSub$CTmin < 0 | dfSub$CTmax > 50,]
      minAIC <- min(dfSub$AIC)
      isMinAICfiltered <- function(df, minAIC){
        if(minAIC %in% df$AIC){
          return("yes")
        }else{
          return("no")
        }
      }
      return(
        data.frame(
          species = speciesX,
          stage = stageX,
          nDelta10 = nDelta10,
          nDelta2 = nDelta2,
          nCTmin0 = nCTmin0,
          nCTmax50 = nCTmax50,
          isMinAICfiltered = isMinAICfiltered(df = dfBio, minAIC = minAIC)
        )
      )
    })
    dfStage <- do.call(rbind, dfStage)
  })
  dfAll <- do.call(rbind, dfSpecies)
  write.csv(x = dfAll, file = "./results/nFiltered.csv")
  return(dfAll)
}


#' extractParam: extract parameters values of model fits
#' @param AIOlist list of nls fits
#' @param poolData pooled dataset
extractParam <- function(AIOlist, poolData){
  params <- lapply(seq_along(AIOlist), function(i){
    xxx <- lapply(seq_along(AIOlist[[i]]), function(j){
      nlsList <- AIOlist[[i]][[j]][[1]]
      zzz <- lapply(seq_along(nlsList), function(k){
        modName <- names(nlsList)[k]
        if(!is.null(nlsList[[k]])){
          parX <- summary(nlsList[[k]])$p[, 1:2]
          return(
            data.frame(
              stage = unique(poolData$STAGE)[j],
              sp = unique(poolData$completeSp)[i],
              modName,
              parName = rownames(parX),
              parX
            )
          )
        }else{
          return(
            data.frame(
              stage = unique(poolData$STAGE)[j],
              sp = unique(poolData$completeSp)[i],
              modName,
              parName = NA,
              Estimate = NA,
              Std..Error = NA
            )
          )
        }
      })
      return(do.call(rbind, zzz))
    })
    return(do.call(rbind, xxx))
  })
  params <- do.call(rbind, params)
  rownames(params) <- NULL
  write.csv(x = params, file = "./results/paramValues.csv")
  return(params)
}

#' plotTPC: plot TPC for selected fits
#' @param poolData pooled datasets
#' @param nlsList list of selected nls fits
plotTPC <- function(poolData, nlsList){
  pdf("./results/TPCs_lifeStages.pdf", width = 20, height = 5)
  split.screen(c(1, 4))
  for(i in 1:3){
    screen(n = i)
    par(mar = c(4, 4, 2, 2), cex = 1.5, xpd = TRUE)
    plot(NULL,
         xlim = c(0, 50),
         ylim = c(0,
                  max(poolData[poolData$STAGE == unique(poolData$STAGE)[i],
                               "DEVRATE"])),
         xlab = "Temperature",
         ylab = "Development rate",
         bty = "l"
    )
    par(xpd = FALSE)
    grid()
    mtext(paste0(letters[i], ") ", unique(poolData$STAGE)[i]),
          side = 3,
          at = 0,
          line = 0.5,
          cex = 1.5)
    for(j in 1:4){
      nlsX <- nlsList[[j]][[i]]
      points(
        poolData$DEVRATE[
          poolData$completeSp == unique(poolData$completeSp)[j] &
            poolData$STAGE == unique(poolData$STAGE)[i]
        ] ~
          poolData$TEMP[
            poolData$completeSp == unique(poolData$completeSp)[j] &
              poolData$STAGE == unique(poolData$STAGE)[i]
          ],
        col = j+1,
        pch = 16,
        cex = 0.7
      )
      seqTemp <- seq(from = 0, to = 50, by = 0.1)
      pred <- predict(nlsX[[1]],
                      newdata = list(T = seqTemp))
      pred[pred < 0] <- 0
      pred[
        seqTemp < coef(nlsX[[1]])[
          grep(pattern = "Tmin",
               x = names(coef(nlsX[[1]])))
        ]
      ] <- 0
      pred[
        seqTemp > coef(nlsX[[1]])[
          grep(pattern = "Tmax",
               x = names(coef(nlsX[[1]])))
        ]
      ] <- 0
      lines(pred ~ seqTemp, col = j+1, lwd = 2)
    }
  }
  screen(n = 4)
  par(mar = c(0, 0, 4, 0), xpd = TRUE)
  plot(
    NULL,
    xlim = c(0, 1),
    ylim = c(0, 1),
    type = "n",
    axes = FALSE,
    ann = FALSE
  )
  legend(
    "topleft",
    pch = 16,
    lty = 1,
    col = c(2:5),
    legend = unique(poolData$completeSp),
    bty = "n",
    cex = 1.5,
    text.font = 3
  )
  dev.off()
}

#' extractBioVal: extract biological parameters values of selected model fits
#' @param nlsList list of nls fits
#' @param poolData pooled dataset
extractBioVal <- function(nlsList, poolData){
  require("nlstools")
  bioVal <- lapply(seq_along(nlsList), function(i){
    bioX <- lapply(seq_along(nlsList[[i]]), function(j){
      if(!is.null(nlsList[[i]][[j]])){
        if(names(nlsList[[i]][[j]]) == "analytis_77"){
          parX <- summary(nlsList[[i]][[j]][[1]])$parameters[4:5, 1:2]
          confX <- nlstools::confint2(nlsList[[i]][[j]][[1]])
          modName <- names(nlsList[[i]][[j]])
          return(data.frame(
            stage = unique(poolData$STAGE)[j],
            sp = unique(poolData$completeSp)[i],
            modName,
            parName = rownames(parX),
            parX,
            confX[4:5, ]
          ))
        }
        if(names(nlsList[[i]][[j]]) == "ratkowsky_83" |
           names(nlsList[[i]][[j]]) == "beta_95" |
           names(nlsList[[i]][[j]]) == "beta_16" |
           names(nlsList[[i]][[j]]) == "briere1_99" |
           names(nlsList[[i]][[j]]) == "briere2_99" |
           names(nlsList[[i]][[j]]) == "kontodimas_04"){
          parX <- summary(nlsList[[i]][[j]][[1]])$parameters[2:3, 1:2]
          confX <- nlstools::confint2(nlsList[[i]][[j]][[1]], )
          modName <- names(nlsList[[i]][[j]])
          return(data.frame(
            stage = unique(poolData$STAGE)[j],
            sp = unique(poolData$completeSp)[i],
            modName,
            parName = rownames(parX),
            parX,
            confX[2:3, ]
          ))
        }
        if(names(nlsList[[i]][[j]]) == "hilbertLogan_83" |
           names(nlsList[[i]][[j]]) == "regniere_12"){
          parX <- summary(nlsList[[i]][[j]][[1]])$parameters[3:4, 1:2]
          confX <- nlstools::confint2(nlsList[[i]][[j]][[1]])
          modName <- names(nlsList[[i]][[j]])
          return(data.frame(
            stage = unique(poolData$STAGE)[j],
            sp = unique(poolData$completeSp)[i],
            modName,
            parName = rownames(parX),
            parX,
            confX[3:4, ]
          ))
        }
        if(names(nlsList[[i]][[j]]) == "shi_11"){
          parX <- summary(nlsList[[i]][[j]][[1]])$parameters[c(3, 5), 1:2]
          confX <- nlstools::confint2(nlsList[[i]][[j]][[1]])
          modName <- names(nlsList[[i]][[j]])
          return(data.frame(
            stage = unique(poolData$STAGE)[j],
            sp = unique(poolData$completeSp)[i],
            modName,
            parName = rownames(parX),
            parX,
            confX[c(3, 5)]
          ))
        }
        if(names(nlsList[[i]][[j]]) == "perf2_11"){
          parX <- summary(nlsList[[i]][[j]][[1]])$parameters[c(2, 4), 1:2]
          confX <- nlstools::confint2(nlsList[[i]][[j]][[1]])
          modName <- names(nlsList[[i]][[j]])
          return(data.frame(
            stage = unique(poolData$STAGE)[j],
            sp = unique(poolData$completeSp)[i],
            modName,
            parName = c("Tmin", "Tmax"),
            parX,
            confX[c(2, 4), ]
          ))
        }
      }else{
        return(NA)
      }
    })
    bioX <- do.call(rbind, bioX)
    rownames(bioX) <- NULL
    return(bioX)
  })
  bioVal <- do.call(rbind, bioVal)
  return(bioVal)
}

#' plotBioVal: plots biological parameters values with precision measurement in
#' function of stage and species
#' @param bioVal
#' @param IC if TRUE plot 95% confidence interval, else plot Standard Errors
plotBioVal <- function(bioVal, IC = TRUE){
  pdf("./results/bioValPlot.pdf", width = 5, height = 4)
  par(mar = c(4, 4, 5, 0))
  plot(NULL, xlim = c(0, 12),
       ylim = c(0, 50),
       panel.first = grid(),
       axes = FALSE,
       xlab = "Stage",
       ylab = "Temperature")
  box(bty = "l")
  axis(side = 1,
       at = c(2, 6, 10),
       labels = c("egg", "larva", "pupa"),
       col = NA, col.ticks = 1)
  axis(side = 2)
  legend(x = 0, y = 80, legend = unique(bioVal$sp),
         pch = 15,
         col = 2:5,
         xpd = TRUE, bty = "n", text.font = 3)
  pchX <- c(16, 17)
  legend(x = 5, y = 80, legend = c("CTmin", "CTmax"),
         pch = pchX,
         xpd = TRUE, bty = "n")
  u <- 0
  for(parNameX in c("Tmin", "Tmax")){
    u <- u + 1
    dfSub <- bioVal[bioVal$parName == parNameX, ]
    k <- 0
    for(i in c(2, 6, 10)){
      k <- k + 1
      z <- 0
      for(j in c(-1, -1/3, 1/3, 1)){
        z <- z + 1
        bioValX <- dfSub[dfSub$stage == unique(dfSub$stage)[k] &
                           dfSub$sp == unique(dfSub$sp)[z], "Estimate"]
        points(x = i+j,
               y = bioValX, pch = pchX[u], col = z+1)
        if(IC == TRUE){
          sub <- dfSub[dfSub$stage == unique(dfSub$stage)[k] &
                         dfSub$sp == unique(dfSub$sp)[z], "X2.5.."]
          sus <- dfSub[dfSub$stage == unique(dfSub$stage)[k] &
                         dfSub$sp == unique(dfSub$sp)[z], "X97.5.."]
        }else{
          sub <- bioValX - dfSub[dfSub$stage == unique(dfSub$stage)[k] &
                                   dfSub$sp == unique(dfSub$sp)[z], "Std..Error"]
          sus <- bioValX + dfSub[dfSub$stage == unique(dfSub$stage)[k] &
                                   dfSub$sp == unique(dfSub$sp)[z], "Std..Error"]
        }
        segments(x0 = i + j,
                 x1 = i + j,
                 y0 = bioValX,
                 y1 = sub,
                 col = z+1)
        segments(x0 = i + j,
                 x1 = i + j,
                 y0 = bioValX,
                 y1 = sus,
                 col = z+1)
      }
    }
  }
  dev.off()
}