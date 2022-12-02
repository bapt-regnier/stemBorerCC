#' get temperature data for region of interest and mean of monthly temperatures
#' for present-day conditions (1990-2014) and future conditions (2081-2100)
#' @param source_dir path to temperature datasets
#' @param region name of region of interest (acronym from IPCC reference
#' regions, e.g. "MED", "SEAF", or "ENA")
getTdata <- function(source_dir, region){
  files <- list.files(source_dir)
  TvarList <- lapply(files, function(fileX){
    print(fileX)
    splitFileName <- strsplit(split = "_", x = fileX)[[1]]
    dfX <- read.csv(paste0(source_dir, fileX), h = TRUE, sep = ",", skip = 15)
    tempRegX <- dfX[, c("date", region)]
    
    if(splitFileName[3] == "historical"){
      nRow <- nrow(tempRegX)
      xxx <- lapply(seq(nRow-299, nRow, 12), function(i){ # get month Tmean for 1990-2014
        sub1 <- tempRegX[i:c(i+11), ]
        dateX <- strsplit(sub1$date, "-")
        dateX <- do.call(rbind, dateX)
        return(data.frame(year = dateX[, 1],
                          month = dateX[, 2],
                          Tmean = sub1[, 2]
        ))
      })
      xxx <- do.call(rbind, xxx)
      tempDf <- aggregate(xxx$Tmean, by = list(xxx$month), mean) # mean for 1990-2014
      colnames(tempDf) <- c("month", "tempMed")
      tempDf$model <- splitFileName[2]
      tempDf$scenario <- splitFileName[3]    
    }else{
      xxx <- lapply(seq(793, nrow(tempRegX), 12), function(i){ # get Tmean month for 2081-2100
        sub1 <- tempRegX[i:c(i+11), ]
        dateX <- strsplit(sub1$date, "-")
        dateX <- do.call(rbind, dateX)
        return(data.frame(year = dateX[, 1],
                          month = dateX[, 2],
                          Tmean = sub1[, 2]
        ))
      })
      xxx <- do.call(rbind, xxx)
      tempDf <- aggregate(xxx$Tmean, by = list(xxx$month), mean)  # mean for the period 2081-2100
      colnames(tempDf) <- c("month", "tempMed")
      tempDf$model <- splitFileName[2]
      tempDf$scenario <- splitFileName[3]
    }
    return(tempDf)
  })
  TvarList <- do.call(rbind, TvarList)
  return(TvarList)
}

#' plot results on impacts of global warming on pest development
#' @param tempDataList temperature data
#' @param nlsList lsit of selected nls
#' @param statList list of median development time and 95% prediction interval
plotRTpropDiff <- function(tempDataList, nlsList, statList){
  sp <- c("Chilo partellus",
          "Busseola fusca",
          "Ostrinia nubilalis",
          "Sesamia nonagrioides")
  colX4 <- c(rgb(200, 150, 16, 255*0.7, maxColorValue = 255),
             rgb(205, 11, 188, 255*0.7, maxColorValue = 255))
  colX2 <- c(rgb(200, 150, 16, 255/4, maxColorValue = 255),
             rgb(205, 11, 188, 255/4, maxColorValue = 255))
  colX <- c(rgb(0, 0, 0, 255, maxColorValue = 255),
            rgb(200, 150, 16, 255, maxColorValue = 255),
            rgb(205, 11, 188, 255, maxColorValue = 255))
  colX3 <- c(rgb(0, 0, 0, 255/6, maxColorValue = 255),
             rgb(200, 150, 16, 255/3, maxColorValue = 255),
             rgb(205, 11, 188, 255/6, maxColorValue = 255))
  colSp <- c(rgb(223, 83, 107, 255, maxColorValue = 255),
             rgb(97, 208, 79, 255, maxColorValue = 255),
             rgb(34, 151, 230, 255, maxColorValue = 255),
             rgb(40, 226, 229, 255, maxColorValue = 255))
  pdf("./results/CCimpacts_vBoxplot_rPropDiff.pdf")
  
  split.screen(matrix(data = c(0, 0.48, 0.5, 0.9,
                               0.5, 0.98, 0.5, 0.9,
                               0, 0.48, 0, 0.4,
                               0.5, 0.98, 0, 0.4,
                               
                               0, 0.48, 0.9, 1,
                               0.5, 0.98, 0.9, 1,
                               0, 0.48, 0.4, 0.5,
                               0.5, 0.98, 0.4, 0.5,
                               
                               0.48, 0.5, 0.9, 1,
                               0.98, 1, 0.9, 1,
                               0.48, 0.5, 0.4, 0.5,
                               0.98, 1, 0.4, 0.5), ncol = 4, byrow = TRUE))
  sink("./results/rawRes_CCimpacts.csv")
  cat("minimum;Q1;median;Q3;maximum;mean;month;scenario;species\n")
  lapply(seq_along(tempDataList), function(i){
    screen(i)
    par(cex = 0.8,
        mar = c(4, 4, 0, 0))
    tempData <- tempDataList[[i]]
    xxx <- statList[[i]]
    xxx <- xxx[!is.na(xxx$med),]
    Tmin <- xxx[1, 1]
    Tmax <- xxx[nrow(xxx), 1]
    Topt <- xxx[xxx$med == min(xxx$med, na.rm = TRUE), 1]
    rTmax <- 1/xxx[xxx$med == min(xxx$med, na.rm = TRUE), 2]
    tempDataH <- tempData[tempData$scenario == "historical", ]
    diffList <- lapply(c("ssp126", "ssp585"), function(scenX){
      sub1 <- tempData[tempData$scenario == scenX, ]
      xxx <- lapply(unique(sub1$model), function(modX){
        sub2 <- sub1[sub1$model == modX, ]
        subH <- tempDataH[tempDataH$model == modX, ]
        sub2$month <- as.numeric(sub2$month)
        subH$month <- as.numeric(subH$month)
        nlsX <- list(nlsList[[i]])
        dTpredH <- computeDT(nlsList = nlsX, nInd = 5000, tempSeq = subH$tempMed)
        dTpredF <- computeDT(nlsList = nlsX, nInd = 5000, tempSeq = sub2$tempMed)
        
        getMeanT <- function(temp, dfX, fun){
          return(
            sapply(X = temp, FUN = function(x){
              meanX <- eval(fun(1/dfX$days[dfX$addTemp == x], na.rm = TRUE))
              return(meanX)
            })
          )
        }
        
        rTmeanH <- getMeanT(temp = subH$tempMed, dfX = dTpredH[[1]], fun = mean)
        rTsdH <- getMeanT(temp = subH$tempMed, dfX = dTpredH[[1]], fun = sd)
        rTmeanF <- getMeanT(temp = sub2$tempMed, dfX = dTpredF[[1]], fun = mean)
        rTsdF <- getMeanT(temp = sub2$tempMed, dfX = dTpredF[[1]], fun = sd)
        
        rTmeanH[subH$tempMed <= Tmin | subH$tempMed >= Tmax] <- 0
        rTmeanF[sub2$tempMed <= Tmin | sub2$tempMed >= Tmax] <- 0
        
        rTsdH[subH$tempMed <= Tmin | subH$tempMed >= Tmax] <- 0
        rTsdF[sub2$tempMed <= Tmin | sub2$tempMed >= Tmax] <- 0
        
        rTsdDiff <- rTsdF - rTsdH
        rTmeanDiff <- rTmeanF - rTmeanH
        diffToptF <- (rTmeanF / rTmax) - (rTmeanH / rTmax)
        
        colX <- vector()
        for(k in 1:nrow(sub2)){
          if(sub2$tempMed[k] >= Topt){
            colX[k] <- 2
          }else{
            colX[k] <- 1
          }
        }
        
        return(data.frame(
          month = subH$month,
          rTmeanDiff = rTmeanDiff,
          rTsdDiff = rTsdDiff,
          rTmean = rTmeanF,
          rTsd = rTsdF,
          diffToptF = diffToptF,
          colX = colX,
          temp = sub2$tempMed,
          tempDiff = sub2$tempMed - subH$tempMed,
          model = modX,
          scenario = scenX))
      })
      xxx <- do.call(rbind, xxx)
      return(xxx)
    })
    xxx2 <- do.call(rbind, diffList)
    y <- xxx2$diffToptF[order(xxx2$temp)]
    
    xaxp <- c(1, 12, 11)
    yaxp <-  c(floor(range(y)), 5)
    plot(NULL, xlim = c(0.8, 12.2), ylim = range(y, na.rm = TRUE),
         xlab = "Month", ylab = "r / rmax difference", axes = FALSE,
         frame = TRUE, mgp = c(3, 0.5, 0))
    axis(1, at = seq(xaxp[1], xaxp[2], (xaxp[2]-xaxp[1])/xaxp[3]),
         labels = c("jan",
                    "feb",
                    "mar",
                    "apr",
                    "may",
                    "jun",
                    "jul",
                    "aug",
                    "sep",
                    "oct",
                    "nov",
                    "dec"), las = 2)
    a <- axis(2, mgp = c(3, 1, 0), las = 2)
    abline(v=seq(xaxp[1], xaxp[2], (xaxp[2]-xaxp[1])/xaxp[3]), lty=3, col = "grey")
    abline(h=a, lty=3,  col = "grey")
    abline(h = 0, col = "darkgrey", lty = 2)
    
    lapply(seq_along(diffList), function(j){
      
      if(j == 1) addX <- -0.2
      if(j == 2) addX <- 0.2
      
      meanX <- vector()
      c <- 0
      for(m in unique(diffList[[j]]$month)){
        c <- c + 1
        meanX[c] <- mean(diffList[[j]]$diffToptF[diffList[[j]]$month == m])
        b <- boxplot(diffList[[j]]$diffToptF[diffList[[j]]$month == m],
                     at = as.numeric(m)+addX,
                     col = colX[j+1], add = TRUE, cex = 0.5, range = 0,
                     boxwex=0.7, outcol=colX[j+1], axes = FALSE,
                     medlwd = 2, boxcol = NA, whiskcol=colX[j+1], whisklty = 1,
                     staplecol = colX[j+1])
        cat(paste0(paste(c(t(b$stats),meanX[c],m,unique(diffList[[j]]$scenario),sp[i]),
                         collapse=";"), "\n"))
        
        nX <- length(diffList[[j]]$diffToptF[diffList[[j]]$month == m])
      }
      if(j == 1) colMean <- "#FFCC33"
      if(j == 2) colMean <- "#FF99FF"
      lines(meanX ~ c(c(1:12)+addX), col = colMean, lty = 1, lwd = 1.5)
      points(meanX ~ c(c(1:12)+addX), pch = 4, col = colMean, cex = 0.8)
    })
    if(i == 3 | i == 4){
      rect(xleft = 0, xright = 3.5, ybottom = -100,ytop = 100,
           col = rgb(0, 0, 0, 255/10, maxColorValue = 255), border = NA)
      rect(xleft = 9.5, xright = 13, ybottom = -100,ytop = 100,
           col = rgb(0, 0, 0, 255/10, maxColorValue = 255),border = NA)
    }else{
      rect(xleft = 4.5, xright = 9.5, ybottom = -100,ytop = 100,
           col = rgb(0, 0, 0, 255/10, maxColorValue = 255), border = NA)
    }
    tempXXX <- aggregate(tempDataList[[i]]$tempMed,
                         by = list(tempDataList[[i]]$month,
                                   tempDataList[[i]]$scenario),
                         quantile, probs = c(0, 0.5, 1))
    j <- 0
    screen(i + 4)
    par(cex = 0.8, mar = c(0.2, 4, 0.2, 0))
    plot(NULL, xlim = c(0.8, 12.2),
         ylim = c(-6, 39),
         xlab = NA, ylab = "T°C", mgp = c(3, 0.5, 0), axes = F)
    xaxp <- c(1, 12, 11)
    yaxp <-  c(-5, 40, 2)
    abline(v=seq(xaxp[1], xaxp[2], (xaxp[2]-xaxp[1])/xaxp[3]), lty=3, col = "grey")
    abline(h=seq(yaxp[1], yaxp[2], (yaxp[2]-yaxp[1])/yaxp[3]), lty=3,  col = "grey")
    box()
    for(scenX in c("historical", "ssp126", "ssp585")){
      j <- j + 1
      sub1 <- tempXXX[tempXXX$Group.2 == scenX, ]
      sub1$Group.1 <- as.numeric(sub1$Group.1)
      lines(sub1$x[, 2] ~ sub1$Group.1, col = colX[j], lwd = 1.5)
      polygon(y = c(sub1$x[, 1][sub1$Group.1 == 1],
                    sub1$x[, 1],
                    sub1$x[, 1][sub1$Group.1 == 12],
                    rev(sub1$x[, 3]),
                    sub1$x[, 3][sub1$Group.1 == 1]),
              x = c(min(sub1$Group.1),
                    sub1$Group.1,
                    max(sub1$Group.1),
                    rev(sub1$Group.1),
                    min(sub1$Group.1)), col = colX3[j], lty = 0)
      axis(side = 2, yaxp = yaxp, mgp = c(1, 1, 0), las = 2)
    }
    if(i == 3 | i == 4){
      rect(xleft = 0, xright = 3.5, ybottom = -100,ytop = 100,
           col = rgb(0, 0, 0, 255/10, maxColorValue = 255), border = NA)
      rect(xleft = 9.5, xright = 13, ybottom = -100,ytop = 100,
           col = rgb(0, 0, 0, 255/10, maxColorValue = 255),border = NA)
    }else{
      rect(xleft = 4.5, xright = 9.5, ybottom = -100,ytop = 100,
           col = rgb(0, 0, 0, 255/10, maxColorValue = 255), border = NA)
    }
    abline(h = Topt, col = 3, lty = 2, lwd = 1.5)
    abline(h = Tmax, col = 2, lty = 2, lwd = 1.5)
    abline(h = Tmin, col = 4, lty = 2, lwd = 1.5)
    screen(i + 8)
    par(mar = c(0.2, 0.2, 0.2, 0.1), xpd = TRUE)
    plot(NULL, 
         ylim = c(-6, 39),
         xlim = c(0, max(1/statList[[i]]$med, na.rm = TRUE)), axes = FALSE)
    lines(statList[[i]]$temp ~ c(1/statList[[i]]$med), col = colSp[i], lwd = 1.5)
  })
  sink()
  dev.off()
}