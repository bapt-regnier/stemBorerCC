#' readRawData: read raw data from CSV file
#' 
#' Load the raw data CSV file and split the development rate end temperature
#' fields into lists.
#' @param rawDataFile Path to raw data file
readRawData <- function(rawDataFile){
  data <- read.table(rawDataFile, header = TRUE, sep = ";")
  getTEMP <- lapply(as.character(data$TEMP), function(i){
    j <- strsplit(i, split = "\\,")[[1]]
    return(unname(sapply(j, function(x) eval(parse(text = x)))))
  })
  getDEVRATE <- lapply(as.character(data$DEVRATE), function(i){
    j <- strsplit(i, split = "\\,")[[1]]
    return(unname(sapply(j, function(x) eval(parse(text = x)))))
  })
  data$TEMP <- getTEMP
  data$DEVRATE <- getDEVRATE
  rm(getTEMP, getDEVRATE)
  data$completeSp <- paste(data$GENUS, data$SPECIES)
  return(data)
}

#' computeAddStage: compute additional life stages by regrouping all larval 
#' stage into "larva" and the when all life stages are described, adding "all"
#' @param data dataset
computeAddStage <-function(data){
  dfNoLar <- data[data$REF == "Calvin et al. 1991" |
                    data$REF == "Andreadis et al. 2013" |
                    data$REF == "Hilal 1981", ] # only three papers
  addedLar <- lapply(unique(dfNoLar$region), function(x){
    dfX <- dfNoLar[dfNoLar$region == x, ]
    dfX <- dfX[grep(x = dfX$STAGE, pattern = "larva "), ]
    nLar <- nrow(dfX)
    if(length(unique(dfX$TEMP)) == 1){
      tempX <- unlist(unique(dfX$TEMP))
    } else {
      nTemp <- sapply(unique(dfX$TEMP), function(x){
        length(x)
      })
      tempX <- unique(dfX$TEMP)[[which(nTemp == min(nTemp))]]
    }
    devRateX <- lapply(seq_along(tempX), function(k){
      xxx <- lapply(1:nLar, function(z) {
        dfX$DEVRATE[[z]][dfX$TEMP[[z]] == tempX[k]]
      })
      unlist(xxx)
    })
    devRateX <- sapply(devRateX, function(x){sum(1/x)})
    devRateX <- 1 / devRateX
    dfAddedLar <- data.frame(id = NA,
                             dfX[1, 2:7],
                             STAGE = "larva",
                             TEMP = NA,
                             DEVRATE = NA,
                             dfX[1, 11:ncol(dfX)])
    dfAddedLar$TEMP <- list(tempX)
    dfAddedLar$DEVRATE <- list(devRateX)
    return(dfAddedLar)
  })
  addedLar <- do.call(rbind, addedLar)
  data <- rbind(data, addedLar)
  return(data)
}

#' selectStages: select egg larva pupa only
#' @param database filtered database
selectStages <- function(data){
  data <- data[data$STAGE == "egg" |
                 data$STAGE == "larva" |
                 data$STAGE == "pupa", ]
  return(data)
}

#' poolData: pool the data
#' @param data dataset
poolData <- function(data){
  poolData <- lapply(unique(data$completeSp), function(speciesX){
    print(speciesX)
    lapply(unique(data$STAGE), function(stageX){
      print(stageX)
      tempX <- unname(
        unlist(
          data[data$completeSp == speciesX &
                 data$STAGE == stageX, ]$TEMP))
      devX <- unname(
        unlist(
          data[data$completeSp == speciesX &
                 data$STAGE == stageX, ]$DEVRATE))
      return(data.frame(STAGE = stageX,
                        completeSp = speciesX,
                        TEMP = tempX,
                        DEVRATE = devX))
    })
  })
  poolData <- do.call(rbind,
                      lapply(seq_along(poolData),
                             function(i){
                               do.call(rbind, poolData[[i]])
                             }))
  return(poolData)
}

