# _targets.R file
library("targets")
library("future")
library("future.callr")
plan("callr")
source("scripts/functions/functions.R")

tar_option_set(
  packages = c(
    "devRate",
    "nlstools",
    "truncnorm",
    "parallel"
  )
) # packages

list(
  tar_target(
    rawDataFile,
    "./data/littData_phenoPestCC_05.csv", 
    format = "file"
  ),
  tar_target(rawData, readRawData(rawDataFile)),
  tar_target(rawDataC, computeAddStage(rawData)),
  tar_target(rawDataCS, selectStages(rawDataC)),
  tar_target(pooledData, poolData(rawDataCS)),
  tar_target(AIOlist, fitNLS(pooledData)),
  tar_target(paramVal, extractParam(AIOlist, pooledData)),
  tar_target(nlsList, selectNLS(AIOlist, pooledData)),
  tar_target(
    TPCplots, plotTPC(
      poolData = pooledData,
      nlsList = nlsList
    )
  ),
  tar_target(bioVal, extractBioVal(nlsList, pooledData)),
  tar_target(bioValPlot, plotBioVal(bioVal, IC = FALSE)),
  tar_target(statSelect, selectStats(pooledData)),
  tar_target(nFiltered, extractNfiltered(statSelect)),
  tar_target(
    ibmList, computeDT(
      nlsList = nlsList,
      nInd = 5000,
      tempSeq = seq(0, 50, 0.002)
    )
  ),
  tar_target(statList, getDTstats(ibmList = ibmList)),
  tar_target(
    plotEffectStudy,
    plotStudyEffect(
      data = rawDataCS,
      path = "./results/studyEffect.pdf"
    )
  ),
  tar_target(
    plotResidualsVsFitted,
    plotResiduals(
      nlsList,
      pooledData,
      path = "./results/residualsVsFitted.pdf"
    )
  ),
  tar_target(
    tempDataList, list(
      getTdata(source_dir = "./data/CMIP6_tas_land/", region = "SEAF"),
      getTdata(source_dir = "./data/CMIP6_tas_land/", region = "SEAF"),
      getTdata(source_dir = "./data/CMIP6_tas_land/", region = "ENA"),
      getTdata(source_dir = "./data/CMIP6_tas_land/", region = "MED")
    )
  ),
  tar_target(plotRpropDiff, plotRTpropDiff(tempDataList, nlsList, statList)),
  tar_target(
    plotCompleteDevAllSp, plotDT(
      nlsList, statList, bioVal, maxLim = 182
    )
  )
)
