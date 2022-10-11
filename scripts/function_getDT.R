getDTx <- function(Tx, eqName, models, nInd){
  require("truncnorm")
  if(eqName[1] == "perf2_11"){
    Tmin <- coef(models[[1]])[grep(pattern = "T1",
                                   x = names(coef(models[[1]])))]
    Tmax <- coef(models[[1]])[grep(pattern = "T2",
                                   x = names(coef(models[[1]])))]
  }else{
    Tmin <- coef(models[[1]])[grep(pattern = "Tmin",
                                   x = names(coef(models[[1]])))]
    Tmax <- coef(models[[1]])[grep(pattern = "Tmax",
                                   x = names(coef(models[[1]])))]
  }
  if(Tx > Tmin & Tx < Tmax){
    rM <- predict(models[[1]], newdata = list(T = Tx))
    rS <- 0.15 * predict(models[[1]], newdata = list(T = Tx))
    De <- rtruncnorm(n = nInd,
                     a = 0,
                     b = Inf,
                     mean = rM,
                     sd = rS)
  }else{
    De <- NA
  }
  if(eqName[2] == "perf2_11"){
    Tmin <- coef(models[[2]])[grep(pattern = "T1",
                                   x = names(coef(models[[2]])))]
    Tmax <- coef(models[[2]])[grep(pattern = "T2",
                                   x = names(coef(models[[2]])))]
  }else{
    Tmin <- coef(models[[2]])[grep(pattern = "Tmin",
                                   x = names(coef(models[[2]])))]
    Tmax <- coef(models[[2]])[grep(pattern = "Tmax",
                                   x = names(coef(models[[2]])))]
  }
  if(Tx > Tmin & Tx < Tmax){
    rM <- predict(models[[2]], newdata = list(T = Tx))
    rS <- 0.15 * predict(models[[2]], newdata = list(T = Tx))
    Dl <- rtruncnorm(n = nInd,
                     a = 0,
                     b = Inf,
                     mean = rM,
                     sd = rS)
  }else{
    Dl <- NA
  }
  if(eqName[3] == "perf2_11"){
    Tmin <- coef(models[[3]])[grep(pattern = "T1",
                                   x = names(coef(models[[3]])))]
    Tmax <- coef(models[[3]])[grep(pattern = "T2",
                                   x = names(coef(models[[3]])))]
  }else{
    Tmin <- coef(models[[3]])[grep(pattern = "Tmin",
                                   x = names(coef(models[[3]])))]
    Tmax <- coef(models[[3]])[grep(pattern = "Tmax",
                                   x = names(coef(models[[3]])))]
  }
  if(Tx > Tmin & Tx < Tmax){
    rM <- predict(models[[3]], newdata = list(T = Tx))
    rS <- 0.15 * predict(models[[3]], newdata = list(T = Tx))
    Dp <- rtruncnorm(n = nInd,
                     a = 0,
                     b = Inf,
                     mean = rM,
                     sd = rS)
  }else{
    Dp <- NA
  }
  D <- 1/De + 1/Dl + 1/Dp
  tempX <- rep(Tx, length(D))
  return(data.frame(days = D, addTemp = tempX))
}