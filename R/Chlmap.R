#'@name chlmap
#'@title Create chlorophyll maps
#'@description Create a chlorophyll concentration surface using streaming data and regression against extracted chlorophyll
#'@details
#'A new interpolation is run after calculating an extracted chlorophyll for all streaming observations. Calculated values that exceed the maximum observed grab sample concentration are discarded.
#'@param yearmon numeric date in yyyymm format
#'@param subgroup subgroup of grab sample stations if creating separate regressions based on area
#'@param remove.flags logical Use QA'd grab data?
#'@param stream.qa logical Use QA'd streaming data?
#'@param fdir file.path to data folder
#'@import stats
#'@return An extracted chlorophyll surface and an updated FullDataSet file
#'@author Jemma Stachelek and Kira Allen
#'@importFrom utils read.table
#'@export

chlmap <- function(yearmon, subgroup=NA, stream.qa=T, remove.flags = TRUE, fdir = getOption("fdir")) {
  #find coefficients that match yearmon####
  coeflist <- read.csv(file.path(fdir, "DF_GrabSamples", "extractChlcoef.csv"), header = T, na.strings = "NA")[,-1]
  names(coeflist) <- tolower(names(coeflist))
  model <- coeflist[coeflist$yearmon == yearmon, "model"]
  if(length(model) == 0){
    stop("No chl model fit for this survey")
  }

  if (all(is.na(subgroup))) {
    group <- coeflist[(coeflist$yearmon == yearmon & coeflist$subgroup == "full"),]
    if (nrow(group) > 1) {
      stop("Remove duplicate entries")
    }
  } else {
    groups <- coeflist[(coeflist$yearmon == yearmon & coeflist$subgroup != "full"),]
    if (nrow(groups) > 2) {
      stop("Remove duplicate entries")
    }
    group1 <- groups[1,]
    group2 <- groups[2,]
  }

  calc_chl <- function(group) {
    group_coeffs <- group[, (which(colnames(group)=="subgroup")+1):which(colnames(group)=="intercept")]
    namelist <- names(group_coeffs)[which(!is.na(group_coeffs))]
    group_coeffs <- group_coeffs[!is.na(group_coeffs)]
    namelist_sq <- namelist[grep("2", namelist)]

    namesalias <- read.table(text = "
                       c6chla c6chl
                       chla chlaiv
                         ")

    for(n in 1:length(namelist)){
      if(any(namelist[n] == namesalias[,2])){
        namelist[n] <- as.character(namesalias[which(namelist[n] == namesalias[,2]), 1])
      }
    }

    namelist_temp <- namelist
    for(n in 1:length(namelist)){
      if(any(namelist[n] == namesalias[,1])){
        namelist_temp[n] <- as.character(namesalias[which(namelist[n] == namesalias[,1]),2])
      }
    }

    if(length(namelist_sq) > 0){
      namelist_sq <- sapply(namelist_sq, function(x) substring(x, 1, (nchar(x) - 1)))
      for(n in 1:length(namelist_sq)){
        if(any(namelist_sq[n] == namesalias[,2])){
          namelist_sq[n] <- as.character(namesalias[which(namelist_sq[n] == namesalias[,2]), 1])
        }
      }
    }

    grabs <-  grabget(yearmon)

    if(all(!(is.na(subgroup)))) {
      subnames <- grabs[which(grabs$fathomlocation %in% subgroup),"name"]
      if(group$subgroup == group1$subgroup) {
        grab <- subset(grabs, fathomlocation %in% subgroup)
      } else {
        grab <- subset(grabs, !(fathomlocation %in% subgroup))
      }
    } else {
      grab <- grabs
    }

    if(length(namelist_sq) > 0) {
      fit <- lm(chla ~ poly(namelist_sq, 2), data=grab)
    } else {
      fit <- lm(as.formula(paste("chla ~ ", paste(namelist_temp[1:(length(namelist_temp) - 1)], collapse = "+"))), data = grab)
    }

    if(stream.qa == TRUE){
      dt <- streamget(yearmon, qa = TRUE)
    }else{
      dt <- streamget(yearmon, qa = FALSE)
    }
    if(all(!(is.na(subgroup)))) {
      if(group$subgroup == group1$subgroup) {
        dt <- subset(dt, name %in% subnames)
      } else {
        dt <- subset(dt, !(name %in% subnames))
      }
    }

    dt_temp <- dt[,namelist[1:(length(namelist) - 1)]]
    if(!(length(ncol(dt_temp)) >= 1)){ #handle one variable namelist
      dt_temp <- data.frame(dt_temp)
    }
    names(dt_temp) <- namelist_temp[1:(length(namelist_temp)-1)]

    for (i in 1:ncol(dt_temp)) {
      if (typeof(dt_temp[,i]) != "numeric") {
        dt_temp[,i] <- as.numeric(dt_temp[,i])
      }
    }

    chlext <- predict(fit, dt_temp)
    chlext_low <- predict(fit, dt_temp, se.fit = TRUE)$fit - predict(fit, dt_temp, se.fit = TRUE)$se.fit
    chlext_hi <- predict(fit, dt_temp, se.fit = TRUE)$fit + predict(fit, dt_temp, se.fit = TRUE)$se.fit

    bad_chl <- chlext > range(grabget(yearmon, remove.flags = remove.flags)$chla, na.rm = TRUE)[2]
    chlext[bad_chl] <- NA
    chlext_low[bad_chl] <- NA
    chlext_hi[bad_chl] <- NA

    dt$chlext <- chlext
    dt$chlext_low <- chlext_low
    dt$chlext_hi <- chlext_hi

    dt

  }

  if (all(!(is.na(subgroup)))) {
    dt1 <- calc_chl(group1)
    dt2 <- calc_chl(group2)
    dt <- rbind(dt1, dt2)
    dt <- dt[order(dt$datetime),]
  } else {
    dt <- calc_chl(group)
  }


  if(stream.qa==T) {
    write.csv(dt, file.path(fdir, "DF_FullDataSets", "QA datasets", paste(yearmon, "j_qa.csv", sep = "")), row.names=F)
  } else {
    write.csv(dt, file.path(fdir, "DF_FullDataSets", paste(yearmon, "j.csv", sep = "")), row.names=F)
  }

  if(file.exists(file.path(fdir,paste0("/DF_Subsets/chlext",yearmon,".csv"),fsep=""))){
    file.remove(file.path(fdir,paste0("/DF_Subsets/chlext",yearmon,".csv"),fsep=""))
    file.remove(file.path(fdir,paste0("/DF_Validation/chlext",yearmon,".csv"),fsep=""))
  }

  streaminterp(dt, paramlist = c("chlext", "chlext_low", "chlext_hi"), yearmon = yearmon, tname = file.path(fdir, paste0("/DF_Subsets/chlext", yearmon, ".csv"), fsep = "") ,vname = file.path(fdir, paste0("/DF_Validation/chlext", yearmon, ".csv"), fsep = ""), missprop = (1/3), trim_negative = TRUE)

}





#'@name cdommap
#'@title Create cdom maps
#'@description Create a cdom concentration surface using streaming data and regression against extracted cdom
#'@details
#'A new interpolation is run after calculating an extracted cdom for all streaming observations. Calculated values that exceed the maximum observed grab sample concentration are discarded.
#'@param yearmon numeric date in yyyymm format
#'@param subgroup subgroup of grab sample stations if creating separate regressions based on area
#'@param remove.flags logical Use QA'd grab data?
#'@param stream.qa logical Use QA'd streaming data?
#'@param fdir file.path to data folder
#'@import stats
#'@return An extracted chlorophyll surface and an updated FullDataSet file
#'@author Jemma Stachelek and Kira Allen
#'@importFrom utils read.table
#'@export

cdommap <- function(yearmon, subgroup=NA, stream.qa=T, remove.flags = TRUE, fdir = getOption("fdir")) {
  #find coefficients that match yearmon####
  coeflist <- read.csv(file.path(fdir, "DF_GrabSamples", "extractCdomcoef.csv"), header = T, na.strings = "NA")[,-1]
  names(coeflist) <- tolower(names(coeflist))
  model <- coeflist[coeflist$yearmon == yearmon, "model"]
  if(length(model) == 0){
    stop("No cdom model fit for this survey")
  }

  if (all(is.na(subgroup))) {
    group <- coeflist[(coeflist$yearmon == yearmon & coeflist$subgroup == "full"),]
    if (nrow(group) > 1) {
      stop("Remove duplicate entries")
    }
  } else {
    groups <- coeflist[(coeflist$yearmon == yearmon & coeflist$subgroup != "full"),]
    if (nrow(groups) > 2) {
      stop("Remove duplicate entries")
    }
    group1 <- groups[1,]
    group2 <- groups[2,]
  }

  calc_cdom <- function(group) {
    group_coeffs <- group[, (which(colnames(group)=="subgroup")+1):which(colnames(group)=="intercept")]
    namelist <- names(group_coeffs)[which(!is.na(group_coeffs))]
    group_coeffs <- group_coeffs[!is.na(group_coeffs)]
    namelist_sq <- namelist[grep("2", namelist)]

    namesalias <- read.table(text = "
                       c6chla c6chl
                       chla chlaiv
                         ")

    for(n in 1:length(namelist)){
      if(any(namelist[n] == namesalias[,2])){
        namelist[n] <- as.character(namesalias[which(namelist[n] == namesalias[,2]), 1])
      }
    }

    namelist_temp <- namelist
    for(n in 1:length(namelist)){
      if(any(namelist[n] == namesalias[,1])){
        namelist_temp[n] <- as.character(namesalias[which(namelist[n] == namesalias[,1]),2])
      }
    }

    if(length(namelist_sq) > 0){
      namelist_sq <- sapply(namelist_sq, function(x) substring(x, 1, (nchar(x) - 1)))
      for(n in 1:length(namelist_sq)){
        if(any(namelist_sq[n] == namesalias[,2])){
          namelist_sq[n] <- as.character(namesalias[which(namelist_sq[n] == namesalias[,2]), 1])
        }
      }
    }

    grabs <-  grabget(yearmon)

    if(all(!(is.na(subgroup)))) {
      subnames <- grabs[which(grabs$fathomlocation %in% subgroup),"name"]
      if(group$subgroup == group1$subgroup) {
        grab <- subset(grabs, fathomlocation %in% subgroup)
      } else {
        grab <- subset(grabs, !(fathomlocation %in% subgroup))
      }
    } else {
      grab <- grabs
    }

    if(length(namelist_sq) > 0) {
      fit <- lm(cdom ~ poly(namelist_sq, 2), data=grab)
    } else {
      fit <- lm(as.formula(paste("cdom ~ ", paste(namelist_temp[1:(length(namelist_temp) - 1)], collapse = "+"))), data = grab)
    }

    if(stream.qa == TRUE){
      dt <- streamget(yearmon, qa = TRUE)
    }else{
      dt <- streamget(yearmon, qa = FALSE)
    }
    if(all(!(is.na(subgroup)))) {
      if(group$subgroup == group1$subgroup) {
        dt <- subset(dt, name %in% subnames)
      } else {
        dt <- subset(dt, !(name %in% subnames))
      }
    }

    dt_temp <- dt[,namelist[1:(length(namelist) - 1)]]
    if(!(length(ncol(dt_temp)) >= 1)){ #handle one variable namelist
      dt_temp <- data.frame(dt_temp)
    }
    names(dt_temp) <- namelist_temp[1:(length(namelist_temp)-1)]

    for (i in 1:ncol(dt_temp)) {
      if (typeof(dt_temp[,i]) != "numeric") {
        dt_temp[,i] <- as.numeric(dt_temp[,i])
      }
    }

    cdomext <- predict(fit, dt_temp)
    cdomext_low <- predict(fit, dt_temp, se.fit = TRUE)$fit - predict(fit, dt_temp, se.fit = TRUE)$se.fit
    cdomext_hi <- predict(fit, dt_temp, se.fit = TRUE)$fit + predict(fit, dt_temp, se.fit = TRUE)$se.fit

    bad_cdom <- cdomext > range(grabget(yearmon, remove.flags = remove.flags)$cdom, na.rm = TRUE)[2]
    cdomext[bad_cdom] <- NA
    cdomext_low[bad_cdom] <- NA
    cdomext_hi[bad_cdom] <- NA

    dt$cdomext <- cdomext
    dt$cdomext_low <- cdomext_low
    dt$cdomext_hi <- cdomext_hi

    dt

  }

  if (all(!(is.na(subgroup)))) {
    dt1 <- calc_cdom(group1)
    dt2 <- calc_cdom(group2)
    dt <- rbind(dt1, dt2)
    dt <- dt[order(dt$datetime),]
  } else {
    dt <- calc_cdom(group)
  }


  if(stream.qa==T) {
    write.csv(dt, file.path(fdir, "DF_FullDataSets", "QA datasets", paste(yearmon, "j_qa.csv", sep = "")), row.names=F)
  } else {
    write.csv(dt, file.path(fdir, "DF_FullDataSets", paste(yearmon, "j.csv", sep = "")), row.names=F)
  }

  if(file.exists(file.path(fdir,paste0("/DF_Subsets/cdomext",yearmon,".csv"),fsep=""))){
    file.remove(file.path(fdir,paste0("/DF_Subsets/cdomext",yearmon,".csv"),fsep=""))
    file.remove(file.path(fdir,paste0("/DF_Validation/cdomext",yearmon,".csv"),fsep=""))
  }

  streaminterp(dt, paramlist = c("cdomext", "cdomext_low", "cdomext_hi"), yearmon = yearmon, tname = file.path(fdir, paste0("/DF_Subsets/cdomext", yearmon, ".csv"), fsep = "") ,vname = file.path(fdir, paste0("/DF_Validation/cdomext", yearmon, ".csv"), fsep = ""), missprop = (1/3), trim_negative = TRUE)

}













