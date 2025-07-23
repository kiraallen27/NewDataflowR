#'@name streamclean
#'@title Cleaning raw streaming Dataflow output
#'@description Cleaning raw streaming Dataflow, C6, Eureka-Manta, and YSI-Exo output
#'@param yearmon numeric designation of survey date formatted as yyyymm
#'@param gps character dataset to use for GPS alignment of other data streams. Choice of "df", "eu", or "exo".
#'@param dfmmin integer optional minimum df measurement frequency (# measurements/min)
#'@param c6mmin integer optional minimum c6 measurement frequency
#'@param eummin integer optional minimum eureka (manta) measurement frequency
#'@param exommin integer optional minimum exo measurement frequency
#'@param tofile logical save cleaned output to DF_FullDataSets?
#'@param sep character optional predesignation of item seperation character in raw data files
#'@param fdir character file path to local data directory
#'@export
#'@importFrom zoo zoo na.approx
#'@importFrom sp coordinates CRS spTransform over
#'@importFrom sf read_sf
#'@importFrom methods as
#'@details Dataflow cleaning drops all minutes that have less measurements than "mmin". C6 data is interpolated to match Dataflow.  Automatically compares salinity against conducitivty/temperature recalculated salinity and replaces if slope of fit is not close to 1. Bad DO columns must sometimes be removed manually. TODO - Add check the make sure that the year of the data (not just the filename) matches the year of yearmon
#'@examples \dontrun{
#'dt <- streamclean(yearmon = 201512, gps = "df", c6mmin = 6, dfmmin = 7)
#'dt <- streamclean(yearmon = 201601, gps = "eu", eummin = 12)
#'dt <- streamclean(yearmon = 201603, gps = "exo", exommin = 40, c6mmin = 12)
#'}
streamclean <- function (yearmon, gps, dfmmin = NA, c6mmin = NA, eummin = NA,
                         exommin = NA, tofile = FALSE, sep = ",", fdir = getOption("fdir"))
{
  options(warn = -1)
  fdir_fd <- file.path(fdir, "DF_FullDataSets", "Raw", "InstrumentOutput")
  flist <- list.files(fdir_fd, include.dirs = T, full.names = T)
  flist <- flist[substring(basename(flist), 1, 6) == yearmon]
  dflist <- list.files(flist, pattern = c(".*.txt"), include.dirs = T,
                       full.names = T)
  if (length(dflist) == 0) {
    dflist <- list.files(flist, pattern = c(".*.TXT"), include.dirs = T,
                         full.names = T)
  }
  if (length(dflist) == 0) {
    dflist <- list.files(flist, pattern = c(".*DF.csv"),
                         include.dirs = T, full.names = T)
  }
  c6list <- list.files(flist, pattern = c(".*C6.csv"), include.dirs = T,
                       full.names = T)
  if (length(c6list) == 0) {
    c6list <- list.files(flist, pattern = c(".*C6.CSV"),
                         include.dirs = T, full.names = T)
  }
  eulist <- list.files(flist, pattern = c(".*eu.csv"), include.dirs = T,
                       full.names = T)
  if (length(eulist) == 0) {
    eulist <- list.files(flist, pattern = c(".*eu.CSV"),
                         include.dirs = T, full.names = T)
  }
  exolist <- list.files(flist, pattern = c(".*exo.csv"), include.dirs = T,
                        full.names = T)
  if (length(exolist) == 0) {
    exolist <- list.files(flist, pattern = c(".*exo.CSV"),
                          include.dirs = T, full.names = T)
  }
  survey_days <- unique(sapply(basename(c(dflist, eulist,
                                          exolist, c6list)), function(x) substring(basename(x),
                                                                                   1, 8)))
  iterate_days_load <- function(survey_days, file_listing,
                                reading_function, cleaning_function) {
    reslist <- list()
    for (i in 1:length(survey_days)) {
      day_data <- reading_function(file_listing[i])
      day_data <- cleaning_function(day_data)
      reslist[[i]] <- day_data
    }
    do.call("rbind", reslist)
  }
  read_df <- function(dfpath) {
    sep <- ","
    dt <- read.csv(dfpath, skip = 0, header = F, sep = sep,
                   strip.white = TRUE)
    if (suppressWarnings(any(nchar(gsub("\t", "", dt[1, ])) <
                             nchar(as.character(dt[1, ]))))) {
      sep <- "\t"
      dt <- read.csv(dfpath, skip = 0, header = F, sep = sep,
                     stringsAsFactors = FALSE)
    }
    if (ncol(dt)==1) {
      sep <- "\t"
      dt <- read.csv(dfpath, skip = 0, header = F, sep = sep,
                     strip.white = TRUE)
    }
    if (yearmon==199602) {
      dt[,1] <- as.integer(dt[,1])
    }#added because wasn't recognizing date for this file as integer
    fskip <- 1
    while (!(!(class(dt[, 1]) != "integer") | !(class(dt[,
                                                         1]) != "numeric")) | dt[1,1] ==0) {
      dt <- read.csv(dfpath, skip = fskip, header = F,
                     sep = sep, stringsAsFactors = FALSE)
      if (!any(!is.na(dt[, 1])) | mean(nchar(as.character(na.omit(dt[,
                                                                     1])))) < 1 | sum(is.na(dt[, 1])) > (nrow(dt)/2) |
          sum(nchar(gsub("_", "", as.character(na.omit(dt[, 1])))) -
              nchar(as.character(na.omit(dt[, 1])))) != 0) {
        dt <- dt[, -1]
      }
      fskip <- fskip + 1
      if (fskip > 20) {
        stop(paste("Cannot find beginning of measurements!",
                   dfpath))
      }
    }
    sep <- ","
    dt
  }
  clean_df <- function(dt){
    # print(names(dt))
    #remove bad columns
    if(class(dt[,3]) == "integer"){
      dt <- dt[,-3]
      print("removing existing seconds column")
    }#remove existing sec column

    #remove bad columns of all 0 or NA
    if(ncol(dt) > 14){
      dtno.na <- dt[complete.cases(dt[,1:12]),]
      if (nrow(dtno.na) >0 ) {
        dt <- dt[,apply(dtno.na, 2, function(x) abs(sum(as.numeric(x), na.rm = T)) != 0)]
      }
    }

    #temp should never be less than 10, these are likely 'bad' DO columns?
    if(mean(as.numeric(dt[,4]), na.rm = T) < 10 & mean(as.numeric(dt[,5]), na.rm = T) < 10){
      dt <- dt[,-4:-5]
    }

    #print(names(dt))
    rm1 <- as.integer(which(apply(dt, 2, function(x) abs(sum(as.numeric(x), na.rm = T)) > 38)==F))
    if (length(rm1) > 0) {
      badcol <- dt[,rm1]
      dt <- dt[,-rm1]#take out all 0 (22 or 38 is an arbitrary "tolerance" value)
    }
    ones <- apply(dt, 2, function(x) sd(as.numeric(x)[as.numeric(x) != 0 & !is.na(as.numeric(x))])) != 0
    ones[is.na(ones)] <- TRUE
    ones[1:2] <- TRUE
    dt <- dt[,ones]
    dt[,ncol(dt)] <- trimws(dt[,ncol(dt)])
    dt <- dt[,apply(dt, 2, function(x) mean(nchar(x), na.rm = T)) >= 3.0] #(3 is an arbitrary "tolerance" value; accounts for a 3 digit timestamp?)
    dt <- dt[,apply(dt[,3:ncol(dt)], 2, function(x) length(unique(x)) != 3)]

    if (ncol(dt)==10) {
      names(dt) <- c("date", "time", "chla", "temp", "cond", "sal", "trans", "cdom", "lat_dd", "lon_dd")
    } else if (ncol(dt) < 10) {
      names <- c("date", "time", "chla", "temp", "cond", "sal", "trans", "cdom", "lat_dd", "lon_dd")
      names(dt) <- names[-rm1]
      dt[,names[rm1]] <- badcol
      dt <- dt[, names]
    } else {
      if (yearmon==200302) {
        names(dt) <- c("date", "time", "fluor", "ph", "temp", "cond", "sal", "flow", "domuwpar", "turb", "lat_dd", "lon_dd")
      } else if (yearmon==199705) {
        names(dt) <- c("date", "time", "fluor", "ph", "temp", "cond", "sal", "flow", "par", "turb", "time2", "lat_dd", "lon_dd")
      } else if(yearmon==199609 | yearmon==199606 | yearmon==199507 | yearmon==199505) {
        names(dt) <- c("date", "time", "fluor", "ph", "temp", "cond", "sal", "flow", "time2", "lat_dd", "lon_dd")
      } else if (yearmon==199602 | yearmon==199508) {
        names(dt) <- c("date", "time", "fluor", "ph", "temp", "cond", "sal", "flow", "par", "uwpar", "ratio", "time2", "lat_dd", "lon_dd")
      } else if (yearmon==199504 | yearmon==199412) {
        names(dt) <- c("date", "time", "fluor", "temp", "cond", "sal", "flow", "par", "uwpar", "ratio", "time2", "lat_dd", "lon_dd")
      } else {
        dt <- dt[, 1:10]
        names(dt) <- c("date", "time", "chla", "temp", "cond", "sal", "trans", "cdom", "lat_dd", "lon_dd")
      }
    }

    #convert factors to numeric
    dt <- data.frame(as.matrix(dt))
    factorToNumeric <- function(f) as.numeric(levels(f))[f]
    #check to make sure that there are any factor class columns
    if(any(sapply(dt,class) == "factor")){
      dt <- data.frame(sapply(dt, factorToNumeric))
    }

    if(typeof(dt$lat_dd)=="character") {
      dt$lat_dd <- as.numeric(dt$lat_dd)
    }
    if(typeof(dt$lon_dd)=="character") {
      dt$lon_dd <- as.numeric(dt$lon_dd)
    }

    #fix lon lat formatting
    if(mean(nchar(as.character(round(dt[,"lat_dd"]))), na.rm = TRUE) != 2){
      lat <- dt[,"lat_dd"]
      latdeg <- as.numeric(substr(lat, 0, 2))
      latmin <- as.numeric(substr(lat, 3, 8))
      dt[,"lat_dd"] <- latdeg + latmin / 60
      lon <- dt[,"lon_dd"]
      londeg <- as.numeric(substr(lon, 0, 2))
      lonmin <- as.numeric(substr(lon, 3, 8))
      dt[,"lon_dd"] <- (londeg + lonmin / 60) * -1
    }

    dt$time <- as.numeric(dt$time)
    dt$date <- as.numeric(dt$date)

    #remove rows of all NA values
    dt <- dt[as.numeric(rowSums(is.na(dt))) < ncol(dt) - 1,]
    dt <- dt[as.numeric(rowSums(is.na(dt[,c("lat_dd", "lon_dd")]))) < 2,]

    #remove unrealistic coordinates
    dt <- dt[abs(dt$lat_dd) > 24.5 & abs(dt$lat_dd) < 25.5, ]
    dt <- dt[abs(dt$lon_dd) > 80.1 & abs(dt$lon_dd) < 82, ]

    dt <- dt[which(!is.na(dt$date)),]

    #check for incomplete minutes
    datelist <- unique(dt$date)
    reslist2 <- list()
    for(j in 1:length(datelist)){
      #j<-1
      curdat <- dt[dt$date == datelist[j],]
      gdata <- data.frame(table(curdat$time))
      fdata <-as.numeric(as.character(gdata[gdata$Freq < dfmmin, 1]))#too few measurements
      odata <- as.numeric(as.character(gdata[gdata$Freq > dfmmin, 1]))#too many measurements
      if(length(odata) > 0){
        for(k in 1:length(odata)){
          #k<-1
          leng <- nrow(curdat[curdat$time == odata[k],])
          remo <- sample(as.numeric(row.names(curdat[curdat$time == odata[k],])), leng - dfmmin)
          curdat <- curdat[-match(remo, as.numeric(row.names(curdat))),]
        }
      }
      curdat <- curdat[!curdat$time %in% fdata,]
      curdat <- curdat[!is.na(curdat$time),]
      curdat$sec <- rep(round(seq(from = 0, to = 60 - 60 / dfmmin, by = 60 / dfmmin), 3), times = nrow(data.frame(table(curdat$time))))
      #curdat$sec <- rep(round(seq(from = 0, to = 60 - 60 / dfmmin, by = 60 / dfmmin), 3), times = nrow(curdat)/dfmmin)

      reslist2[[j]] <- curdat
    }
    dt <- do.call("rbind", reslist2)

    #detect when dt is measured at fractional seconds
    if(!identical(round(dt$sec), dt$sec)){
      dt$sec <- round(dt$sec)
    }

    #create POSIXct datetime column
    yr <- substring(dt$date, nchar(dt$date) - 1, nchar(dt$date))
    day <- substring(dt$date, nchar(dt$date) - 3, nchar(dt$date) - 2)
    mon <- substring(dt$date, 1, nchar(dt$date) - 4)
    hr <- substring(dt$time, 1, nchar(dt$time) - 2)
    min <- substring(dt$time, nchar(dt$time) - 1, nchar(dt$time))

    if(mean(nchar(mon)) == 1){mon <- paste("0", mon, sep = "")}
    dt$datetime <- paste(yr, "-", mon, "-", day, "-", hr, "-", min, "-", dt$sec, sep = "")
    rm(min)
    dt$datetime <- as.POSIXct(strptime(dt$datetime, format = "%y-%m-%d-%H-%M-%S"))

    #clean data frame
    #trim beginning and end based on when data is all zeros
    trimdt <- function(dt){
      j <- 1
      for(i in 1:nrow(dt)){
        if(dt[i,1:9][order(dt[i,1:9])][2] > 0){
          break
        }
        j <- i + 1
      }
      k <- nrow(dt)
      for(i in nrow(dt):1){
        if(!is.na(min(dt[i, 1:9])) > 0){
          break
        }
        k <- i - 1
      }
      dt[j:k,]
    }

    #issues with trimdt function and not really necessary, so not applying it
    #dt <- trimdt(dt)

    #check for correct cond to salinity calculations
    if(typeof(dt$cond)=="character") {
      dt$cond <- as.numeric(dt$cond)
    }
    if(typeof(dt$temp)=="character") {
      dt$temp <- as.numeric(dt$temp)
    }
    corsal <- NewDataflowR::cond2sal(dt$cond * 1000, dt$temp)
    if((lm(corsal ~ dt$sal)$coefficients[2] - 1) > 0.02){
      dt$sal <- corsal
    }

    #print(paste(basename(dflist[i]), "processed", sep = " "))
    dt
  }
  read_c6 <- function(c6path) {
    c6 <- read.csv(c6path, header = F)
    c6
  }
  clean_c6 <- function(c6) {
    #Added below because newer c6 has parameters listed in different order
    rmvrows <- which(c6[,1] == "Units")
    cnames <- c6[(rmvrows-1),] #Pull out column names
    c6 <- c6[-(1:rmvrows),] #Removes rows without data
    #Rename columns
    cnames[which(cnames == "Channel Name")] <- "datetime"
    cnames[which(cnames == "Chlorophyll a")] <- "c6chla"
    cnames[which(cnames == "Red Chlorophyll")] <- "c6chlaR"
    cnames[which(cnames == "Phycoerythrin")] <- "phycoe"
    cnames[which(cnames == "Phycocyanin")] <- "phycoc"
    cnames[which(cnames == "CDOM")] <- "c6cdom"
    cnames[which(cnames == "Turbidity")] <- "c6turbidity"
    cnames[which(cnames == "Depth")] <- "depth"
    cnames[which(cnames == "Temp")] <- "c6temp"
    colnames(c6) <- cnames
    c6 <- c6[,-(which(colnames(c6)=="depth"))] #Remove depth column

    if (nrow(c6)==0) {
      c6[1,] <- NA
    }

    if (all(is.na(as.POSIXct(strptime(c6$datetime, "%m/%d/%y %H:%M:%S"))))) {
      c6sec <- unlist(lapply(rle(sapply(c6$datetime, function(x) strftime(strptime(x,
                                                                                   format = "%m/%d/%Y %H:%M"), format = "%M")))$lengths,
                             function(x) seq(0, 60 - (60/x), length.out = x)))
      c6$datetime <- as.POSIXct(strptime(paste0(c6$datetime,
                                                ":", c6sec), "%m/%d/%Y %H:%M:%S"))
    } else {
      c6$datetime <- as.POSIXct(strptime(c6$datetime,
                                         "%m/%d/%y %H:%M:%S"))
    }
    c6$sec <- as.numeric(format(c6$datetime, "%S"))
    c6freq <- c6$sec[2] - c6$sec[1]
    c6$datetime <- as.POSIXct(c6$datetime)
    return(c6)
  }
  read_eu <- function(eupath) {
    read.csv(eupath, header = TRUE, stringsAsFactors = FALSE)
  }
  clean_eu <- function(eu) {
    names(eu) <- tolower(make.names(names(eu)))
    eu$datetime <- paste(sapply(eu$date, function(x) mdy2mmyyyy(x)),
                         eu$time)
    if (all(is.na(as.POSIXct(strptime(eu$datetime, "%m/%d/%Y %H:%M:%S"))))) {
      eusec <- unlist(lapply(rle(sapply(eu$date, function(x) strftime(strptime(x,
                                                                               format = "%m/%d/%Y %H:%M"), format = "%M")))$lengths,
                             function(x) seq(0, 60 - (60/x), length.out = x)))
      eu$date <- as.POSIXct(strptime(paste0(eu$date, ":",
                                            eusec), "%m/%d/%Y %H:%M:%S"))
    }
    else {
      eu$datetime <- as.POSIXct(strptime(eu$datetime,
                                         "%m/%d/%Y %H:%M:%S"))
      eu$datetime <- as.POSIXct(strptime(eu$datetime,
                                         "%Y-%m-%d %H:%M:%S"))
      eu <- eu[!is.na(eu$datetime) & nchar(eu$datetime) ==
                 10, ]
    }
    return(eu)
  }
  read_exo <- function(exopath) {
    exo <- read.csv(exopath, header = F, stringsAsFactors = FALSE)
    rmvrows <- which(exo[,1]=="Date")
    cnames <- exo[rmvrows,]
    exo <- exo[-(1:rmvrows),]
    colnames(exo) <- cnames
    exo
  }
  clean_exo <- function(exo) {
    names(exo) <- tolower(gsub("\\.", "", names(exo)))
    exo <- exo[, !(names(exo) %in% c("timefractsec", "sitename",
                                     "x"))]
    exo$datetime <- as.POSIXct(strptime(paste(exo[, grep("date",
                                                         names(exo))], exo[, grep("time", names(exo))]),
                                        format = "%m/%d/%Y %H:%M:%S"))
    exo$sec <- strftime(exo$datetime, "%S")
    exo <- exo[which(!duplicated(exo$datetime)), ]
    if (any(exo[, "lon"] > 0)==TRUE) {
      exo[, "lon"] <- exo[, "lon"] *
        -1
    }
    exo <- exo[exo$lat > 24, ]
    exo
  }
  if (!is.na(dfmmin)) {
    df <- iterate_days_load(survey_days, dflist, read_df,
                            clean_df)
  }
  if (!is.na(c6mmin)) {
    c6 <- iterate_days_load(survey_days, c6list, read_c6,
                            clean_c6)
  }
  if (!is.na(eummin)) {
    eu <- iterate_days_load(survey_days, eulist, read_eu,
                            clean_eu)
  }
  if (!is.na(exommin)) {
    exo <- iterate_days_load(survey_days, exolist, read_exo,
                             clean_exo)
  }
  streams <- c("df", "c6", "eu", "exo")
  streams <- streams[sapply(streams, function(x) exists(x))]
  streams <- streams[sapply(streams, function(x) is.data.frame(eval(as.symbol(x))))]
  stream_mmin <- c("dfmmin", "c6mmin", "eummin", "exommin")
  target <- eval(as.symbol(gps))
  target_mmin <- eval(as.symbol(stream_mmin[grep(gps, stream_mmin)]))
  detect_mmin <- function(dt) {
    round(60/Mode(diff(as.numeric(format(dt$datetime, "%S")))))
  }
  check_correct_mmin <- function(dt, target, target_mmin) {
    dt_mmin <- detect_mmin(dt)
    begin_date <- as.POSIXct(as.Date(unique(strftime(min(na.omit(dt$datetime)),
                                                     format = "%Y-%m-%d"))[1]))
    end_date <- as.POSIXct(as.Date(unique(strftime(max(na.omit(dt$datetime)),
                                                   format = "%Y-%m-%d"))[1]) + 1)
    target_dates <- as.Date(seq(range(target$datetime)[1],
                                range(target$datetime)[2], 86400))
    if (abs(as.Date(begin_date) - min(target_dates)) > 2) {
      begin_date <- as.POSIXct(as.Date(min(target_dates)) -
                                 2)
      dt <- dt[dt$datetime > (begin_date - 86400), ]
    }
    if (abs(as.Date(end_date) - max(target_dates)) > 2) {
      end_date <- as.POSIXct(as.Date(max(target_dates)) +
                               2)
      dt <- dt[dt$datetime < (end_date + 86400), ]
    }
    dt_zoo <- zoo::zoo(dt, dt$datetime)
    dt_zoo_full <- data.frame(seq(begin_date, end_date,
                                  1))
    names(dt_zoo_full) <- "datetime"
    dt_zoo_full <- zoo::zoo(dt_zoo_full, dt_zoo_full$datetime)
    #New code to remove any duplicate rows
    if (any(duplicated(dt_zoo$datetime)==TRUE)) {
      dt_zoo <- dt_zoo[-c(anyDuplicated(dt_zoo$datetime)),]
    }
    if (any(duplicated(dt_zoo_full$datetime)==TRUE)) {
      dt_zoo_full <- dt_zoo_full[-c(anyDuplicated(dt_zoo_full$datetime)),]
    }
    ##
    dt_zoo_full <- merge(dt_zoo_full, dt_zoo)
    dt_zoo_full <- dt_zoo_full[min(which(!is.na(dt_zoo_full[,
                                                            2]))):max(which(!is.na(dt_zoo_full[, 2]))), ]
    is_interp_column <- colSums(!is.na(dt_zoo_full)) < nrow(dt_zoo_full)
    is_interp_column[c(1:2)] <- FALSE
    dt_zoo_full <- zoo::na.approx(dt_zoo_full[, is_interp_column])
    dt <- data.frame(dt_zoo_full, zoo::index(dt_zoo_full),
                     row.names = NULL, stringsAsFactors = FALSE)
    names(dt)[ncol(dt)] <- "datetime"
    dt[, 1:(ncol(dt) - 1)] <- apply(dt[, 1:(ncol(dt) - 1)],
                                    2, function(x) as.numeric(as.character(x)))
    dt_names <- names(dt)
    dt <- merge(target, dt, by = "datetime", all.x = T)
    dt[, dt_names[dt_names %in% names(dt)]]
  }
  if (!is.na(c6mmin) & (all(gps != "c6"))) {
    c6 <- c6[which(!is.na(c6$datetime)),]
    c6 <- check_correct_mmin(c6, target, target_mmin)
  }
  if (!is.na(dfmmin) & (all(gps != "df"))) {
    df <- check_correct_mmin(df, target, target_mmin)
  }
  if (!is.na(eummin) & (all(gps != "eu"))) {
    eu <- check_correct_mmin(eu, target, target_mmin)
  }
  if (!is.na(exommin) & (all(gps != "exo"))) {
    exo <- check_correct_mmin(exo, target, target_mmin)
  }
  contributing_streams <- streams[!(streams %in% gps)]
  dt <- target
  for (i in contributing_streams) {
    dt <- merge(dt, eval(as.symbol(i)), all.x = TRUE)
  }

  detect_coord_names <- function(x) {
    lat_name <- names(x)[grep("lat", names(x))]
    lon_name <- names(x)[grep("lon", names(x))]
    c(lat_name, lon_name)
  }
  coord_names <- detect_coord_names(eval(as.symbol(gps)))
  #need to manually enter coord_names (lon_dd, lat_dd) for streamparse data
  create_basin_labels <- function(dt, coord_names) {
    projstr <- "+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
    latlonproj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    #### Added to make sure CRS were compatible ###
    fathombasins <- sf::read_sf(file.path(fdir, "DF_Basefile/fathom_basins_proj_updated/fathom_basins_proj_updated.shp"))
    cerpbasins <- sf::read_sf(file.path(fdir, "DF_Basefile/fbfs_zones.shp"))
    selectiongrid <- sf::read_sf(file.path(fdir, "DF_Basefile/testgrid9/testgrid9.shp"))

    fathombasins <- as(fathombasins, "Spatial")
    cerpbasins <- as(cerpbasins, "Spatial")
    selectiongrid <- as(selectiongrid, "Spatial")

    fathombasins <- sp::spTransform(fathombasins, projstr)
    cerpbasins <- sp::spTransform(cerpbasins, projstr)
    selectiongrid <- sp::spTransform(selectiongrid, projstr)
    #######################################################
    xy <- cbind(dt[, coord_names[2]], dt[, coord_names[1]])
    xy <- data.frame(xy)
    latcol <- which(colnames(dt)==coord_names[1])
    if (typeof(dt[,latcol])!="numeric") {
      dt[,latcol] <- as.numeric(dt[,latcol])
    }
    loncol <- which(colnames(dt)==coord_names[2])
    if (typeof(dt[,loncol])!="numeric") {
      dt[,loncol] <- as.numeric(dt[,loncol])
    }
    fulldataset <- coordinatize(dt, latname = coord_names[1],
                                lonname = coord_names[2])
    fulldataset.over <- sp::over(fulldataset, selectiongrid)
    fulldataset.over2 <- sp::over(fulldataset, fathombasins[,
                                                            1:2])
    fulldataset.over3 <- sp::over(fulldataset, cerpbasins[,
                                                          2])
    fulldataset.over <- cbind(data.frame(fulldataset), data.frame(fulldataset.over),
                              data.frame(fulldataset.over2), data.frame(fulldataset.over3))
    fulldataset.over$lon_dd <- xy[, 1]
    fulldataset.over$lat_dd <- xy[, 2]
    fulldataset.over[, names(fulldataset.over) != "NA."]
  }
  dt <- create_basin_labels(dt, coord_names)
  names(dt) <- tolower(names(dt))
  if (tofile == TRUE) {
    dtname <- file.path(fdir, .Platform$file.sep, "DF_FullDataSets",
                        .Platform$file.sep, substring(survey_days[1], 1,
                                                      6), "j.csv", fsep = "")
    if (!file.exists(dtname)) {
      write.csv(dt, dtname, row.names = FALSE)
    }
    else {
      stop("overwrite file?")
    }
  }
  options(warn = 0)
  dt
}


#'@name streamget
#'@title Retrieve previously cleaned full streaming datasets
#'@description Retrieve previously cleaned full streaming datasets
#'@param yearmon numeric date in yyyymm format
#'@param fdir character file path to local data directory
#'@param qa logical strip flagged data?
#'@export
#'@examples \dontrun{
#'yearmon <- 201212
#'dt <- streamget(yearmon)
#'}
streamget <- function(yearmon, qa = TRUE, fdir = getOption("fdir")){
  #Read QA file if it exists
  if(qa == TRUE && file.exists(file.path(fdir, "DF_FullDataSets", "QA datasets", paste(yearmon, "j_qa.csv", sep = "")))){
    dt <- read.csv(file.path(fdir, "DF_FullDataSets", "QA datasets", paste(yearmon, "j_qa.csv", sep = "")))
    if(!any(names(dt) == "chlext")){
      dt$chlext <- NA
    }
  } else {
    #Read original streamcleaned file if no QA
    fdir_fd <- file.path(fdir, "DF_FullDataSets")
    flist <- list.files(fdir_fd, include.dirs = T, full.names = T)
    flist <- flist[substring(basename(flist),1,6) == yearmon]
    dt <- read.csv(flist, stringsAsFactors = FALSE)
  }
  dt
}






