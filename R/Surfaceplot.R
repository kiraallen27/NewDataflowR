#'@name surfplot
#'@title Plotting Interpolated Surfaces
#'@description accounts for corner case parameter spellings, variable specific contour breaks are (need to be) defined
#'@author Jemma Stachelek and Kira Allen
#'@param rnge numeric string of no more than two dates in yyyymm format
#'@param params character. string of parameter fields to plot
#'@param fdir character file path to local data directory
#'@param yext numeric length 2 y extent
#'@param xext numeric length 2 x extent
#'@return output plots to plot window
#'@importFrom sf st_read
#'@importFrom rasterVis levelplot
#'@importFrom latticeExtra layer
#'@importFrom sp spplot sp.polygons
#'@importFrom raster raster stack reclassify calc writeRaster
#'@export
#'@examples \dontrun{
#'surfplot(rnge = c(200707), params = c("cdom"))
#'surfplot(201513, c("ph", "c6turbidity", "c6chl", "c6cdom"),
#' yext = c(2786102, 2797996), xext = c(557217, 567415))
#'}
surfplot <- function(rnge, params, fdir = getOption("fdir"), yext = c(2772256, 2798000), xext = c(518000.2, 566000)){

  if(length(rnge) == 1){
    rnge <- c(rnge, rnge)
  }
  namesalias <- read.table(text = "
                       chlorophyll.a c6chl
                       c6chla c6chl
                       ")
  #define breaks
  brks <- read.table(text = "
        sal list(seq(0,40,2))
        salinity.pss list(seq(0,40,2))
        salpsu list(seq(0,40,2))
        c6chl list(seq(50,200,10))
        chlext list(seq(0,5,0.5),seq(10,30,5))
        chlext_hi list(seq(0,5,0.5),seq(10,30,5))
        chlext_low list(seq(0,5,0.5),seq(10,30,5))
        temp list(seq(14,36,2))
        c6temp list(seq(14,36,2))
        c6cdom list(seq(80,360,40))
        ph list(seq(7.3,8.1,0.1))
        c6turbidity list(seq(0,45,5))")

  dirlist <- list.dirs(file.path(fdir, "DF_Surfaces"), recursive = F)

  minrnge <- min(which(substring(basename(dirlist), 1, 6) >= rnge[1]))
  maxrnge <- max(which(substring(basename(dirlist), 1, 6) <= rnge[2]))
  rlist <- list.files(dirlist[minrnge:maxrnge], full.names = T, include.dirs = T, pattern = "\\.tif$")
  plist <- tolower(sub("[.][^.]*$", "", basename(rlist)))

  for(n in 1:length(plist)){
    if(any(plist[n] == namesalias[,1])){
      plist[n] <- as.character(namesalias[which(plist[n] == namesalias[,1]), 2])
      #print(names(dt)[n])
    }
  }

  rlist <- rlist[which(!is.na(
    match(plist, params)
  ))]
  plist <- plist[which(!is.na(match(plist, params)))]

  for(i in 1:length(rlist)){
    #New
    r <- raster::raster(rlist[i])
    ymin <- r@extent[3]
    ymax <- r@extent[4]
    xmin <- r@extent[1]
    xmax <- r@extent[2]

    # Plot
    print(
      rasterVis::levelplot(
        r,
        ylim = c(ymin, ymax),
        xlim = c(xmin, xmax),
        par.settings = rasterVis::BuRdTheme(),
        margin = FALSE,
        auto.key = FALSE,
        scales = list(draw = FALSE),
        main = paste(
          as.character(plist[i]),
          unlist(strsplit(rlist[i], "/"))[length(unlist(strsplit(rlist[i], "/"))) - 1]
        )
      )
      + latticeExtra::layer({
        sp::SpatialPolygonsRescale(
          sp::layout.north.arrow(),
          offset = c(563000, 2775000),
          scale = 4400
        )
      })
      + latticeExtra::layer(sp::sp.polygons(as(sf::st_read(dsn = file.path(getOption("fdir"), "DF_Basefile", "FBcoast_big", "FBcoast_big.shp"), quiet = TRUE), "Spatial"),
                                            fill = "seagreen", alpha = 0.6))
    )

  }
}
