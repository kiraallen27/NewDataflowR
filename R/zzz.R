.onAttach <- function(libname = find.package("NewDataflowR"), pkgname = "NewDataflowR"){

  localpath <- system.file("localpath", package = "NewDataflowR")
  fdir <- matrix(utils::read.delim(localpath, header = FALSE,
                                   stringsAsFactors = FALSE))[[1]][1]

  packageStartupMessage(paste("NewDataflowR Data Directory:",
                              fdir, "\n", "To change, modify:", localpath))

  options("fdir" = fdir)
}
