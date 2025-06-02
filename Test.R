# Test #

#Package creation resources:
#https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/

#R markdown info
#https://www.geeksforgeeks.org/how-to-create-an-r-markdown-document/


#Packages needed for package creation
library(devtools)
library(roxygen2)

#Create package documentation
setwd("\\\\ad.sfwmd.gov/dfsroot/userdata/kiallen/Docs/NewDataflowR")
#document()

#Install package
setwd("..")
#install("NewDataflowR")
library(NewDataflowR)

#Test functions

#Clean streaming data
dt <- streamclean(yearmon = 202503, gps = "exo", exommin = 12, c6mmin = 12, tofile = T)

#Load streamcleaned data
dt <- streamget(yearmon=202503, qa=T)

#Interpolate streaming data
streaminterp(streamget(yearmon = 202503, qa = TRUE), paramlist = c("salpsu"), 202503)

#Quick map of interpolated data
surfplot(rnge=c(202503), params=c("salpsu"))

#Clean grab data
g <- grabclean(yearmon = 202412, tofile = F)

#Load cleaned grab data
grabs <- grabget(rnge = c(201308))

#Create chl or cdom regressions
chlcoef(yearmon=202412, varlist=c("bgaperfu", "chlrfu", "fdomrfu", "c6chl", "bgapcrfu",
                                  "c6chlar", "phycoe", "phycoc", "c6cdom"), subgroup=c("Middle", "Monroe Lake", "Taylor River", "Pond 5", "Seven Palm Lake"),
        remove.flags=T, overwrite=F)

cdomcoef(yearmon=202412, varlist=c("chlrfu", "fdomrfu", "c6chl",
                                   "c6chlar", "c6cdom"), subgroup=c("Middle", "Monroe Lake", "Taylor River", "Seven Palm Lake"),
         remove.flags=T, overwrite=F)
cdomcoef(yearmon=202412, varlist=c("chlrfu", "fdomrfu", "c6chl",
                                   "c6chlar", "c6cdom"),
         remove.flags=T, overwrite=F)

#Interpolate maps of chl and cdom
chlmap(yearmon = 202412, subgroup=c("Middle", "Monroe Lake", "Taylor River", "Seven Palm Lake", "Pond 5"))

cdommap(yearmon = 202412, subgroup=c("Middle", "Monroe Lake", "Taylor River", "Seven Palm Lake", "Pond 5"))

#Quick map of interpolate chl
surfplot(rnge=c(202412), params=c("chlext"))



