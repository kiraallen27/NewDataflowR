# Test #

#Package creation resources:
#https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/

#R markdown info
#https://www.geeksforgeeks.org/how-to-create-an-r-markdown-document/

n <- 1
pkginstall <- function(package) {
  while(system.file(package=package)=="") {
    install.packages(package)
    if(n > 10) {
      break
    }
    n <- n+1
  }
}

package <- "xfun"
pkginstall(package)


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
dt <- streamclean(yearmon = 202602, gps = "exo", exommin = 12, c6mmin = 12, tofile = T)

dt <- streamclean(yearmon=200704, gps="df", dfmmin=7, tofile=T)

#Load streamcleaned data
dt <- streamget(yearmon=202602, qa=F)

#Interpolate streaming data
streaminterp(streamget(yearmon = 202602, qa = FALSE), paramlist = c("temp.c","sal.psu", "chlorophyll.rfu"), 202602)

#Quick map of interpolated data
surfplot(rnge=c(202602), params=c("chlorophyll.rfu"))

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

#Test QAdash
QAdash(yearmon=202510, original.data=F)


#### Install package from GitHub #########

devtools::install_github("kiraallen27/NewDataflowR", dependencies=T)
library(NewDataflowR)





## Old QA flags functions ##

#Test QA functions
d <- QAflags(yearmon=202602, param="turbidity.fnu", bad_min=0, bad_max=500, sus_min=0.01, sus_max=200, step_threshold=10,
             plots=c("flag map", "step flag points"), original.data=T)
data <- d[[1]]
d[[2]]

df <- replace_data(df=data, yearmon=202602, value_var="tal.pc.rfu", min_value=0, replace_val=0, save=T)
df <- replace_data(df=data, yearmon=202602, value_var="tal.pc.rfu", max_value=6, replace_val=NA, save=T)

d <- QAflags(yearmon=202510, param="tal.pc.rfu", bad_min=0, bad_max=10,   sus_min=0.01, sus_max=4, step_threshold=2, plots=c("flag points", "step flag map"), original.data=T)
d[[2]]

d <- QAflags(yearmon=202510, param="tal.pc.rfu", bad_min=0, bad_max=10,   sus_min=0.01, sus_max=4, step_threshold=2, plots=c("flag points", "flag map"), original.data=TRUE)
#View data output
data <- d[[1]]
#View data plots
d[[2]]

#Removes values of tal.pc.rfu that are flagged as “bad” and when #step_flag=TRUE
df <- replace_data(df=data, #data is output from QAflags
                   yearmon=202510, value_var="tal.pc.rfu", flag="Bad", step_flag=TRUE, replace_val=NA, logic="and", save=TRUE)

#Replaces negative tal.pc.rfu values with 0
df <- replace_data(df=data, #data is output from QAflags
                   yearmon=202510, value_var="tal.pc.rfu", min_value=0, replace_val=0, save=TRUE)

#Removes values of tal.pc.rfu that occur between a range of datetimes
df <- replace_data(df=data, #data is output from QAflags
                   yearmon=202510, value_var="tal.pc.rfu", start_date="2025-10-21 12:57:42", end_date="2025-10-21 13:13:13", save=TRUE)

d <- QAflags(yearmon=202510, param="tal.pc.rfu", bad_min=0, bad_max=10,   sus_min=0.01, sus_max=4, step_threshold=2, plots=c("flag points", "flag map"), original.data=FALSE)
data <- d[[1]]
d[[2]]
