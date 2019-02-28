setwd("~/Dropbox/Univ/masterThesis/analysis/20180918_final/")

library(xlsx)
library(ggplot2)

data <- read.xlsx("measurements.xlsx", 1)

cellid = gsub("(.*)/cell_", "\\1/", as.character(data$NAME))
data$DATE = gsub("(\\d+)[AB]?/(\\d+)_\\d{4}.abf$", "\\1", cellid)
data$CELLNO = as.numeric(gsub("(.+)/(\\d+)_\\d{4}.abf$", "\\2", cellid))
data$CELLID = paste0(data$DATE, "_", data$CELLNO)

data$is.base = grepl("Base", data$PROTOCOL)
data$is.nbqx = grepl("NBQX", data$PROTOCOL)
data$is.quinpirole = grepl("Quinpirole", data$PROTOCOL)
data$is.sulpiride = grepl("Sulpiride", data$PROTOCOL)
data$is.bicuculline = grepl("Bicuculline", data$PROTOCOL)
data$is.baclofen = grepl("Baclofen", data$PROTOCOL)
data$is.sch23390 = grepl("SCH23390", data$PROTOCOL)

data = data[data$DATE != "20180903", ]
data = data[data$CELLID != "20180830_115", ]
data = data[data$CELLID != "20180803_71", ]
data = data[data$CELLID != "20180803_59", ]
data = data[data$CELLID != "20180911_142", ]
data = data[data$BASELINE.1.0.pA > -400, ]
# data = data[!(duplicated(paste(data$CELLID, data$PROTOCOL), fromLast = T) & !data$is.base), ]

# a = aggregate(data$is.quinpirole, list(CELLID = data$CELLID), sum)
# data = data[data$CELLID %in% a$CELLID[a$x == 1], ]
# a = aggregate(data$is.nbqx, list(CELLID = data$CELLID), sum)
# data = data[data$CELLID %in% a$CELLID[a$x == 1], ]

data$AMPLITUDE_REGION_00_SERIES_00_1.0.pA = abs(data$AMPLITUDE_REGION_00_SERIES_00_1.0.pA)
data$AMPLITUDE_REGION_01_SERIES_00_1.0.pA = abs(data$AMPLITUDE_REGION_01_SERIES_00_1.0.pA)

res = data.frame(CELLID = data$CELLID[!duplicated(data$CELLID)])

res$amplitude1 = data$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[!duplicated(data$CELLID)]
res$amplitude2 = data$AMPLITUDE_REGION_01_SERIES_00_1.0.pA[!duplicated(data$CELLID)]

res$charge = data$CHARGE_REGION_00_SERIES_00_1.0.s.pA[!duplicated(data$CELLID)]

res$risetime = as.numeric(as.character(data$RISE_TIME_REGION_00_SERIES_00_1.0.ms[!duplicated(data$CELLID)]))
res$npeaks = as.numeric(as.character(data$N_PEAKS_REGION_00_SERIES_00_1.0.dimensionless[!duplicated(data$CELLID)]))

res = res[res$risetime > 0 & res$charge > 0, ]

week = format(as.Date(gsub("(\\d+)_.+", "\\1", as.character(res$CELLID)), format = "%Y%m%d"), "%U")
res$CULTURE = factor(week)
res$CULTURE[res$CULTURE == 31] = 30

fields = c("nbqx", "sulpiride", "quinpirole", "bicuculline", "baclofen", "sch23390")
for (j in 1:length(fields)) {
  field = fields[j]
  res[field] = NA
  res[paste0(field, ".amplitude")] = NA
  res[paste0(field, ".ppf")] = NA
  res[paste0(field, ".relppf")] = NA
  for (i in 1:nrow(res)) {
    index = which(res$CELLID[i] == data$CELLID & data[paste0("is.", field)])
    if (length(index) > 0) {
      if (length(index) > 1) {
        index = index[length(index)]
      }
      if (data$CELLID[index+1] == res$CELLID[i] &&
          data$CELLID[index-1] == res$CELLID[i]) {
        reference = (data$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index-1] + 
                       data$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index+1]) / 2
        res[i, field] = data$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index] / reference
        
        
        reference0 = data$AMPLITUDE_REGION_01_SERIES_00_1.0.pA[index-1] / data$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index-1]
        reference1 = data$AMPLITUDE_REGION_01_SERIES_00_1.0.pA[index+1] / data$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index+1]
        reference = (reference0 + 
                       reference1) / 2
        res[i, paste0(field, ".relppf")] = 
          data$AMPLITUDE_REGION_01_SERIES_00_1.0.pA[index] / data$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index] / reference
      } else if (data$CELLID[index+1] == res$CELLID[i] ) {
        reference = data$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index+1]
        res[i, field] = data$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index] / reference
        
        reference = data$AMPLITUDE_REGION_01_SERIES_00_1.0.pA[index+1] / data$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index+1]
        res[i, paste0(field, ".relppf")] = 
          data$AMPLITUDE_REGION_01_SERIES_00_1.0.pA[index] / data$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index] / reference
      }else if (data$CELLID[index-1] == res$CELLID[i] ) {
        reference = data$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index-1]
        res[i, field] = data$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index] / reference
        
        reference = data$AMPLITUDE_REGION_01_SERIES_00_1.0.pA[index-1] / data$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index-1]
        res[i, paste0(field, ".relppf")] = 
          data$AMPLITUDE_REGION_01_SERIES_00_1.0.pA[index] / data$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index] / reference
      }
      res[i, paste0(field, ".amplitude")] = data$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index]
      # res[i, paste0(field, ".amplitude2")] = data$AMPLITUDE_REGION_01_SERIES_00_1.0.pA[index]
      res[i, paste0(field, ".ppf")] = data$AMPLITUDE_REGION_01_SERIES_00_1.0.pA[index] / res[i, paste0(field, ".amplitude")]
    }
  }
}

fields = c("nbqx", "sulpiride", "quinpirole", "bicuculline", "baclofen", "sch23390")
for (j in 1:length(fields)) {
  field = fields[j]
  res[paste0("c", field)] = NA
  for (i in 1:nrow(res)) {
    index = which(res$CELLID[i] == data$CELLID & data[paste0("is.", field)])
    if (length(index) > 0) {
      if (length(index) > 1) {
        index = index[length(index)]
      }
      if (data$CELLID[index+1] == res$CELLID[i] &&
          data$CELLID[index-1] == res$CELLID[i]) {
        reference = (data$CHARGE_REGION_00_SERIES_00_1.0.s.pA[index-1] + 
                       data$CHARGE_REGION_00_SERIES_00_1.0.s.pA[index+1]) / 2
        res[i, paste0("c", field)] = data$CHARGE_REGION_00_SERIES_00_1.0.s.pA[index] / reference
      } else if (data$CELLID[index+1] == res$CELLID[i] ) {
        reference = data$CHARGE_REGION_00_SERIES_00_1.0.s.pA[index+1]
        res[i, paste0("c", field)] = data$CHARGE_REGION_00_SERIES_00_1.0.s.pA[index] / reference
      }else if (data$CELLID[index-1] == res$CELLID[i] ) {
        reference = data$CHARGE_REGION_00_SERIES_00_1.0.s.pA[index-1]
        res[i, paste0("c", field)] = data$CHARGE_REGION_00_SERIES_00_1.0.s.pA[index] / reference
      }
    }
  }
}

res$ppf = NA
for (i in 1:nrow(res)) {
  first = data[!duplicated(data$CELLID) & data$is.base & res$CELLID[i] == data$CELLID, ]
  first = na.omit(first)
  if (nrow(first) > 0) {
    ppf = first$AMPLITUDE_REGION_01_SERIES_00_1.0.pA / first$AMPLITUDE_REGION_00_SERIES_00_1.0.pA
    res$ppf[i] = ppf
  }
}

log.trans.mean_cl_normal <- function(x, ...) {
  print("log.trans.mean_cl_normal")
  print(x)
  log10(mean_cl_normal(10 ^ x, ...))
}


silent = res$amplitude1 < 10 
glut = res$nbqx < 0.1 & !silent
res$glut = factor(c("Other", "Glu")[(glut&!silent)+1], levels = c("Glu", "Other"))

fix.names <- function(dat) {
  dat$variable = as.character(dat$variable)
  dat$variable[dat$variable == "nbqx"] = "NBQX"
  dat$variable[dat$variable == "mix"] = "NBQX\nSCH23390\nSulpiride"
  dat$variable[dat$variable == "sulp"] = "Sulpiride"
  dat$variable[dat$variable == "sch"] = "SCH23390"
  dat$variable[dat$variable == "quin"] = "Quinpirole"
  dat$variable[dat$variable == "baclofen"] = "Baclofen"
  dat$variable = factor(dat$variable, levels = c("SCH23390", "Sulpiride", "NBQX\nSCH23390\nSulpiride", "Quinpirole", "Baclofen"))
  return(dat)
}


# hippocampus <- read.xlsx("hippocampus.xlsx", 1)
# 
# cellid = gsub("(.*)/cell_", "\\1/", as.character(hippocampus$NAME))
# hippocampus$DATE = gsub("(\\d+)[AB]?/(\\d+)_\\d{4}.abf$", "\\1", cellid)
# hippocampus$CELLNO = as.numeric(gsub("(.+)/(\\d+)_\\d{4}.abf$", "\\2", cellid))
# hippocampus$CELLID = paste0(hippocampus$DATE, "_", hippocampus$CELLNO)
# 
# hippocampus = hippocampus[hippocampus$BASELINE.1.0.pA > -400, ]
# 
# hippocampus$is.base = grepl("Base", hippocampus$PROTOCOL)
# hippocampus$is.nbqx = grepl("NBQX", hippocampus$PROTOCOL)
# 
# hres = data.frame(CELLID = hippocampus$CELLID[!duplicated(hippocampus$CELLID)])
# 
# hres$amplitude1 = as.numeric(as.character(hippocampus$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[!duplicated(hippocampus$CELLID)]))
# hres$amplitude2 = as.numeric(as.character(hippocampus$AMPLITUDE_REGION_01_SERIES_00_1.0.pA[!duplicated(hippocampus$CELLID)]))
# 
# hres$risetime = as.numeric(as.character(hippocampus$RISE_TIME_REGION_00_SERIES_00_1.0.ms[!duplicated(hippocampus$CELLID)]))
# hres$npeaks = as.numeric(as.character(hippocampus$N_PEAKS_REGION_00_SERIES_00_1.0.dimensionless[!duplicated(hippocampus$CELLID)]))
# 
# field = "nbqx"
# hres[field] = NA
# for (i in 1:nrow(hres)) {
#   index = which(hres$CELLID[i] == hippocampus$CELLID & hippocampus[paste0("is.", field)])
#   if (length(index) > 0) {
#     if (hippocampus$CELLID[index+1] == hres$CELLID[i] && hippocampus$CELLID[index-2] == hres$CELLID[i]) {
#       reference = (hippocampus$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index-2] + 
#                      hippocampus$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index+1]) / 2
#       hres[i, field] = hippocampus$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index] / reference
#     } else if (hippocampus$CELLID[index-2] == hres$CELLID[i]) {
#       reference = hippocampus$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index-2] 
#       hres[i, field] = hippocampus$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index] / reference
#     }
#   }
# }
# 
# hres$ppf = NA
# for (i in 1:nrow(res)) {
#   first = hippocampus[!duplicated(hippocampus$CELLID) & hippocampus$is.base & hres$CELLID[i] == hippocampus$CELLID, ]
#   first = na.omit(first)
#   if (nrow(first) > 0) {
#     ppf = first$AMPLITUDE_REGION_01_SERIES_00_1.0.pA / first$AMPLITUDE_REGION_00_SERIES_00_1.0.pA
#     hres$ppf[i] = ppf
#   }
# }
# 
# # silent = hres$amplitude1 < 10
# # hres = hres[hres$nbqx > 0,]
# # hres = hres[hres$risetime > 0,]
# 
# week = format(as.Date(gsub("(\\d+)_.+", "\\1", as.character(hres$CELLID)), format = "%Y%m%d"), "%U")
# hres$CULTURE = factor(week)
# hres$CULTURE[hres$CULTURE == 31] = 30

# source("figure1.R")
# source("figure2.R")
# source("figure3.R")
# source("figure4.R")
# source("figure5.R")
# source("figure6.R")
# source("figure7.R")
