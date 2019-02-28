

fiveap <- read.xlsx("5ap.xlsx", 1)

cellid = gsub("(.*)/cell_", "\\1/", as.character(fiveap$NAME))
fiveap$DATE = gsub("(\\d+)[AB]?/(\\d+)_\\d{4}.abf$", "\\1", cellid)
fiveap$CELLNO = as.numeric(gsub("(.+)/(\\d+)_\\d{4}.abf$", "\\2", cellid))
fiveap$CELLID = paste0(fiveap$DATE, "_", fiveap$CELLNO)

fiveap = fiveap[fiveap$CELLID %in% res$CELLID,]
fiveap$CULTURE = res$CULTURE[match(fiveap$CELLID, res$CELLID)]

fiveap$amplitude1 = res$amplitude1[match(fiveap$CELLID, res$CELLID)]
fiveap$amplitude2 = res$amplitude2[match(fiveap$CELLID, res$CELLID)]
fiveap$nbqx = res$nbqx[match(fiveap$CELLID, res$CELLID)]

fiveap = fiveap[fiveap$BASELINE.1.0.pA > -400, ]

# fiveap$sucrose = sucrose$CHARGE_SERIES_00_1.0.s.pA[match(fiveap$cellid, sucrose$cellid)]

# fiveap = fiveap[fiveap$amplitude1 > 10 & fiveap$amplitude2 > 10 & fiveap$nbqx < 0.1, ]
fiveap$PROTOCOL = as.character(fiveap$PROTOCOL)

res$mix = NA

unique(fiveap$CELLID) %in% res$CELLID

field = "mix"
for (i in 1:nrow(fiveap)) {
  name_ = gsub("(.+)\\d{4}.abf", "\\1", as.character(fiveap$NAME)[i])
  index_ = as.numeric(gsub("(.+)(\\d{4}).abf", "\\2", as.character(fiveap$NAME)[i]))
  index0 = which(data$NAME == sprintf("%s%04d.abf", name_, index_+1))
  index1 = which(data$NAME == sprintf("%s%04d.abf", name_, index_-1))
  index2 = which(data$NAME == sprintf("%s%04d.abf", name_, index_-2))
  index3 = which(data$NAME == sprintf("%s%04d.abf", name_, index_-3))
  index4 = which(data$NAME == sprintf("%s%04d.abf", name_, index_-4))
  
  index5 = which(res$CELLID == fiveap$CELLID[i])
  
  logical0 = data$CELLID[index0] == fiveap$CELLID[i]
  logical1 = data$CELLID[index1] == fiveap$CELLID[i]
  logical2 = data$CELLID[index2] == fiveap$CELLID[i]
  logical3 = data$CELLID[index3] == fiveap$CELLID[i]
  logical4 = data$CELLID[index4] == fiveap$CELLID[i]
  
  if (length(logical0) == 0) logical0 = F
  if (length(logical1) == 0) logical1 = F
  if (length(logical2) == 0) logical2 = F
  if (length(logical3) == 0) logical3 = F
  if (length(logical4) == 0) logical4 = F
  
  if (logical0 && logical1) {
    reference = (data$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index0] + 
                   data$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index1]) / 2
  } else if (logical0) {
    reference = data$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index0]
  } else if (logical1) {
    reference = data$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index1]
  } else if (logical2) {
    reference = data$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index2]
  } else if (logical3) {
    reference = data$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index3]
  } else if (logical4) {
    reference = data$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[index4]
  } else if (is.na(res[index5, field])) {
    print(i)
    print(as.character(fiveap$NAME)[i])
  }
  res[index5, field] = fiveap$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[i] / reference
  res[index5, paste0(field, ".amplitude")] = fiveap$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[i] 
  res[index5, paste0(field, ".ppf")] = fiveap$AMPLITUDE_REGION_01_SERIES_00_1.0.pA[i] / fiveap$AMPLITUDE_REGION_00_SERIES_00_1.0.pA[i]
}


field = "cmix"
for (i in 1:nrow(fiveap)) {
  name_ = gsub("(.+)\\d{4}.abf", "\\1", as.character(fiveap$NAME)[i])
  index_ = as.numeric(gsub("(.+)(\\d{4}).abf", "\\2", as.character(fiveap$NAME)[i]))
  index0 = which(data$NAME == sprintf("%s%04d.abf", name_, index_+1))
  index1 = which(data$NAME == sprintf("%s%04d.abf", name_, index_-1))
  index2 = which(data$NAME == sprintf("%s%04d.abf", name_, index_-2))
  index3 = which(data$NAME == sprintf("%s%04d.abf", name_, index_-3))
  index4 = which(data$NAME == sprintf("%s%04d.abf", name_, index_-4))
  
  index5 = which(res$CELLID == fiveap$CELLID[i])
  
  logical0 = data$CELLID[index0] == fiveap$CELLID[i]
  logical1 = data$CELLID[index1] == fiveap$CELLID[i]
  logical2 = data$CELLID[index2] == fiveap$CELLID[i]
  logical3 = data$CELLID[index3] == fiveap$CELLID[i]
  logical4 = data$CELLID[index4] == fiveap$CELLID[i]
  
  if (length(logical0) == 0) logical0 = F
  if (length(logical1) == 0) logical1 = F
  if (length(logical2) == 0) logical2 = F
  if (length(logical3) == 0) logical3 = F
  if (length(logical4) == 0) logical4 = F
  
  if (logical0 && logical1) {
    reference = (data$CHARGE_REGION_00_SERIES_00_1.0.s.pA[index0] + 
                   data$CHARGE_REGION_00_SERIES_00_1.0.s.pA[index1]) / 2
    res[index5, field] = fiveap$CHARGE_REGION_00_SERIES_00_1.0.s.pA[i] / reference
  } else if (logical0) {
    reference = data$CHARGE_REGION_00_SERIES_00_1.0.s.pA[index0]
    res[index5, field] = fiveap$CHARGE_REGION_00_SERIES_00_1.0.s.pA[i] / reference
  } else if (logical1) {
    reference = data$CHARGE_REGION_00_SERIES_00_1.0.s.pA[index1]
    res[index5, field] = fiveap$CHARGE_REGION_00_SERIES_00_1.0.s.pA[i] / reference
  } else if (logical2) {
    reference = data$CHARGE_REGION_00_SERIES_00_1.0.s.pA[index2]
    res[index5, field] = fiveap$CHARGE_REGION_00_SERIES_00_1.0.s.pA[i] / reference
  } else if (logical3) {
    reference = data$CHARGE_REGION_00_SERIES_00_1.0.s.pA[index3]
    res[index5, field] = fiveap$CHARGE_REGION_00_SERIES_00_1.0.s.pA[i] / reference
  } else if (logical4) {
    reference = data$CHARGE_REGION_00_SERIES_00_1.0.s.pA[index4]
    res[index5, field] = fiveap$CHARGE_REGION_00_SERIES_00_1.0.s.pA[i] / reference
  } else if (is.na(res[index5, field])) {
    print(i)
    print(as.character(fiveap$NAME)[i])
  }
}





# in_field = "AMPLITUDE"
# traceScr = fiveap[grepl("20Hz_5AP_NBQX", fiveap$PROTOCOL), grepl(in_field, names(fiveap))]
# traceScr[] <- lapply(traceScr, function(x) {
#   as.numeric(as.character(x)) 
# })
# sapply(traceScr, class)
# traceScr = traceScr[traceScr[,1] > 0,]
# 
# traceKd = fiveap[grepl("20Hz_5AP_Mix", fiveap$PROTOCOL), grepl(in_field, names(fiveap))]
# traceKd[] <- lapply(traceKd, function(x) {
#   as.numeric(as.character(x)) 
# })
# sapply(traceKd, class)
# traceKd = traceKd[traceKd[,1] > 0,]
# 
# normTraceScr = traceScr
# normTraceKd = traceKd
# # normTraceScr = sweep(traceScr, 1, apply(traceScr, 1, min), `-`)
# # normTraceKd = sweep(traceKd, 1, apply(traceKd, 1, min), `-`)
# 
# # normTraceScrMax = apply(normTraceScr, 1, max)
# # normTraceScrMax[normTraceScrMax < 1e-9] = 1
# # normTraceScr = sweep(normTraceScr, 1, normTraceScr[,1], `/`)
# # normTraceKdMax = apply(normTraceKd, 1, max)
# # normTraceKdMax[normTraceKdMax < 1e-9] = 1
# # normTraceKd = sweep(normTraceKd, 1, normTraceKd[,1], `/`)
# # 
# # if (in_field == "RISE_TIME") {
# #   normTraceScr = sweep(traceScr, 1, 1, `/`)
# #   normTraceKd = sweep(traceKd, 1, 1, `/`)
# # }
# 
# meanTraceScr = apply(
#   normTraceScr,
#   2, mean, na.rm = T)
# meanTraceKd = apply(
#   normTraceKd,
#   2, mean, na.rm = T)
# 
# seTraceScr = apply(
#   normTraceScr,
#   2, sd, na.rm = T) / sqrt(nrow(traceScr)/2)
# seTraceKd = apply(
#   normTraceKd,
#   2, sd, na.rm = T) / sqrt(nrow(traceKd))
# 
# ylab = "Amplitude [pA]"
# 
# require(ggplot2)
# int = 1
# dfPlot = data.frame(x = rep(int * 1:length(meanTraceScr), 2),
#                     y = c(meanTraceScr, meanTraceKd),
#                     se = c(seTraceScr, seTraceKd),
#                     Group = c(rep("Control", length(meanTraceScr)), rep("-DA", length(meanTraceScr))))
# p <- ggplot(dfPlot, aes(x = x, y = y, group = Group, color = Group))
# p <- p + geom_line()
# p <- p + geom_ribbon(aes(ymin = y - se, ymax = y + se, fill = Group),
#                      color = "transparent", alpha = 0.15, show.legend = F)
# p <- p + labs(fill = "", x = "Stimulus #", y = ylab)
# p <- p + theme_classic()
# 
# pdf(paste0("figure7b", ".pdf"), width = 4, height = 3)
# p
# dev.off()
# 
# 
# in_field = "AMPLITUDE"
# traceScr = fiveap[grepl("50Hz_5AP_NBQX", fiveap$PROTOCOL), grepl(in_field, names(fiveap))]
# traceScr[] <- lapply(traceScr, function(x) {
#   as.numeric(as.character(x)) 
# })
# sapply(traceScr, class)
# 
# traceKd = fiveap[grepl("50Hz_5AP_Mix", fiveap$PROTOCOL), grepl(in_field, names(fiveap))]
# traceKd[] <- lapply(traceKd, function(x) {
#   as.numeric(as.character(x)) 
# })
# sapply(traceKd, class)
# 
# normTraceScr = traceScr
# normTraceKd = traceKd
# # normTraceScr = sweep(traceScr, 1, apply(traceScr, 1, min), `-`)
# # normTraceKd = sweep(traceKd, 1, apply(traceKd, 1, min), `-`)
# 
# # normTraceScrMax = apply(normTraceScr, 1, max)
# # normTraceScrMax[normTraceScrMax < 1e-9] = 1
# # normTraceScr = sweep(normTraceScr, 1, normTraceScr[,1], `/`)
# # normTraceKdMax = apply(normTraceKd, 1, max)
# # normTraceKdMax[normTraceKdMax < 1e-9] = 1
# # normTraceKd = sweep(normTraceKd, 1, normTraceKd[,1], `/`)
# # 
# # if (in_field == "RISE_TIME") {
# #   normTraceScr = sweep(traceScr, 1, 1, `/`)
# #   normTraceKd = sweep(traceKd, 1, 1, `/`)
# # }
# 
# meanTraceScr = apply(
#   normTraceScr,
#   2, mean, na.rm = T)
# meanTraceKd = apply(
#   normTraceKd,
#   2, mean, na.rm = T)
# 
# seTraceScr = apply(
#   normTraceScr,
#   2, sd, na.rm = T) / sqrt(nrow(traceScr)/2)
# seTraceKd = apply(
#   normTraceKd,
#   2, sd, na.rm = T) / sqrt(nrow(traceKd))
# 
# ylab = "Amplitude [pA]"
# 
# require(ggplot2)
# int = 1
# dfPlot = data.frame(x = rep(int * 1:length(meanTraceScr), 2),
#                     y = c(meanTraceScr, meanTraceKd),
#                     se = c(seTraceScr, seTraceKd),
#                     Group = c(rep("Control", length(meanTraceScr)), 
#                               rep("-DA", length(meanTraceKd))))
# p <- ggplot(dfPlot, aes(x = x, y = y, group = Group, color = Group))
# p <- p + geom_line()
# p <- p + geom_ribbon(aes(ymin = y - se, ymax = y + se, fill = Group),
#                      color = "transparent", alpha = 0.15, show.legend = F)
# p <- p + labs(fill = "", x = "Stimulus #", y = ylab)
# p <- p + theme_classic()
# 
# pdf(paste0("figure7c", ".pdf"), width = 4, height = 3)
# p
# dev.off()

