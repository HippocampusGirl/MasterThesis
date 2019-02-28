setwd("~/Dropbox/Univ/masterThesis/analysis/20181106/")

silent = res$amplitude1 < 10 
glut = res$nbqx < 0.1 & !silent
nonglut = !silent & !glut

cc = read.csv("cumtrapz.csv", sep = " ")

cellid = gsub("(.*)/cell_", "\\1/", as.character(cc$FILENAME))
cc$DATE = gsub("(\\d+)[AB]?/(\\d+)_\\d{4}.abf$", "\\1", cellid)
cc$CELLNO = as.numeric(gsub("(.+)/(\\d+)_\\d{4}.abf$", "\\2", cellid))
cc$RECNO = as.numeric(gsub("(.+)/(\\d+)_(\\d{4}).abf$", "\\3", cellid))
cc$CELLID = paste0(cc$DATE, "_", cc$CELLNO)

cc = cc[cc$CELLID %in% res$CELLID[glut],]
# cc = cc[cc$DATA_1000 > -100,]

# hcc = read.csv("hcumtrapz.csv", sep = " ")
# 
# cellid = gsub("(.*)/cell_", "\\1/", as.character(hcc$FILENAME))
# hcc$DATE = gsub("(\\d+)[AB]?/(\\d+)_\\d{4}.abf$", "\\1", cellid)
# hcc$CELLNO = as.numeric(gsub("(.+)/(\\d+)_\\d{4}.abf$", "\\2", cellid))
# hcc$RECNO = as.numeric(gsub("(.+)/(\\d+)_(\\d{4}).abf$", "\\3", cellid))
# hcc$CELLID = paste0(hcc$DATE, "_", hcc$CELLNO)
# 
# hcc = hcc[hcc$CELLID %in% hres$CELLID[hres$amplitude1 > 10 & hres$nbqx < 0.1],]
# hcc = hcc[hcc$DATA_1000 > -100,]

get_mean_sd = function(data) {
  in_field = "DATA"
  traceScr = data[, grepl(in_field, names(data))][,1:nt]
  traceScr[] <- lapply(traceScr, function(x) {
    as.numeric(as.character(x)) 
  })
  # sapply(traceScr, class)
  
  traceScr = -traceScr
  traceScr[traceScr < 0] = 0
  
  # normTraceScr = sweep(traceScr, 1, apply(traceScr, 1, min), `-`)
  # normTraceKd = sweep(traceKd, 1, apply(traceKd, 1, min), `-`)
  # # apply(normTraceKd, 1, min)
  # 
  # normTraceScrMax = apply(normTraceScr, 1, max)
  # normTraceScrMax[normTraceScrMax < 1e-9] = 1
  # normTraceScr = sweep(normTraceScr, 1, normTraceScrMax, `/`)
  
  meanTraceScr = apply(
    traceScr,
    2, mean, na.rm = T)
  
  seTraceScr = apply(
    traceScr,
    2, sd, na.rm = T) / sqrt(nrow(traceScr))
  return(list(mean = meanTraceScr, sd = seTraceScr))
}

# ccbase = cc[grepl("Base", cc$PROTOCOL) & !(duplicated(paste(cc$CELLID, cc$PROTOCOL))),]
# plot((50:150) / 1000, -ccbase[8,][, grepl(in_field, names(ccbase))][,50:150])

nt = 1000

cc$CULTURE = res$CULTURE[match(cc$CELLID, res$CELLID)]

ccBase = cc[grepl("Base", cc$PROTOCOL) & !(duplicated(paste(cc$CELLID, cc$PROTOCOL), fromLast = F)),]
dataBase = get_mean_sd(ccBase)
dataNBQX = get_mean_sd(cc[grepl("NBQX", cc$PROTOCOL),])
dataSulpiride = get_mean_sd(cc[grepl("Sulpiride", cc$PROTOCOL),])
dataQuinpirole = get_mean_sd(cc[grepl("Quinpirole", cc$PROTOCOL),])
dataSCH23390 = get_mean_sd(cc[grepl("SCH23390", cc$PROTOCOL),])

ttext = sprintf("%d/%d",nrow(ccBase),length(unique(na.omit(ccBase$CULTURE))))

# dataMix = get_mean_sd(cc[grepl("Mix", cc$PROTOCOL),])

# hdataBase = get_mean_sd(hcc[grepl("Base", hcc$PROTOCOL) & !(duplicated(paste(hcc$CELLID, hcc$PROTOCOL))),])
# hdataNBQX = get_mean_sd(hcc[grepl("NBQX", hcc$PROTOCOL),])

ylab = "Cumulative charge [pC]"

# require(ggplot2)
# int = 1
# dfPlot = data.frame(x = rep(int * 1:nt, 7) / 1000,
#                     y = c(dataBase$mean, dataNBQX$mean, dataSulpiride$mean, dataQuinpirole$mean, dataSCH23390$mean),
#                     se = c(dataBase$sd, dataNBQX$sd, dataSulpiride$sd, dataQuinpirole$sd, dataSCH23390$sd,
#                            hdataBase$sd, hdataNBQX$sd),
#                     Type = c(rep("DA",nt), rep("DA NBQX", nt), 
#                              rep("DA Sulpiride", nt), rep("DA Quinpirole", nt), rep("DA SCH23390", nt),
#                             rep("HC", nt), rep("HC NBQX", nt)),
#                     Group = factor(c(rep("DA",nt*5),rep("HC",nt*2))))
# p <- ggplot(dfPlot, aes(x = x, y = y, color = Type))
# p <- p + geom_line()
# p <- p + geom_ribbon(aes(ymin = y - se, ymax = y + se, fill = Type),
#                      color = "transparent", alpha = 0.15, show.legend = F)
# p <- p + labs(fill = "", x = "Time [s]", y = ylab)
# p <- p + scale_color_manual(values = c("turquoise", "blue", "darkslateblue", "lightblue", "purple", "gray30", "gray50"))
# p <- p + scale_fill_manual(values =c("turquoise", "blue", "darkslateblue", "lightblue", "purple", "gray30", "gray50"))
# p <- p + theme_classic()
# p <- p + facet_grid(. ~ Group)
# 
# pdf(paste0("figure8a", ".pdf"), width = 6, height = 3.5)
# p
# dev.off()


# int = 1
# dfPlot = data.frame(x = rep(int * (1:nt) / 1000, 5),
#                     y = c(dataBase$mean, dataNBQX$mean, dataSCH23390$mean, dataSulpiride$mean, dataQuinpirole$mean),
#                     se = c(dataBase$sd, dataNBQX$sd, dataSCH23390$sd, dataSulpiride$sd, dataQuinpirole$sd),
#                     Type = factor(c(rep("Control",nt), rep("NBQX", nt), 
#                              rep("SCH23390", nt), rep("Sulpiride", nt), rep("Quinpirole", nt)), 
#                              levels = c("Control", "NBQX", "SCH23390", "Sulpiride", "Quinpirole", "Baclofen")))
# p <- ggplot(dfPlot, aes(x = x, y = y, color = Type))
# p <- p + geom_line(show.legend = F)
# p <- p + geom_ribbon(aes(ymin = y - se, ymax = y + se, fill = Type),
#                      color = "transparent", alpha = 0.15, show.legend = F)
# p <- p + geom_text(data = subset(dfPlot, x == 1), aes(label = Type, colour = Type, x = x, y = y), hjust = -.1, show.legend = F) 
# p <- p + labs(fill = "", x = "Time [s]", y = ylab)
# p <- p + scale_color_manual(values = c("turquoise", "blue", "darkslateblue", "lightblue", "purple", "gray30", "gray50"))
# p <- p + scale_fill_manual(values =c("turquoise", "blue", "darkslateblue", "lightblue", "purple", "gray30", "gray50"))
# p <- p + theme_classic()
# p <- p + theme(plot.margin = unit(c(1,3.5,1,1), "lines")) 
# p <- p + coord_cartesian(clip = "off") 
# 
# pdf(paste0("figure8a_", ".pdf"), width = 4, height = 3.5)
# p
# dev.off()


nt = 1000
ccBase = cc[grepl("Base", cc$PROTOCOL) & !(duplicated(paste(cc$CELLID, cc$PROTOCOL), fromLast = F)),]
dataBase = get_mean_sd(ccBase)

ttext = sprintf("%d/%d",nrow(ccBase),length(unique(na.omit(ccBase$CULTURE))))

int = 1
dfPlot = data.frame(x = rep(int * (1:nt) / 1000, 1),
                    y = c(dataBase$mean),
                    se = c(dataBase$sd),
                    Type = factor(c(rep("Control",nt)), 
                                  levels = c("Control")))

dftext = data.frame(Type = factor(c("DA")), label = c(ttext), stringsAsFactors = F)
dftext$x = mean(dfPlot$x)
dftext$y = max(dfPlot$y + dfPlot$se) + diff(range(dfPlot$y + dfPlot$se)) / 20

p <- ggplot(dfPlot, aes(x = x, y = y, color = Type))
p <- p + geom_line(show.legend = F)
p <- p + geom_ribbon(aes(ymin = y - se, ymax = y + se, fill = Type),
                     color = "transparent", alpha = 0.15, show.legend = F)
p <- p + geom_label(aes(label = label), 
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, colour = "black", 
                    show.legend = F) + expand_limits(y = dftext$value[1])
p <- p + labs(fill = "", x = "Time [s]", y = ylab)
p <- p + scale_color_manual(values = c("turquoise", "blue", "darkslateblue", "lightblue", "purple", "gray30", "gray50"))
p <- p + scale_fill_manual(values =c("turquoise", "blue", "darkslateblue", "lightblue", "purple", "gray30", "gray50"))
p <- p + theme_classic()
# p <- p + theme(plot.margin = unit(c(1,3.5,1,1), "lines")) 
# p <- p + coord_cartesian(clip = "off") 

pdf(paste0("figure8a", ".pdf"), width = 2, height = 3.5)
p
dev.off()


nt = 250
ccBase = cc[grepl("Base", cc$PROTOCOL) & !(duplicated(paste(cc$CELLID, cc$PROTOCOL), fromLast = F)),]
dataBase = get_mean_sd(ccBase)

ttext = sprintf("%d/%d",nrow(ccBase),length(unique(na.omit(ccBase$CULTURE))))

int = 1
dfPlot = data.frame(x = rep(int * (1:nt), 1),
                    y = c(dataBase$mean),
                    se = c(dataBase$sd),
                    Type = factor(c(rep("Control",nt)), 
                                  levels = c("Control")))

dftext = data.frame(Type = factor(c("DA")), label = c(ttext), stringsAsFactors = F)
dftext$x = mean(dfPlot$x)
dftext$y = max(dfPlot$y + dfPlot$se) + diff(range(dfPlot$y + dfPlot$se)) / 20

p <- ggplot(dfPlot, aes(x = x, y = y, color = Type))
p <- p + geom_line(show.legend = F)
p <- p + geom_ribbon(aes(ymin = y - se, ymax = y + se, fill = Type),
                     color = "transparent", alpha = 0.15, show.legend = F)
p <- p + geom_label(aes(label = label), 
                    data = dftext, stat = "identity", alpha = 0.5, colour = "black", 
                    show.legend = F) + expand_limits(y = dftext$value[1])
p <- p + labs(fill = "", x = "Time [ms]", y = ylab)
p <- p + scale_color_manual(values = c("turquoise", "blue", "darkslateblue", "lightblue", "purple", "gray30", "gray50"))
p <- p + scale_fill_manual(values =c("turquoise", "blue", "darkslateblue", "lightblue", "purple", "gray30", "gray50"))
p <- p + theme_classic()
# p <- p + theme(plot.margin = unit(c(1,3.5,1,1), "lines")) 
# p <- p + coord_cartesian(clip = "off") 

pdf(paste0("figure8a_", ".pdf"), width = 2, height = 3.5)
p
dev.off()




# get_mean_sd_norm = function(data) {
#   in_field = "DATA"
#   traceScr = data[, grepl(in_field, names(data))][,1:nt]
#   traceScr[] <- lapply(traceScr, function(x) {
#     as.numeric(as.character(x)) 
#   })
#   # sapply(traceScr, class)
#   
#   traceScr = -traceScr
#   traceScr[traceScr < 0] = 0
#   
#   normTraceScr = sweep(traceScr, 1, apply(traceScr, 1, min), `-`)
#   normTraceKd = sweep(traceKd, 1, apply(traceKd, 1, min), `-`)
#   # apply(normTraceKd, 1, min)
# 
#   normTraceScrMax = apply(normTraceScr, 1, max)
#   normTraceScrMax[normTraceScrMax < 1e-9] = 1
#   normTraceScr = sweep(normTraceScr, 1, normTraceScrMax, `/`)
#   
#   meanTraceScr = apply(
#     normTraceScr,
#     2, mean, na.rm = T)
#   
#   seTraceScr = apply(
#     normTraceScr,
#     2, sd, na.rm = T) / sqrt(nrow(traceScr))
#   return(list(mean = meanTraceScr, sd = seTraceScr))
# }
# 
# # ccbase = cc[grepl("Base", cc$PROTOCOL) & !(duplicated(paste(cc$CELLID, cc$PROTOCOL))),]
# # plot((50:150) / 1000, -ccbase[8,][, grepl(in_field, names(ccbase))][,50:150])
# 
# nt = 1000
# 
# dataBase = get_mean_sd_norm(cc[grepl("Base", cc$PROTOCOL) & !(duplicated(paste(cc$CELLID, cc$PROTOCOL), fromLast = T)),])
# dataNBQX = get_mean_sd_norm(cc[grepl("NBQX", cc$PROTOCOL),])
# dataSulpiride = get_mean_sd_norm(cc[grepl("Sulpiride", cc$PROTOCOL),])
# dataQuinpirole = get_mean_sd_norm(cc[grepl("Quinpirole", cc$PROTOCOL),])
# dataSCH23390 = get_mean_sd_norm(cc[grepl("SCH23390", cc$PROTOCOL),])
# # dataMix = get_mean_sd(cc[grepl("Mix", cc$PROTOCOL),])
# 
# hdataBase = get_mean_sd_norm(hcc[grepl("Base", hcc$PROTOCOL) & !(duplicated(paste(hcc$CELLID, hcc$PROTOCOL))),])
# hdataNBQX = get_mean_sd_norm(hcc[grepl("NBQX", hcc$PROTOCOL),])
# 
# ylab = "Charge / max charge"
# 
# require(ggplot2)
# int = 1
# dfPlot = data.frame(x = rep(int * 1:nt, 7) / 1000,
#                     y = c(dataBase$mean, dataNBQX$mean, dataSulpiride$mean, dataQuinpirole$mean, dataSCH23390$mean,
#                           hdataBase$mean, hdataNBQX$mean),
#                     se = c(dataBase$sd, dataNBQX$sd, dataSulpiride$sd, dataQuinpirole$sd, dataSCH23390$sd,
#                            hdataBase$sd, hdataNBQX$sd),
#                     Type = c(rep("DA",nt), rep("DA NBQX", nt), 
#                              rep("DA Sulpiride", nt), rep("DA Quinpirole", nt), rep("DA SCH23390", nt),
#                              rep("HC", nt), rep("HC NBQX", nt)),
#                     Group = factor(c(rep("DA",nt*5),rep("HC",nt*2))))
# p <- ggplot(dfPlot, aes(x = x, y = y, color = Type))
# p <- p + geom_line()
# p <- p + geom_ribbon(aes(ymin = y - se, ymax = y + se, fill = Type),
#                      color = "transparent", alpha = 0.15, show.legend = F)
# p <- p + labs(fill = "", x = "Time [s]", y = ylab)
# p <- p + scale_color_manual(values = c("turquoise", "blue", "darkslateblue", "lightblue", "slateblue1", "gray30", "gray50"))
# p <- p + scale_fill_manual(values =c("turquoise", "blue", "darkslateblue", "lightblue", "slateblue1", "gray30", "gray50"))
# p <- p + theme_classic()
# p <- p + facet_grid(. ~ Group)
# 
# pdf(paste0("figure8b", ".pdf"), width = 6, height = 3.5)
# p
# dev.off()


# cc2 = read.csv("20cumtrapz.csv", sep = " ")
# 
# cellid = gsub("(.*)/cell_", "\\1/", as.character(cc2$FILENAME))
# cc2$DATE = gsub("(\\d+)[AB]?/(\\d+)_\\d{4}.abf$", "\\1", cellid)
# cc2$CELLNO = as.numeric(gsub("(.+)/(\\d+)_\\d{4}.abf$", "\\2", cellid))
# cc2$RECNO = as.numeric(gsub("(.+)/(\\d+)_(\\d{4}).abf$", "\\3", cellid))
# cc2$CELLID = paste0(cc2$DATE, "_", cc2$CELLNO)
# 
# cc2 = cc2[cc2$CELLID %in% res$CELLID[glut],]
# cc2 = cc2[cc2$DATA_1000 > -100,]
# 
# nt = 1000
# dataNBQX = get_mean_sd(cc2[grepl("NBQX", cc2$PROTOCOL),])
# dataMix = get_mean_sd(cc2[grepl("Mix", cc2$PROTOCOL),])
# 
# ylab = "Charge [pC]"
# 
# require(ggplot2)
# int = 1
# dfPlot = data.frame(x = rep(int * 1:nt, 2) / 1000,
#                     y = c(dataNBQX$mean, dataMix$mean),
#                     se = c(dataNBQX$sd, dataMix$sd),
#                     Type = c(rep("DA NBQX", nt), rep("DA NBQX Sulpiride SCH23390", nt)))
# p <- ggplot(dfPlot, aes(x = x, y = y, color = Type))
# p <- p + geom_line(show.legend = F)
# p <- p + geom_ribbon(aes(ymin = y - se, ymax = y + se, fill = Type),
#                      color = "transparent", alpha = 0.15, show.legend = F)
# p <- p + labs(fill = "", x = "Time [s]", y = ylab)
# p <- p + scale_color_manual(values = c("turquoise", "blue", "darkslateblue", "lightblue", "slateblue1", "gray30", "gray50"))
# p <- p + scale_fill_manual(values =c("turquoise", "blue", "darkslateblue", "lightblue", "slateblue1", "gray30", "gray50"))
# p <- p + theme_classic()
# # p <- p + facet_grid(. ~ Group)
# 
# pdf(paste0("figure8c", ".pdf"), width = 2, height = 3.5)
# p
# dev.off()
# 
# 
# 
# 
# cc5 = read.csv("50cumtrapz.csv", sep = " ")
# 
# cellid = gsub("(.*)/cell_", "\\1/", as.character(cc5$FILENAME))
# cc5$DATE = gsub("(\\d+)[AB]?/(\\d+)_\\d{4}.abf$", "\\1", cellid)
# cc5$CELLNO = as.numeric(gsub("(.+)/(\\d+)_\\d{4}.abf$", "\\2", cellid))
# cc5$RECNO = as.numeric(gsub("(.+)/(\\d+)_(\\d{4}).abf$", "\\3", cellid))
# cc5$CELLID = paste0(cc5$DATE, "_", cc5$CELLNO)
# 
# cc5 = cc5[cc5$CELLID %in% res$CELLID[glut],]
# cc5 = cc5[cc5$DATA_1000 > -100,]
# 
# nt = 1000
# dataNBQX = get_mean_sd(cc5[grepl("NBQX", cc5$PROTOCOL),])
# dataMix = get_mean_sd(cc5[grepl("Mix", cc5$PROTOCOL),])
# 
# ylab = "Charge [pC]"
# 
# require(ggplot2)
# int = 1
# dfPlot = data.frame(x = rep(int * 1:nt, 2) / 1000,
#                     y = c(dataNBQX$mean, dataMix$mean),
#                     se = c(dataNBQX$sd, dataMix$sd),
#                     Type = c(rep("DA NBQX", nt), rep("DA NBQX Sulpiride SCH23390", nt)))
# p <- ggplot(dfPlot, aes(x = x, y = y, color = Type))
# p <- p + geom_line(show.legend = F)
# p <- p + geom_ribbon(aes(ymin = y - se, ymax = y + se, fill = Type),
#                      color = "transparent", alpha = 0.15, show.legend = F)
# p <- p + labs(fill = "", x = "Time [s]", y = ylab)
# p <- p + scale_color_manual(values = c("turquoise", "blue", "darkslateblue", "lightblue", "slateblue1", "gray30", "gray50"))
# p <- p + scale_fill_manual(values =c("turquoise", "blue", "darkslateblue", "lightblue", "slateblue1", "gray30", "gray50"))
# p <- p + theme_classic()
# # p <- p + facet_grid(. ~ Group)
# 
# pdf(paste0("figure8d", ".pdf"), width = 2, height = 3.5)
# p
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# cc2 = read.csv("20cumtrapz.csv", sep = " ")
# 
# cellid = gsub("(.*)/cell_", "\\1/", as.character(cc2$FILENAME))
# cc2$DATE = gsub("(\\d+)[AB]?/(\\d+)_\\d{4}.abf$", "\\1", cellid)
# cc2$CELLNO = as.numeric(gsub("(.+)/(\\d+)_\\d{4}.abf$", "\\2", cellid))
# cc2$RECNO = as.numeric(gsub("(.+)/(\\d+)_(\\d{4}).abf$", "\\3", cellid))
# cc2$CELLID = paste0(cc2$DATE, "_", cc2$CELLNO)
# 
# cc2 = cc2[cc2$CELLID %in% res$CELLID[glut],]
# cc2 = cc2[cc2$DATA_1000 > -100,]
# 
# nt = 1000
# dataNBQX = get_mean_sd(cc2[grepl("NBQX", cc2$PROTOCOL),])
# dataMix = get_mean_sd(cc2[grepl("Mix", cc2$PROTOCOL),])
# 
# ylab = "Charge [pC]"
# 
# require(ggplot2)
# int = 1
# dfPlot = data.frame(x = rep(int * 1:nt, 2) / 1000,
#                     y = c(dataNBQX$mean, dataMix$mean),
#                     se = c(dataNBQX$sd, dataMix$sd),
#                     Type = c(rep("DA NBQX", nt), rep("DA NBQX Sulpiride SCH23390", nt)))
# p <- ggplot(dfPlot, aes(x = x, y = y, color = Type))
# p <- p + geom_line(show.legend = F)
# p <- p + geom_ribbon(aes(ymin = y - se, ymax = y + se, fill = Type),
#                      color = "transparent", alpha = 0.15, show.legend = F)
# p <- p + labs(fill = "", x = "Time [s]", y = ylab)
# p <- p + scale_color_manual(values = c("turquoise", "blue", "darkslateblue", "lightblue", "slateblue1", "gray30", "gray50"))
# p <- p + scale_fill_manual(values =c("turquoise", "blue", "darkslateblue", "lightblue", "slateblue1", "gray30", "gray50"))
# p <- p + theme_classic()
# # p <- p + facet_grid(. ~ Group)
# 
# pdf(paste0("figure8c", ".pdf"), width = 2, height = 3.5)
# p
# dev.off()




cc1 = read.csv("1apcumtrapz.csv", sep = " ")

cellid = gsub("(.*)/cell_", "\\1/", as.character(cc1$FILENAME))
cc1$DATE = gsub("(\\d+)[AB]?/(\\d+)_\\d{4}.abf$", "\\1", cellid)
cc1$CELLNO = as.numeric(gsub("(.+)/(\\d+)_\\d{4}.abf$", "\\2", cellid))
cc1$RECNO = as.numeric(gsub("(.+)/(\\d+)_(\\d{4}).abf$", "\\3", cellid))
cc1$CELLID = paste0(cc1$DATE, "_", cc1$CELLNO)

cc1 = cc1[cc1$CELLID %in% res$CELLID[!silent],]

# cc1 = cc1[cc1$DATA_1000 > -100,]

# cc1[c("FILENAME", "DATA_1000")]

plot(1:991, -cc1[1, 10:1000], type = "l" ,ylim = c(-5, 20))
# for (j in 2:nrow(cc1)) {
#   lines(1:991, -cc1[j, 10:1000], type = "l")
# }

nt = 1000
data = get_mean_sd(cc1)

cc1$CULTURE = res$CULTURE[match(cc1$CELLID, res$CELLID)]
ttext = sprintf("%d/%d",nrow(cc1),length(unique(na.omit(cc1$CULTURE))))

ylab = "Cumulative charge [pC]"

require(ggplot2)
int = 1
dfPlot = data.frame(x = rep(int * 1:nt, 1) / 1000,
                    y = c(data$mean),
                    se = c(data$sd),
                    Type = c(rep("DA", nt)))

dftext = data.frame(Type = factor(c("DA")), label = c(ttext), stringsAsFactors = F)
dftext$x = mean(dfPlot$x)
dftext$y = max(dfPlot$y + dfPlot$se) + diff(range(dfPlot$y + dfPlot$se)) / 20

p <- ggplot(dfPlot, aes(x = x, y = y, color = Type))
p <- p + geom_line(show.legend = F)
p <- p + geom_ribbon(aes(ymin = y - se, ymax = y + se, fill = Type),
                     color = "transparent", alpha = 0.15, show.legend = F)
p <- p + geom_label(aes(label = label), 
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, colour = "black", 
                    show.legend = F) + expand_limits(y = dftext$value[1])
p <- p + labs(fill = "", x = "Time [s]", y = ylab)
p <- p + scale_color_manual(values = c("turquoise", "blue", "darkslateblue", "lightblue", "slateblue1", "gray30", "gray50"))
p <- p + scale_fill_manual(values =c("turquoise", "blue", "darkslateblue", "lightblue", "slateblue1", "gray30", "gray50"))
p <- p + theme_classic()
# p <- p + facet_grid(. ~ Group)

pdf(paste0("figure8e", ".pdf"), width = 2, height = 3.5)
p
dev.off()




nt = 250
data = get_mean_sd(cc1)

cc1$CULTURE = res$CULTURE[match(cc1$CELLID, res$CELLID)]
ttext = sprintf("%d/%d",nrow(cc1),length(unique(na.omit(cc1$CULTURE))))

# ylab = "Charge [pC]"

require(ggplot2)
int = 1
dfPlot = data.frame(x = rep(int * 1:nt, 1),
                    y = c(data$mean),
                    se = c(data$sd),
                    Type = c(rep("DA", nt)))

dftext = data.frame(Type = factor(c("DA")), label = c(ttext), stringsAsFactors = F)
dftext$x = mean(dfPlot$x)
dftext$y = max(dfPlot$y + dfPlot$se) + diff(range(dfPlot$y + dfPlot$se)) / 20

p <- ggplot(dfPlot, aes(x = x, y = y, color = Type))
p <- p + geom_line(show.legend = F)
p <- p + geom_ribbon(aes(ymin = y - se, ymax = y + se, fill = Type),
                     color = "transparent", alpha = 0.15, show.legend = F)
p <- p + geom_label(aes(label = label), 
                    data = dftext, stat = "identity", alpha = 0.5, colour = "black", 
                    show.legend = F) + expand_limits(y = dftext$value[1])
p <- p + labs(fill = "", x = "Time [ms]", y = ylab)
p <- p + scale_color_manual(values = c("gray50", "blue", "darkslateblue", "lightblue", "slateblue1", "gray30", "gray50"))
p <- p + scale_fill_manual(values =c("gray50", "blue", "darkslateblue", "lightblue", "slateblue1", "gray30", "gray50"))
p <- p + scale_x_continuous(breaks = c(50,150,250))
p <- p + theme_classic()
# p <- p + facet_grid(. ~ Group)

pdf(paste0("figure8e_", ".pdf"), width = 1.75, height = 3.5)
p
dev.off()

