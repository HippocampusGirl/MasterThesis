setwd("~/Dropbox/Univ/masterThesis/analysis/20180918_final/")

library(ggplot2)

silent = fres$amplitude1 < 10 & fres$amplitude2 < 10
glut = fres$nbqx < 0.1 & !silent
# 
# x = na.omit(res[glut & res$amplitude1 > 800 & res$amplitude1 < 1200 & res$ppf > 1.3,]$CELLID)
# x = sample(x, 1)
# s = data[!duplicated(data$CELLID) & data$CELLID == x,]
# s = as.character(s$NAME)
# # x = na.omit(hres[hres$nbqx < 0.1 & hres$amplitude1 > 800 & hres$amplitude1 < 1200 & hres$ppf < 1.1,]$CELLID)
# # x = sample(x, 1)
# # q = hippocampus[!duplicated(hippocampus$CELLID) & hippocampus$CELLID == x,]
# # q = as.character(q$NAME)
# write.table(c(s, q), "figure4a.txt", quote = F, row.names = F, col.names= F)

dat = melt(fres[glut, c("CELLID", "CULTURE", "ppf")])
dat = na.omit(dat)
dat$Type = "DA"

# hdat = melt(hres[!(hres$amplitude1 < 10 & hres$amplitude2 < 10) & hres$nbqx < 0.1, c("CELLID", "CULTURE", "ppf")])
# hdat = na.omit(hdat)
# hdat$Type = "HC"
# 
# dat = rbind(dat, hdat)
# table(dat$Type)

dat$variable = as.character(dat$variable)
dat$variable[dat$variable == "ppf"] = "Paired pulse ratio"
dat$variable = factor(dat$variable)

dat$Type = factor(dat$Type, levels = c("HC", "DA"))

dftext = (dat %>% group_by(Type) %>% 
            summarise(n = n(), m = n_distinct(CULTURE)) %>% 
            mutate(label = sprintf("%d/%d", n, m)))[c("Type", "label")]
dftext$value = max(dat$value) + diff(range(dat$value)) / 20

# p <- ggplot(dat, aes(x = value, color = Type)) 
# # p <- p + geom_histogram(bins = 20, aes(y =..count../sum(..count..)))
# p <- p + stat_ecdf(geom = "step", pad = FALSE)
# p <- p + labs(fill = "", x = "Paired pulse ratio", y = "%")
# p <- p + scale_y_continuous(labels=scales::percent) + scale_x_log10(breaks = c(0.5, 1, 2, 4))
# p <- p + geom_vline(xintercept = 1, color = "coral")
# p <- p + scale_color_manual(values = c("gray50", "darkblue"))
# p <- p + theme_classic() 
# 
# pdf(paste0("figure4b", ".pdf"), width = 4, height = 3)
# p
# dev.off()

p <- ggplot(dat, aes(x = Type, y = value, color = Type)) 
# p <- p + geom_bar(stat = "summary", fun.y = "mean", show.legend = F)
# p <- p + stat_summary(geom = "bar", fun.y = mean, position = "dodge", show.legend = F) 
# p <- p + stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", show.legend = F, width = 0.2)
# p <- p + geom_dotplot(binaxis = "y", stackdir = "center",
#                       position = position_dodge(0.8), show.legend = F, color = "transparent", binwidth = 0.15)
p <- p + geom_quasirandom(show.legend = F)
p <- p + stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult=1),
                      geom = "pointrange", color = "black",
                      position = position_dodge(0.8), show.legend = F)
p <- p + geom_label(aes(x = Type, label = label), 
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, colour = "black", 
                    show.legend = F) + expand_limits(y = dftext$value[1])
p <- p + scale_y_continuous(breaks = (1:8)/2) 
p <- p + geom_hline(yintercept = 1, color = "black")
p <- p + labs(fill = "", x = "", y = "Paired pulse ratio")
# p <- p + scale_y_log10(breaks = c(0.5, 1, 2, 4)) + annotation_logticks(sides = "l")
p <- p + scale_color_manual(values = c("turquoise", "gray50")) 
p <- p + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p <- p + facet_grid(. ~ Group)

pdf(paste0("figure4b_", ".pdf"), width = 2, height = 3.5)
p
dev.off()

describe(dat)

wilcox.test(dat$value ~ dat$Type)

train = read.xlsx("train10.xlsx", 1)

cellid = gsub("(.*)/cell_", "\\1/", as.character(train$NAME))
train$DATE = gsub("(\\d+)B?/(\\d+)_\\d{4}.abf$", "\\1", cellid)
train$CELLNO = as.numeric(gsub("(.+)/(\\d+)_\\d{4}.abf$", "\\2", cellid))
train$CELLID = paste0(train$DATE, "_", train$CELLNO)

train = train[train$CELLID %in% res$CELLID,]

# htrain = read.xlsx("htrain10.xlsx", 1)
# 
# cellid = gsub("(.*)/cell_", "\\1/", as.character(htrain$NAME))
# htrain$DATE = gsub("(\\d+)B?/(\\d+)_\\d{4}.abf$", "\\1", cellid)
# htrain$CELLNO = as.numeric(gsub("(.+)/(\\d+)_\\d{4}.abf$", "\\2", cellid))
# htrain$CELLID = paste0(htrain$DATE, "_", htrain$CELLNO)
# 
# htrain = htrain[htrain$CELLID %in% hres$CELLID,]

train$amplitude1 = res$amplitude1[match(train$CELLID, res$CELLID)]
train$amplitude2 = res$amplitude2[match(train$CELLID, res$CELLID)]
train$nbqx = res$nbqx[match(train$CELLID, res$CELLID)]

# train = train[!(train$amplitude1 < 10 & train$amplitude2 < 10) & train$nbqx < 0.1, ]
# 
# htrain$amplitude1 = hres$amplitude1[match(htrain$CELLID, hres$CELLID)]
# htrain$amplitude2 = hres$amplitude2[match(htrain$CELLID, hres$CELLID)]
# htrain$nbqx = hres$nbqx[match(htrain$CELLID, hres$CELLID)]
# 
# htrain = htrain[!(htrain$amplitude1 < 10 & htrain$amplitude2 < 10) & htrain$nbqx < 0.1, ]


train$CULTURE = res$CULTURE[match(train$CELLID, res$CELLID)]
# htrain$CULTURE = hres$CULTURE[match(htrain$CELLID, hres$CELLID)]

ttext10 = sprintf("%d/%d", nrow(train), length(unique(na.omit(train$CULTURE))))
# httext = sprintf("%d/%d", length(unique(na.omit(htrain$CULTURE))), nrow(htrain))

in_field = "AMPLITUDE"
traceScr = train[, grepl(in_field, names(train))]
traceScr[] <- lapply(traceScr, function(x) {
  as.numeric(as.character(x)) 
})
sapply(traceScr, class)

# traceKd = htrain[, grepl(in_field, names(train))]
# traceKd[] <- lapply(traceKd, function(x) {
#   as.numeric(as.character(x)) 
# })
# sapply(traceKd, class)

normTraceScr = traceScr

# normTraceScr = sweep(traceScr, 1, apply(traceScr, 1, min), `-`)
# normTraceKd = sweep(traceKd, 1, apply(traceKd, 1, min), `-`)
# apply(normTraceKd, 1, min)

normTraceScrFirst = normTraceScr[,1]
normTraceScrFirst[normTraceScrFirst < 1e-9] = 1
normTraceScr = sweep(normTraceScr, 1, normTraceScrFirst, `/`)

# normTraceKdMax = apply(normTraceKd, 1, max)
# normTraceKdMax[normTraceKdMax < 1e-9] = 1
# normTraceKd = sweep(normTraceKd, 1, normTraceKdMax, `/`)

# apply(normTraceKd, 1, max)

if (in_field == "RISE_TIME") {
  normTraceScr = sweep(traceScr, 1, 1, `/`)
  normTraceKd = sweep(traceKd, 1, 1, `/`)
}

meanTraceScr = apply(
  normTraceScr,
  2, mean, na.rm = T)
# meanTraceKd = apply(
#   normTraceKd,
#   2, mean, na.rm = T)

meanTraceScr10 = meanTraceScr

seTraceScr = apply(
  normTraceScr,
  2, sd, na.rm = T) / sqrt(nrow(traceScr))

seTraceScr10 = seTraceScr
# seTraceKd = apply(
#   normTraceKd,
#   2, sd, na.rm = T) / sqrt(nrow(traceKd))

# print(t.test(apply(normTraceScr[,180:200], 1, mean), apply(normTraceKd[,180:200], 1, mean)))

# ylab = "Evoked PSC amplitude / max amplitude"

# if (in_field == "CHARGE") {
#   ylab = "Charge / Max charge"
# }
# if (in_field == "RISE_TIME") {
#   ylab = "Rise time "
# }
# 
# dftext = data.frame(Type = factor(c("DA", "HC")), label = c(ttext, httext), stringsAsFactors = F)
# dftext$x = 100
# dftext$y = c(meanTraceScr[100] - 0.1, meanTraceKd[100] + 0.1)
# 
# require(ggplot2)
# int = 1
# dfPlot = data.frame(x = rep(int * 1:length(meanTraceKd), 2),
#                     y = c(meanTraceScr, meanTraceKd),
#                     se = c(seTraceScr, seTraceKd),
#                     Type = c(rep("DA", length(meanTraceScr)), rep("HC", length(meanTraceScr))))
# p <- ggplot(dfPlot, aes(x = x, y = y, color = Type))
# p <- p + ggtitle("10Hz")
# p <- p + geom_line(show.legend = F)
# p <- p + geom_ribbon(aes(ymin = y - se, ymax = y + se, fill = Type),
#                      color = "transparent", alpha = 0.15, show.legend = F)
# p <- p + geom_text(aes(label = label), data = dftext, stat = "identity", position = position_dodge(0.9), show.legend = F) 
# p <- p + labs(fill = "", x = "Stimulus #", y = ylab)
# p <- p + scale_color_manual(values = c("turquoise", "gray50"))
# p <- p + scale_fill_manual(values = c("turquoise", "gray50"))
# p <- p + theme_classic()
# 
# pdf(paste0("figure4c", ".pdf"), width = 2, height = 3.5)
# p
# dev.off()






train = read.xlsx("train20.xlsx", 1)

cellid = gsub("(.*)/cell_", "\\1/", as.character(train$NAME))
train$DATE = gsub("(\\d+)[AB]?/(\\d+)_\\d{4}.abf$", "\\1", cellid)
train$CELLNO = as.numeric(gsub("(.+)/(\\d+)_\\d{4}.abf$", "\\2", cellid))
train$CELLID = paste0(train$DATE, "_", train$CELLNO)

train = train[train$CELLID %in% res$CELLID,]

# htrain = read.xlsx("htrain20.xlsx", 1)
# 
# cellid = gsub("(.*)/cell_", "\\1/", as.character(htrain$NAME))
# htrain$DATE = gsub("(\\d+)B?/(\\d+)_\\d{4}.abf$", "\\1", cellid)
# htrain$CELLNO = as.numeric(gsub("(.+)/(\\d+)_\\d{4}.abf$", "\\2", cellid))
# htrain$CELLID = paste0(htrain$DATE, "_", htrain$CELLNO)
# 
# htrain = htrain[htrain$CELLID %in% hres$CELLID,]

train$amplitude1 = res$amplitude1[match(train$CELLID, res$CELLID)]
train$amplitude2 = res$amplitude2[match(train$CELLID, res$CELLID)]
train$nbqx = res$nbqx[match(train$CELLID, res$CELLID)]

train = train[!(train$amplitude1 < 10 & train$amplitude2 < 10) & train$nbqx < 0.1, ]

# htrain$amplitude1 = hres$amplitude1[match(htrain$CELLID, hres$CELLID)]
# htrain$amplitude2 = hres$amplitude2[match(htrain$CELLID, hres$CELLID)]
# htrain$nbqx = hres$nbqx[match(htrain$CELLID, hres$CELLID)]
# 
# htrain = htrain[!(htrain$amplitude1 < 10 & htrain$amplitude2 < 10) & htrain$nbqx < 0.1, ]

train$CULTURE = res$CULTURE[match(train$CELLID, res$CELLID)]
# htrain$CULTURE = hres$CULTURE[match(htrain$CELLID, hres$CELLID)]

ttext20 = sprintf("%d/%d", nrow(train), length(unique(na.omit(train$CULTURE))))
# httext = sprintf("%d/%d", length(unique(na.omit(htrain$CULTURE))), nrow(htrain))

in_field = "AMPLITUDE"
traceScr = train[, grepl(in_field, names(train))]
traceScr[] <- lapply(traceScr, function(x) {
  as.numeric(as.character(x)) 
})
sapply(traceScr, class)

# traceKd = htrain[, grepl(in_field, names(train))]
# traceKd[] <- lapply(traceKd, function(x) {
#   as.numeric(as.character(x)) 
# })
# sapply(traceKd, class)

# traceScr[traceScr < 0] = 0
# traceKd[traceKd < 0] = 0

normTraceScr = traceScr

# normTraceScr = sweep(traceScr, 1, apply(traceScr, 1, min), `-`)
# normTraceKd = sweep(traceKd, 1, apply(traceKd, 1, min), `-`)
# apply(normTraceKd, 1, min)

normTraceScrFirst = normTraceScr[,1]
normTraceScr = sweep(normTraceScr, 1, normTraceScrFirst, `/`)
normTraceScr = normTraceScr[normTraceScrFirst > 1e-9,]

# normTraceKdMax = apply(normTraceKd, 1, max)
# normTraceKdMax[normTraceKdMax < 1e-9] = 1
# normTraceKd = sweep(normTraceKd, 1, normTraceKdMax, `/`)

# apply(normTraceScr, 1, max)

# if (in_field == "RISE_TIME") {
#   normTraceScr = sweep(traceScr, 1, 1, `/`)
#   normTraceKd = sweep(traceKd, 1, 1, `/`)
# }

meanTraceScr = apply(
  normTraceScr,
  2, mean, na.rm = T)
# meanTraceKd = apply(
#   normTraceKd,
#   2, mean, na.rm = T)

meanTraceScr20 = meanTraceScr

seTraceScr = apply(
  normTraceScr,
  2, sd, na.rm = T) / sqrt(nrow(traceScr))
# seTraceKd = apply(
#   normTraceKd,
#   2, sd, na.rm = T) / sqrt(nrow(traceKd))

seTraceScr20 = seTraceScr

# print(t.test(apply(normTraceScr[,180:200], 1, mean), apply(normTraceKd[,180:200], 1, mean)))

# ylab = "Amplitude / Max amplitude"

# if (in_field == "CHARGE") {
#   ylab = "Charge / Max charge"
# }
# if (in_field == "RISE_TIME") {
#   ylab = "Rise time "
# }

dftext = data.frame(Type = factor(c("10Hz", "20Hz")), 
                    label = paste0(c("10 Hz\n", "20 Hz\n"), c(ttext10, ttext20)), stringsAsFactors = F)

require(ggplot2)
int = 1
dfPlot = data.frame(x = rep(int * 1:length(meanTraceScr), 2),
                    y = c(meanTraceScr10, meanTraceScr20),
                    se = c(seTraceScr10, seTraceScr20),
                    Type = c(rep("10Hz", length(meanTraceScr)), rep("20Hz", length(meanTraceScr))))


dftext$x = mean(dfPlot$x) 
dftext$x = dftext$x + c(-25, 25)
dftext$y = max(dfPlot$y+dfPlot$se*0.9) 

p <- ggplot(dfPlot, aes(x = x, y = y, color = Type))
p <- p + geom_line(show.legend = F)
p <- p + geom_ribbon(aes(ymin = y - se, ymax = y + se, fill = Type),
                     color = "transparent", alpha = 0.15, show.legend = F)
p <- p + labs(fill = "", x = "Stimulus #", y = "Norm. evoked PSC amplitude")
p <- p + geom_label(data = dftext, aes(label = label), hjust = 0.5, show.legend = F)
p <- p + scale_color_manual(values = c("turquoise", "slateblue"))
p <- p + scale_fill_manual(values = c("turquoise", "slateblue"))
p <- p + theme_classic()
# p <- p + theme(plot.margin = unit(c(1,3.5,1,1), "lines")) 
# p <- p + coord_cartesian(clip = "off") 

pdf(paste0("figure4d", ".pdf"), width = 4.25, height = 3.5)
p
dev.off()


