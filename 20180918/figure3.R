setwd("~/Dropbox/Univ/masterThesis/analysis/20180918_final/")

library(dplyr)
library(freshape2)
library(ggbeeswarm)
library(xlsx)

suc <- read.xlsx("sucrose.xlsx", 1)

# library(psych)
# plot(suc$CHARGE_SERIES_00_1.0.s.pA, -suc$Axograph)
# cor(suc$CHARGE_SERIES_00_1.0.s.pA, -suc$Axograph, use = "pairwise.complete.obs")
# corr.test(suc$CHARGE_SERIES_00_1.0.s.pA, -suc$Axograph)
# na.omit(cbind(suc$CHARGE_SERIES_00_1.0.s.pA, -suc$Axograph))

cellid = gsub("(.*)/cell_", "\\1/", as.character(suc$NAME))
suc$DATE = gsub("(\\d+)[AB]?/(\\d+)_\\d{4}.abf$", "\\1", cellid)
suc$CELLNO = as.numeric(gsub("(.+)/(\\d+)_\\d{4}.abf$", "\\2", cellid))
suc$CELLID = paste0(suc$DATE, "_", suc$CELLNO)

suc$FNAME = gsub("(\\d+)[AB]?/(\\d+_\\d{4}).abf$", "\\2", cellid)

suc = suc[suc$CELLID %in% fres$CELLID,]
suc["Sucrose"] = suc$CHARGE_SERIES_00_1.0.s.pA 

# include = c(
#   "0005_0007",
#   "0006_0007",
#   "0007_0007",
#   "0009_0007",
#   "0011_0007",
#   "0013_0007",
#   "0014_0007",
#   "0015_0007",
#   "0016_0007",
#   "0018_0007",
#   "0021_0007",
#   "0023_0007",
#   "0028_0007",
#   "0043_0009",
#   "0044_0009",
#   "0059_0009",
#   "0060_0009",
#   "0062_0009",
#   "0069_0009",
#   "0088_0009",
#   "0112_0010",
#   "0118_0009",
#   "0129_0007",
#   "0130_0007",
#   "0132_0007",
#   "0133_0009",
#   "0142_0007",
#   "0145_0007",
#   "0146_0007",
#   "0148_0007",
#   "0151_0007"
# )

# suc = suc[suc$FNAME %in% include,]

suc$CULTURE = fres$CULTURE[match(suc$CELLID, fres$CELLID)]
suc$amplitude1 = fres$amplitude1[match(suc$CELLID, fres$CELLID)]
suc$amplitude2 = fres$amplitude2[match(suc$CELLID, fres$CELLID)]
suc$nbqx = fres$nbqx[match(suc$CELLID, fres$CELLID)]

suc = suc[suc$amplitude1 > 10, ]
suc = suc[!is.na(suc$CELLID),]
suc = suc[suc$Sucrose > 0,]

sglut = suc$nbqx < 0.1

describe(suc[sglut,])
describe(suc[!sglut,])

# suc = suc[!sglut, ]



# hsuc <- read.xlsx("hsucrose.xlsx", 1)
# 
# cellid = gsub("(.*)/cell_", "\\1/", as.character(hsuc$NAME))
# hsuc$DATE = gsub("(\\d+)[AB]?/(\\d+)_\\d{4}.abf$", "\\1", cellid)
# hsuc$CELLNO = as.numeric(gsub("(.+)/(\\d+)_\\d{4}.abf$", "\\2", cellid))
# hsuc$CELLID = paste0(hsuc$DATE, "_", hsuc$CELLNO)
# 
# hsuc$FNAME = gsub("(\\d+)[AB]?/(\\d+_\\d{4}).abf$", "\\2", cellid)
# 
# hsuc = hsuc[hsuc$CELLID %in% hfres$CELLID,]
# hsuc["Sucrose"] = hsuc$CHARGE_SERIES_00_1.0.s.pA
# 
# # hinclude = c(
# #   "20180718/cell_05_0004",
# #   "20180718/cell_07_0004",
# #   "20180718/cell_08_0004",
# #   "20180719/cell_10_0004",
# #   "20180719/cell_11_0004",
# #   "20180719/cell_12_0004",
# #   "20180719/cell_13_0004",
# #   "20180719/cell_14_0004",
# #   "20180719/cell_16_0004",
# #   "0029_0004",
# #   "0030_0004",
# #   "0032_0004",
# #   "0033_0004",
# #   "0034_0004",
# #   "0035_0004",
# #   "0036_0004",
# #   "0037_0004",
# #   "0038_0004",
# #   "0051_0004",
# #   "0032_0004",
# #   "0053_0004",
# #   "0056_0004",
# #   "0057_0004",
# #   "0074_0004",
# #   "0075_0004",
# #   "0076_0004",
# #   "0077_0004",
# #   "0078_0004",
# #   "0079_0004",
# #   "0081_0004",
# #   "0083_0004",
# #   "0084_0004",
# #   "0085_0004",
# #   "0092_0004",
# #   "0093_0004",
# #   "0094_0004",
# #   "0095_0004",
# #   "0096_0004",
# #   "0098_0004",
# #   "0099_0004",
# #   "0100_0004",
# #   "0101_0004",
# #   "0103_0004",
# #   "0104_0004",
# #   "0105_0004",
# #   "0106_0004",
# #   "0108_0004",
# #   "0109_0004",
# #   "0110_0010"
# # )
# 
# # hsuc = hsuc[Reduce(`|`, lapply(hinclude, grepl, x = hsuc$NAME)),]
# 
# hsuc$CULTURE = hfres$CULTURE[match(hsuc$CELLID, hfres$CELLID)]
# hsuc$amplitude1 = hfres$amplitude1[match(hsuc$CELLID, hfres$CELLID)]
# hsuc$amplitude2 = hfres$amplitude2[match(hsuc$CELLID, hfres$CELLID)]
# hsuc$nbqx = hfres$nbqx[match(hsuc$CELLID, hfres$CELLID)]
# 
# hsuc = hsuc[hsuc$amplitude1 > 10 & hsuc$nbqx < 0.10, ]
# hsuc = hsuc[!is.na(hsuc$CELLID),]
# hsuc = hsuc[hsuc$Sucrose > 0,]




dat = melt(suc[sglut,c("CELLID", "CULTURE", "Sucrose")])
dat = na.omit(dat)
dat$Type = "DA"

# hdat = melt(hsuc[c("CELLID", "CULTURE", "Sucrose")])
# hdat = na.omit(hdat)
# hdat$Type = "HC"
# 
# dat = rbind(dat, hdat)

dat$Type = factor(dat$Type, levels = c("DA", "HC"))

dftext = (dat %>% group_by(Type) %>% 
            summarise(n = sum(is.finite(value)), m = n_distinct(CULTURE[is.finite(value)])) %>% 
            mutate(label = sprintf("%d/%d", n, m)))[c("Type", "label")]
dftext$value = max(dat$value) + diff(range(dat$value)) / 20 
# 
# p <- ggplot(dat, aes(x = value, color = Type)) 
# # p <- p + geom_histogram(bins = 20, aes(y =..count../sum(..count..))) 
# p <- p + stat_ecdf(geom = "step", pad = FALSE)
# p <- p + labs(fill = "", x = "Sucrose charge [pC]", y = "%")
# p <- p + scale_y_continuous(labels=scales::percent) 
# p <- p + scale_color_manual(values = c("turquoise", "gray50"))
# p <- p + theme_classic() 
# 
# pdf(paste0("figure3b", ".pdf"), width = 3, height = 3)
# p
# dev.off()

p <- ggplot(dat, aes(x = Type, y = value, color = Type)) 
# p <- p + geom_bar(stat = "summary", fun.y = "mean", show.legend = F)
# p <- p + stat_summary(geom = "bar", fun.y = mean, position = "dodge", show.legend = F) 
# p <- p + stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", show.legend = F, width = 0.2)
# p <- p + geom_dotplot(binaxis = "y", stackdir = "center",
#                       position = position_dodge(0.8), show.legend = F, color = "transparent")
p <- p + geom_quasirandom(show.legend = F)
p <- p + stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult=1),
                      geom = "pointrange", color = "black",
                      position = position_dodge(0.8), show.legend = F)
p <- p + geom_label(aes(x = Type, label = label), 
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, colour = "black", 
                    show.legend = F) + expand_limits(y = dftext$value[1])
p <- p + labs(fill = "", x = "", y = "Sucrose charge [pC]")
# p <- p + scale_y_log10() + annotation_logticks(sides = "l")
# p <- p + scale_y_continuous(labels=scales::percent) + scale_x_log10()
p <- p + scale_color_manual(values = c("turquoise", "gray50"))
p <- p + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p <- p + facet_grid(. ~ Group)

pdf(paste0("figure3b_", ".pdf"), width = 1.75, height = 3.5)
p
dev.off()

describe(dat)


ap <- read.xlsx("sucrose_1ap.xlsx", 1)

cellid = gsub("(.*)/cell_", "\\1/", as.character(ap$NAME))
ap$DATE = gsub("(\\d+)[AB]?/(\\d+)_\\d{4}.abf$", "\\1", cellid)
ap$CELLNO = as.numeric(gsub("(.+)/(\\d+)_\\d{4}.abf$", "\\2", cellid))
ap$CELLID = paste0(ap$DATE, "_", ap$CELLNO)

ap = ap[ap$CELLID %in% suc$CELLID,]
ap = ap[ap$CHARGE_REGION_00_SERIES_00_1.0.s.pA > 0,]

ap$CULTURE = suc$CULTURE[match(ap$CELLID, ap$CELLID)]
ap$Sucrose = suc$Sucrose[match(ap$CELLID, suc$CELLID)]
ap$pr = ap$CHARGE_REGION_00_SERIES_00_1.0.s.pA / ap$Sucrose 
ap$ppf = fres$ppf[match(ap$CELLID, fres$CELLID)]

aglut = sglut[match(ap$CELLID, suc$CELLID)]

library(psych)
corr.test(ap$ppf, ap$pr)
plot(ap$ppf, ap$pr)


# hap <- read.xlsx("hsucrose_1ap.xlsx", 1)
# 
# cellid = gsub("(.*)/cell_", "\\1/", as.character(hap$NAME))
# hap$DATE = gsub("(\\d+)[AB]?/(\\d+)_\\d{4}.abf$", "\\1", cellid)
# hap$CELLNO = as.numeric(gsub("(.+)/(\\d+)_\\d{4}.abf$", "\\2", cellid))
# hap$CELLID = paste0(hap$DATE, "_", hap$CELLNO)
# 
# hap = hap[hap$CELLID %in% hsuc$CELLID,]
# hap = hap[hap$CHARGE_REGION_00_SERIES_00_1.0.s.pA > 0,]
# 
# hap$CULTURE = hsuc$CULTURE[match(hap$CELLID, hsuc$CELLID)]
# hap$Sucrose = hsuc$Sucrose[match(hap$CELLID, hsuc$CELLID)]
# hap$pr = hap$CHARGE_REGION_00_SERIES_00_1.0.s.pA / hap$Sucrose 
# hap$ppf = hfres$ppf[match(hap$CELLID, hfres$CELLID)]

ap = ap[ap$pr < 1,]
# hap = hap[hap$pr < 1,]

# library(psych)
# corr.test(hap$ppf, hap$pr)
# plot(hap$ppf, hap$pr)



dat = melt(ap[aglut, c("CELLID", "CULTURE", "pr")])
dat = na.omit(dat)
dat$Type = "DA"

# hdat = melt(hap[c("CELLID", "CULTURE", "pr")])
# hdat = na.omit(hdat)
# hdat$Type = "HC"
# 
# dat = rbind(dat, hdat)

dat$Type = factor(dat$Type, levels = c("DA", "HC"))

dftext = (dat %>% group_by(Type) %>% 
            summarise(n = sum(is.finite(value)), m = n_distinct(CULTURE[is.finite(value)])) %>% 
            mutate(label = sprintf("%d/%d", n, m)))[c("Type", "label")]
dftext$value = max(dat$value) + diff(range(dat$value)) / 20

# p <- ggplot(dat, aes(x = value, color = Type)) 
# # p <- p + geom_histogram(bins = 20, aes(y =..count../sum(..count..))) 
# p <- p + stat_ecdf(geom = "step", pad = FALSE)
# p <- p + labs(fill = "", x = "Pvr", y = "%")
# p <- p + scale_y_continuous(labels=scales::percent) 
# p <- p + scale_color_manual(values = c("gray50", "turquoise"))
# p <- p + theme_classic() 
# 
# pdf(paste0("figure3c", ".pdf"), width = 3, height = 3)
# p
# dev.off()


p <- ggplot(dat, aes(x = Type, y = value, color = Type)) 
# p <- p + geom_bar(stat = "summary", fun.y = "mean", show.legend = F)
# p <- p + stat_summary(geom = "bar", fun.y = mean, position = "dodge", show.legend = F) 
# p <- p + stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", show.legend = F, width = 0.2)
# p <- p + geom_dotplot(binaxis = "y", stackdir = "center",
#                       position = position_dodge(0.8), show.legend = F, color = "transparent")
p <- p + geom_quasirandom(show.legend = F)
p <- p + stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult=1),
                      geom = "pointrange", color = "black",
                      position = position_dodge(0.8), show.legend = F)
p <- p + geom_label(aes(x = Type, label = label), 
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, colour = "black", 
                    show.legend = F) + expand_limits(y = dftext$value[1])
# p <- p + scale_y_log10() + annotation_logticks(sides = "l")
p <- p + labs(fill = "", x = "", y = "Pvr")
# p <- p + scale_y_continuous(labels=scales::percent) + scale_x_log10()
p <- p + scale_color_manual(values = c("turquoise", "gray50"))
p <- p + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p <- p + facet_grid(. ~ Group)

pdf(paste0("figure3c_", ".pdf"), width = 1.75, height = 3.5)
p
dev.off()


describe(dat)



describe(ap[!aglut, c("CELLID", "CULTURE", "pr")])
