setwd("~/Dropbox/Univ/masterThesis/analysis/20181112/")

library(ggplot2)
library(ggbeeswarm)
library(dplyr)

library(boot)
library(compare)

icc = read.csv("Results4.csv")
iccb = read.csv("Results4b.csv")

icc = rbind(icc, iccb)

icc$cell = gsub(pattern = ".+/(aut_.*)/.*.tif", "\\1", icc$File)

table(icc$cell)

icc1 = icc[grepl("unmix_11.tif", icc$File),]
table(icc1$cell)

hist(icc1$Std_C6)

staining1 = icc1[grepl("Homer1", icc1$cell),]
table(staining1$cell)

ishc = grepl("Hc", staining1$cell)
isstr = grepl("Str", staining1$cell)
isfn = grepl("FACSneg", staining1$cell)
isfm = grepl("FM", staining1$cell)
iskcl = !grepl("noKCl", staining1$cell)
sum(ishc)
sum(isstr)
sum(isfm)
sum(iskcl)
sum(!iskcl)
sum(!ishc&!isstr)

hist(staining1$Mean_C6)

group = factor(c("No", "Yes")[(!grepl("noKCl", icc1$cell))+1])
wilcox.test(icc1$Mean_C6 ~ group)

# nt = 30

breaks = c(-Inf, seq(0, 200, length.out = 20), Inf)

histc <- function(x) {
  h <- hist(x, breaks = breaks, include.lowest = T, plot = F)
  h$count[-length(breaks)+1] / length(x)
}

dns = aggregate(icc1$Mean_C6, by = list(id = icc1$cell), histc)
ishc = grepl("Hc", dns$id)
isstr = grepl("Str", dns$id)
iskcl = !grepl("noKCl", dns$id)

# apply(dns$x[iskcl,1:19], 2, sd, na.rm = T) / nrow(dns$x)
# 
plot(colMeans(dns$x[!iskcl,], na.rm = T), type = "l")
lines(colMeans(dns$x[iskcl,], na.rm = T), col = "red")

plot(colMeans(dns$x[iskcl,], na.rm = T) / colMeans(dns$x[!iskcl,], na.rm = T), type = "l")
# 
# 
# cutoff = quantile(staining1$Mean_C6[iskcl], probs = c(0.95))
# cutoff = breaks[which(colMeans(dns$x[!iskcl&ishc,], na.rm = T) < colMeans(dns$x[iskcl&ishc,], na.rm = T))[2]]

# hist(staining1$Mean_C4[!ishc&!isstr])
# hist(staining1$Mean_C6[!ishc&!isstr&iskcl&staining1$Mean_C6<1500], breaks = 20)
# staining1$File

get_mean_sd = function(traceScr) {
  meanTraceScr = apply(
    traceScr,
    2, mean, na.rm = T)
  
  seTraceScr = apply(
    traceScr,
    2, sd, na.rm = T) / sqrt(nrow(traceScr))
  return(list(mean = meanTraceScr, sd = seTraceScr))
}

nokcl = get_mean_sd(dns$x[!iskcl,])
yeskcl = get_mean_sd(dns$x[iskcl,])

ii = which(nokcl$mean < yeskcl$mean)[-1][1]+1
cutoff = breaks[ii]

rel = colMeans(dns$x[iskcl,], na.rm = T) / colMeans(dns$x[!iskcl,], na.rm = T)
mean(rel[ii:length(rel)])

nt = length(breaks[-1])-1

ttext = sprintf("%d/%d",nrow(dns),1)

ylab = "Count"
dfPlot = data.frame(x = rep(breaks[-length(breaks)][-1], 2),
                    y = c(nokcl$mean, yeskcl$mean),
                    se = c(nokcl$sd, yeskcl$sd),
                    Type = factor(c(rep("Control",nt), rep("+KCl",nt)), 
                                  levels = c("Control", "+KCl")))
# 
dftext = data.frame(Type = factor(c("Control", "+KCl")), label = c("Control", "+KCl"), stringsAsFactors = F)
xrange = diff(range(dfPlot$x[is.finite(dfPlot$x)]), na.rm = T)
xmean = mean(dfPlot$x[is.finite(dfPlot$x)])
dftext$x = c(xmean - xrange / 4, xmean + xrange / 4)
dftext$y = max(dfPlot$y + dfPlot$se) + diff(range(dfPlot$y + dfPlot$se)) / 10

p <- ggplot(dfPlot, aes(x = x, y = y, color = Type))
p <- p + geom_line(show.legend = F)
p <- p + geom_ribbon(aes(ymin = y - se, ymax = y + se, fill = Type),
                     color = "transparent", alpha = 0.15, show.legend = F)
# p <- p + geom_vline(xintercept = cutoff, color = "black")
p <- p + geom_label(aes(label = label),
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, 
                    show.legend = F) + expand_limits(y = dftext$value[1])
p <- p + labs(fill = "", x = "FM4-64FX fluorescence [AU]", y = "Proportion of synaptic puncta")
p <- p + scale_color_manual(values = c("turquoise", "blue", "darkslateblue", "lightblue", "purple", "gray30", "gray50"))
p <- p + scale_fill_manual(values =c("turquoise", "blue", "darkslateblue", "lightblue", "purple", "gray30", "gray50"))
# p <- p + scale_y_log10() + annotation_logticks(sides = "l")
p <- p + theme_classic()
# p <- p + theme(plot.margin = unit(c(1,3.5,1,1), "lines")) 
# p <- p + coord_cartesian(clip = "off") 

pdf(paste0("figure12a", ".pdf"), width = 2.5, height = 3.5)
p
dev.off()


get_mean_sd = function(traceScr) {
  meanTraceScr = apply(
    traceScr,
    2, mean, na.rm = T)
  
  seTraceScr = apply(
    traceScr,
    2, sd, na.rm = T) / sqrt(nrow(traceScr))
  return(list(mean = meanTraceScr, sd = seTraceScr))
}

nokcl = get_mean_sd(dns$x[!iskcl,])
yeskcl = get_mean_sd(dns$x[iskcl,])


ylab = "Count"
dfPlot = data.frame(x = rep(breaks[-length(breaks)][-1], 1),
                    y = c(yeskcl$mean / nokcl$mean),
                    se = c(rep(NA,nt)),
                    Type = factor(c(rep("+KCl",nt)), 
                                  levels = c("Control", "+KCl")))
# 
dftext = data.frame(Type = factor(c("Control", "+KCl")), label = c("Control", "+KCl"), stringsAsFactors = F)
xrange = diff(range(dfPlot$x[is.finite(dfPlot$x)]), na.rm = T)
xmean = mean(dfPlot$x[is.finite(dfPlot$x)])
dftext$x = c(xmean - xrange / 4, xmean + xrange / 4)
dftext$y = max(dfPlot$y + dfPlot$se) + diff(range(dfPlot$y + dfPlot$se)) / 10

p <- ggplot(dfPlot, aes(x = x, y = y, color = Type))
p <- p + geom_line(show.legend = F)
p <- p + geom_ribbon(aes(ymin = y - se, ymax = y + se, fill = Type),
                     color = "transparent", alpha = 0.15, show.legend = F)
p <- p + geom_vline(xintercept = cutoff, color = "black")
p <- p + geom_hline(yintercept = 1.0, color = "black")
# p <- p + geom_label(aes(label = label),
#                     data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5,
#                     show.legend = F) + expand_limits(y = dftext$value[1])
p <- p + labs(fill = "", x = "FM4-64FX fluorescence [AU]", y = "Proportion of synaptic puncta / control")
p <- p + scale_color_manual(values = c("blue", "lightblue", "purple", "gray30", "gray50"))
p <- p + scale_fill_manual(values =c("blue", "lightblue", "purple", "gray30", "gray50"))
# p <- p + scale_y_log10() + annotation_logticks(sides = "l")
p <- p + theme_classic()
# p <- p + theme(plot.margin = unit(c(1,3.5,1,1), "lines")) 
# p <- p + coord_cartesian(clip = "off") 

pdf(paste0("figure12a_", ".pdf"), width = 2.5, height = 3.5)
p
dev.off()


countcutoff = function(x) {
  sum(x > cutoff, na.rm = T) / sum(is.finite(x))
}

cnt = aggregate(staining1$Mean_C5, by = list(id = staining1$cell), countcutoff)
iskcl = !grepl("noKCl", cnt$id)

cnt2 = cnt[!iskcl,]
cnt = cnt[iskcl,]

ishc = grepl("Hc", cnt$id)
isstr = grepl("Str", cnt$id)

dat = data.frame(variable = factor(c("DA", "Striatum", "Hippocampus")[1+(isstr*1)+(ishc*2)]), value = 1 - cnt$x)
dat$CULTURE = 1

aggregate(dat$value, b=list(dat$variable), FUN = mean)
se = aggregate(dat$value, b=list(dat$variable), FUN = sd) 
se$x = se$x / sqrt(aggregate(dat$value, b=list(dat$variable), FUN = length)$x)

dftext = (dat %>% group_by(variable) %>% 
            summarise(n = sum(is.finite(value)), m = n_distinct(CULTURE[is.finite(value)])) %>% 
            mutate(label = sprintf("%d/%d", n, m)))[c("variable", "label")]
dftext$value = max(dat$value) + diff(range(dat$value)) / 10

p <- ggplot(dat, aes(x = variable, y = value, color = variable))
# p <- p + geom_bar(stat = "summary", fun.y = "mean", show.legend = F)
# p <- p + stat_summary(geom = "bar", fun.y = mean, position = "dodge", show.legend = F)
# p <- p + stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", show.legend = F, width = 0.2)
p <- p + geom_quasirandom(show.legend = F)
p <- p + stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult=1),
                      geom = "pointrange", color = "black",
                      position = position_dodge(0.8), show.legend = F)
p <- p + geom_label(aes(label = label), 
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, 
                    show.legend = F) + expand_limits(y = dftext$value[1])
p <- p + scale_color_manual(values = c("darkslateblue", "lightblue", "deeppink", "lightseagreen", "limegreen"))
p <- p + labs(fill = "", x = "", y = "Percentage of silent synapses")
p <- p + scale_y_continuous(labels=scales::percent, breaks = c(0, 0.25, 0.5, 0.75, 1))
p <- p + theme_classic() #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p <- p + facet_grid(. ~ Group)

pdf(paste0("figure12b", ".pdf"), width = 3.7, height = 3.5)
p
dev.off()

p <- ggplot(data.frame(x = log(staining1$Mean_C6), y = log(staining1$Mean_C3)),
       aes(x = x, y = y)) 
p <- p + geom_hex(bins = 128)



staining1$aboveCutoff = staining1$Mean_C6 > cutoff
staining1$active = factor(c("Silent", "Active")[staining1$aboveCutoff+1])

##### MUNC13

cnt = aggregate(staining1$Mean_C4, by = list(cell = staining1$cell, active = staining1$active), mean)
iskcl = !grepl("noKCl", cnt$cell)
ishc = grepl("Hc", cnt$cell)
isstr = grepl("Str", cnt$cell)
isfn = grepl("FACSneg", cnt$cell)

dat = data.frame(cell = cnt$cell, variable = factor(c("DA", "Striatum", "Hippocampus")[1+(isstr*1)+(ishc*2)]), 
                 condition = c("Control", "+KCl")[iskcl+1],
                 type = cnt$active, value = cnt$x)

r1 = aggregate(dat$value, by = list(type = dat$type, variable = dat$variable, condition = dat$condition), mean, na.rm = T, drop = F)
r1 = r1[r1$condition == "Control",]

# for (i in unique(dat$variable)) {
#   for (j in unique(dat$type)) {
#     dat[dat$type == j & dat$variable == i, "value"] =
#       dat[dat$type == j & dat$variable == i, "value"] / r1[r1$type == j & r1$variable == i, "x"]
# 
#   }
# }
dat$Type = dat$type

aggregate(dat$value, by = list(type = dat$type, variable = dat$variable, condition = dat$condition), mean, na.rm = T, drop = F)
aggregate(dat$value, by = list(type = dat$type, variable = dat$variable, condition = dat$condition), sd, na.rm = T, drop = F) 

dat = dat[dat$condition == "+KCl",]

ulen = function(x) {
  length(unique(x))
}
dftext = aggregate(dat$cell, by = list(variable = dat$variable), ulen, drop = F)
dftext$label = sprintf("%d/1", dftext$x)
dftext$value = max(dat$value) + diff(range(dat$value)) / 5
dftext$Type = "Active"

p <- ggplot(dat, aes(x = variable, y = value, colour = Type))
p <- p + geom_quasirandom(show.legend = T, dodge.width = 0.5)
p <- p + stat_summary(aes(fill = Type), fun.data = "mean_cl_normal", fun.args = list(mult=1),
                      geom = "pointrange", color = "black",
                      position = position_dodge(0.5), show.legend = F)
p <- p + geom_label(aes(label = label), colour = "black",
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5,
                    show.legend = F) + expand_limits(y = dftext$value[1])
p <- p + scale_color_manual(values = c("lightseagreen", "gray50"))
p <- p + labs(y = "Munc13 fluorescence [AU]", x = "")
# # p <- p + scale_y_log10()
p <- p + theme_classic() #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p <- p + facet_grid(. ~ Group)

pdf(paste0("figure12c", ".pdf"), width = 4, height = 3.5)
p
dev.off()

dadat_munc13 = dat[dat$variable == "DA",]
wilcox.test(dadat_munc13$value ~ dadat_munc13$type, paired = T)


dadat = dat[dat$variable == "Hippocampus",]
wilcox.test(dadat$value ~ dadat$type)

dadat = dat[dat$variable == "Striatum",]
wilcox.test(dadat$value ~ dadat$type)

#### Homer1

cnt = aggregate(staining1$Mean_C5, by = list(cell = staining1$cell, active = staining1$active), mean)
iskcl = !grepl("noKCl", cnt$cell)
ishc = grepl("Hc", cnt$cell)
isstr = grepl("Str", cnt$cell)
isfn = grepl("FACSneg", cnt$cell)

dat = data.frame(cell = cnt$cell, variable = factor(c("DA", "Striatum", "Hippocampus")[1+(isstr*1)+(ishc*2)]), 
                 condition = c("Control", "+KCl")[iskcl+1],
                 type = cnt$active, value = cnt$x)

r1 = aggregate(dat$value, by = list(type = dat$type, variable = dat$variable, condition = dat$condition), mean, na.rm = T, drop = F)
r1 = r1[r1$condition == "Control",]

# for (i in unique(dat$variable)) {
#   for (j in unique(dat$type)) {
#     dat[dat$type == j & dat$variable == i, "value"] = 
#       dat[dat$type == j & dat$variable == i, "value"] / r1[r1$type == j & r1$variable == i, "x"]
#     
#   }
# }
dat$Type = dat$type

aggregate(dat$value, by = list(type = dat$type, variable = dat$variable, condition = dat$condition), mean, na.rm = T, drop = F)
aggregate(dat$value, by = list(type = dat$type, variable = dat$variable, condition = dat$condition), sd, na.rm = T, drop = F) 

dat = dat[dat$condition == "+KCl",]

dftext = aggregate(dat$cell, by = list(variable = dat$variable), ulen, drop = F)
dftext$label = sprintf("%d/1", dftext$x)
dftext$value = max(dat$value) + diff(range(dat$value)) / 10
dftext$Type = "Active"

p <- ggplot(dat, aes(x = variable, y = value, colour = Type))
p <- p + geom_quasirandom(show.legend = T, dodge.width = 0.5)
p <- p + stat_summary(aes(fill = Type), fun.data = "mean_cl_normal", fun.args = list(mult=1),
                      geom = "pointrange", color = "black",
                      position = position_dodge(0.5), show.legend = F)
p <- p + geom_label(aes(label = label), colour = "black",
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5,
                    show.legend = F) + expand_limits(y = dftext$value[1])
p <- p + scale_color_manual(values = c("lightseagreen", "gray50"))
p <- p + labs(y = "Homer1 fluorescence [AU]", x = "")
# # p <- p + scale_y_log10()
p <- p + theme_classic() #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p <- p + facet_grid(. ~ Group)

pdf(paste0("figure12d", ".pdf"), width = 4, height = 3.5)
p
dev.off()

dadat_homer1 = dat[dat$variable == "DA",]
wilcox.test(dadat_homer1$value ~ dadat_homer1$type)

dadat = dat[dat$variable == "Hippocampus",]
wilcox.test(dadat$value ~ dadat$type)

dadat = dat[dat$variable == "Striatum",]
wilcox.test(dadat$value ~ dadat$type)

staining2 = icc1[grepl("RIM", icc1$cell),]
table(staining2$cell)

staining2$aboveCutoff = staining2$Mean_C6 > cutoff
staining2$active = factor(c("Silent", "Active")[staining2$aboveCutoff+1])

#### RIM

cnt = aggregate(staining2$Mean_C4, by = list(cell = staining2$cell, active = staining2$active), mean)
iskcl = !grepl("noKCl", cnt$cell)
ishc = grepl("Hc", cnt$cell)
isstr = grepl("Str", cnt$cell)
isfn = grepl("FACSneg", cnt$cell)

dat = data.frame(cell = cnt$cell, variable = factor(c("DA", "Striatum", "Hippocampus")[1+(isstr*1)+(ishc*2)]), 
                 condition = c("Control", "+KCl")[iskcl+1],
                 type = cnt$active, value = cnt$x)

r1 = aggregate(dat$value, by = list(type = dat$type, variable = dat$variable, condition = dat$condition), mean, na.rm = T, drop = F)
r1 = r1[r1$condition == "Control",]

# for (i in unique(dat$variable)) {
#   for (j in unique(dat$type)) {
#     dat[dat$type == j & dat$variable == i, "value"] = 
#       dat[dat$type == j & dat$variable == i, "value"] / r1[r1$type == j & r1$variable == i, "x"]
#     
#   }
# }
dat$Type = dat$type

aggregate(dat$value, by = list(type = dat$type, variable = dat$variable, condition = dat$condition), mean, na.rm = T, drop = F)
aggregate(dat$value, by = list(type = dat$type, variable = dat$variable, condition = dat$condition), sd, na.rm = T, drop = F) 

dat = dat[dat$condition == "+KCl",]

dftext = aggregate(dat$cell, by = list(variable = dat$variable), ulen, drop = T)
dftext$label = sprintf("%d/1", dftext$x)
dftext$value = max(dat$value) + diff(range(dat$value)) / 10
dftext$Type = "Active"

p <- ggplot(dat, aes(x = variable, y = value, colour = Type))
p <- p + geom_quasirandom(show.legend = T, dodge.width = 0.5)
p <- p + stat_summary(aes(fill = Type), fun.data = "mean_cl_normal", fun.args = list(mult=1),
                      geom = "pointrange", color = "black",
                      position = position_dodge(0.5), show.legend = F)
p <- p + geom_label(aes(label = label), colour = "black",
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5,
                    show.legend = F) + expand_limits(y = dftext$value[1])
p <- p + scale_color_manual(values = c("lightseagreen", "gray50"))
p <- p + labs(y = "RIM fluorescence [AU]", x = "")
# # p <- p + scale_y_log10()
p <- p + theme_classic() #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p <- p + facet_grid(. ~ Group)

pdf(paste0("figure12e", ".pdf"), width = 2.5, height = 3.5)
p
dev.off()

dadat_rim = dat[dat$variable == "DA",]
wilcox.test(dadat_rim$value ~ dadat_rim$type)


#### Bassoon

cnt = aggregate(staining2$Mean_C5, by = list(cell = staining2$cell, active = staining2$active), mean)
iskcl = !grepl("noKCl", cnt$cell)
ishc = grepl("Hc", cnt$cell)
isstr = grepl("Str", cnt$cell)
isfn = grepl("FACSneg", cnt$cell)

dat = data.frame(cell = cnt$cell, variable = factor(c("DA", "Striatum", "Hippocampus")[1+(isstr*1)+(ishc*2)]), 
                 condition = c("Control", "+KCl")[iskcl+1],
                 type = cnt$active, value = cnt$x)

r1 = aggregate(dat$value, by = list(type = dat$type, variable = dat$variable, condition = dat$condition), mean, na.rm = T, drop = F)
r1 = r1[r1$condition == "Control",]

# for (i in unique(dat$variable)) {
#   for (j in unique(dat$type)) {
#     dat[dat$type == j & dat$variable == i, "value"] = 
#       dat[dat$type == j & dat$variable == i, "value"] / r1[r1$type == j & r1$variable == i, "x"]
#     
#   }
# }
dat$Type = dat$type

aggregate(dat$value, by = list(type = dat$type, variable = dat$variable, condition = dat$condition), mean, na.rm = T, drop = F)
aggregate(dat$value, by = list(type = dat$type, variable = dat$variable, condition = dat$condition), sd, na.rm = T, drop = F) 

dat = dat[dat$condition == "+KCl",]

dftext = aggregate(dat$cell, by = list(variable = dat$variable), ulen, drop = T)
dftext$label = sprintf("%d/1", dftext$x)
dftext$value = max(dat$value) + diff(range(dat$value)) / 5
dftext$Type = "Active"

p <- ggplot(dat, aes(x = variable, y = value, colour = Type))
p <- p + geom_quasirandom(show.legend = T, dodge.width = 0.5)
p <- p + stat_summary(aes(fill = Type), fun.data = "mean_cl_normal", fun.args = list(mult=1),
                      geom = "pointrange", color = "black",
                      position = position_dodge(0.5), show.legend = F)
p <- p + geom_label(aes(label = label), colour = "black",
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5,
                    show.legend = F) + expand_limits(y = dftext$value[1])
p <- p + scale_color_manual(values = c("lightseagreen", "gray50"))
p <- p + labs(y = "Bassoon fluorescence [AU]", x = "")
# # p <- p + scale_y_log10()
p <- p + theme_classic() #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p <- p + facet_grid(. ~ Group)

pdf(paste0("figure12f", ".pdf"), width = 2.5, height = 3.5)
p
dev.off()

dadat_bsn = dat[dat$variable == "DA",]
wilcox.test(dadat_bsn$value ~ dadat_bsn$type)


