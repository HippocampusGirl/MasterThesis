setwd("~/Dropbox/Univ/masterThesis/analysis/20181112/")

library(ggplot2)
library(ggbeeswarm)
library(dplyr)

# icc = read.csv("Results2.csv")

# icc$cell = gsub(pattern = ".+/(aut_.*)/.*.tif", "\\1", icc$File)

table(icc$cell)

icc2 = icc[!grepl("FM", icc$File),]
table(icc2$cell)

icc2 = icc2[log(icc2$Area) < 0 & log(icc2$Area) > -4,]

hist(log(icc2$Area))

staining1 = icc2[grepl("Munc13_555_Homer1", icc2$File),]
hist(log(staining1$Area))

hist(log(staining1$Mean_C3[staining1$Mean_C3 > 1]), breaks = 200)
hist(staining1$Mean_C3[staining1$Mean_C3 > 1])

hist(log(staining1$Mean_C4[staining1$Mean_C4 > 1]), breaks = 200)
hist(staining1$Mean_C4[staining1$Mean_C4 > 1])

staining2 = icc2[grepl("Pclo", icc2$File),]
hist(log(staining2$Area))

hist(log(staining2$Mean_C3[staining2$Mean_C3 > 1]), breaks = 200)
hist(staining2$Mean_C3[staining2$Mean_C3 > 1])

hist(log(staining2$Mean_C4), breaks = 100)
hist(staining2$Mean_C4[staining2$Mean_C4 > 1])

staining3 = icc2[grepl("RIM", icc2$File),]
hist(log(staining2$Area))

hist(log(staining3$Mean_C3[staining3$Mean_C3 > 1]), breaks = 20)
hist(staining3$Mean_C3[staining3$Mean_C3 > 1])

hist(log(staining3$Mean_C4), breaks = 20)
hist(staining3$Mean_C4[staining3$Mean_C4 > 1])

icc2 = icc[!grepl("FM", icc$File),]
table(icc2$cell)

staining1 = icc2[grepl("Munc13_555_Homer1", icc2$File),]
staining1 = staining1[staining1$Mean_C3 > 1,]
# staining1 = icc2[grepl("aut_MAP2_350_Syb2_405_GFP_488_Munc13_555_Homer1_647", icc2$File),]
table(staining1$cell)
hist(log(staining1$Area))


plot(log(staining1$Mean_C3), log(staining1$Mean_C4))

hist(log(staining1$Mean_C3[staining1$Mean_C3 > 1]), breaks = 100)
hist(staining1$Mean_C3[staining1$Mean_C3 > 1], breaks = 400)

histc <- function(x) {
  h <- hist(x, breaks = breaks, include.lowest = T, plot = F)
  h$count[-length(breaks)+1] / length(x)
}

get_mean_sd = function(traceScr) {
  meanTraceScr = apply(
    traceScr,
    2, mean, na.rm = T)
  
  seTraceScr = apply(
    traceScr,
    2, sd, na.rm = T) / sqrt(nrow(traceScr))
  return(list(mean = meanTraceScr, sd = seTraceScr))
}

#### Munc13

breaks = c(-Inf, seq(0, 2500, length.out = 20), Inf)
dns = aggregate(staining1$Mean_C3, by = list(id = staining1$cell), histc)

plot(colMeans(dns$x, na.rm = T), type = "l")

nokcl = get_mean_sd(dns$x)

dfPlot = data.frame(x = rep(breaks[-length(breaks)][-1], 1),
                    y = c(nokcl$mean),
                    se = c(nokcl$sd))

n = length(dns$id)
m = 1
dftext = data.frame(label = sprintf("%d/%d", n, m))
dftext$x = mean(dfPlot$x)
dftext$y = max(dfPlot$y) + max(dfPlot$se)

p <- ggplot(dfPlot, aes(x = x, y = y))
p <- p + geom_line(show.legend = F, color = "darkslateblue")
p <- p + geom_ribbon(aes(ymin = y - se, ymax = y + se),
                     color = "transparent", fill = "darkslateblue", alpha = 0.15, show.legend = F)
p <- p + geom_label(aes(label = label), color = "darkslateblue",
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, 
                    show.legend = F) 
p <- p + labs(fill = "", x = "Munc13 fluorescence [AU]", y = "Proportion of synaptic puncta")
p <- p + theme_classic()

pdf(paste0("figure13b", ".pdf"), width = 2.7, height = 3.5)
p
dev.off()


hist(log(staining1$Mean_C4[staining1$Mean_C4 > 1]), breaks = 100)
hist(staining1$Mean_C4[staining1$Mean_C4 > 1], breaks = 100)


##### HOMER1



breaks = c(-Inf, seq(0, 250, length.out = 20), Inf)
dns = aggregate(staining1$Mean_C4, by = list(id = staining1$cell), histc)

plot(colMeans(dns$x, na.rm = T), type = "l")

nokcl = get_mean_sd(dns$x)

dfPlot = data.frame(x = rep(breaks[-length(breaks)][-1], 1),
                    y = c(nokcl$mean),
                    se = c(nokcl$sd))

n = length(dns$id)
m = 1
dftext = data.frame(label = sprintf("%d/%d", n, m))
dftext$x = mean(dfPlot$x)
dftext$y = max(dfPlot$y) + max(dfPlot$se)

p <- ggplot(dfPlot, aes(x = x, y = y))
p <- p + geom_line(show.legend = F, color = "darkslateblue")
p <- p + geom_ribbon(aes(ymin = y - se, ymax = y + se),
                     color = "transparent", fill = "darkslateblue", alpha = 0.15, show.legend = F)
p <- p + geom_label(aes(label = label), color = "darkslateblue",
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, 
                    show.legend = F) 
p <- p + labs(fill = "", x = "Homer1 fluorescence [AU]", y = "Proportion of synaptic puncta")
p <- p + theme_classic()

pdf(paste0("figure13c", ".pdf"), width = 2.7, height = 3.5)
p
dev.off()




hist(log(staining1$Mean_C5[staining1$Mean_C5 > 1]), breaks = 200)
hist(staining1$Mean_C5[staining1$Mean_C5 > 1], breaks = 200)

hist(log(staining1$Mean_C6[staining1$Mean_C6 > 1]), breaks = 200)
hist(staining1$Mean_C6[staining1$Mean_C6 > 1], breaks = 200)

staining2 = icc2[grepl("Pclo_555_VGluT2_647", icc2$File),]

hist(log(staining2$Mean_C3[staining2$Mean_C3 > 1]), breaks = 100)
hist(staining2$Mean_C3[staining2$Mean_C3 > 1], breaks = 100)

plot(log(staining2$Mean_C3), log(staining2$Mean_C4))

#### Pclo


breaks = c(-Inf, seq(0, 4000, length.out = 20), Inf)
dns = aggregate(staining2$Mean_C3, by = list(id = staining2$cell), histc)

plot(colMeans(dns$x, na.rm = T), type = "l")

nokcl = get_mean_sd(dns$x)

dfPlot = data.frame(x = rep(breaks[-length(breaks)][-1], 1),
                    y = c(nokcl$mean),
                    se = c(nokcl$sd))

n = length(dns$id)
m = 1
dftext = data.frame(label = sprintf("%d/%d", n, m))
dftext$x = mean(dfPlot$x)
dftext$y = max(dfPlot$y) + max(dfPlot$se)

p <- ggplot(dfPlot, aes(x = x, y = y))
p <- p + geom_line(show.legend = F, color = "darkslateblue")
p <- p + geom_ribbon(aes(ymin = y - se, ymax = y + se),
                     color = "transparent", fill = "darkslateblue", alpha = 0.15, show.legend = F)
p <- p + geom_label(aes(label = label), color = "darkslateblue",
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, 
                    show.legend = F) 
p <- p + labs(fill = "", x = "Piccolo fluorescence [AU]", y = "Proportion of synaptic puncta")
p <- p + theme_classic()

pdf(paste0("figure13d", ".pdf"), width = 2.7, height = 3.5)
p
dev.off()

plot(staining2$Mean_C3, staining2$Mean_C4)

plot(staining2$Mean_C3, log(staining2$Mean_C4))

hist(log(staining2$Mean_C4[staining2$Mean_C4 > 1]), breaks = 100)
hist(staining2$Mean_C4[staining2$Mean_C4 > 1], breaks = 400)
hist(staining2$Mean_C3[staining2$Mean_C3 > 1], breaks = 400)


#### vglut2


breaks = c(-Inf, seq(0, 2000, length.out = 20), Inf)
dns = aggregate(staining2$Mean_C4, by = list(id = staining2$cell), histc)

plot(colMeans(dns$x, na.rm = T), type = "l")

nokcl = get_mean_sd(dns$x)

dfPlot = data.frame(x = rep(breaks[-length(breaks)][-1], 1),
                    y = c(nokcl$mean),
                    se = c(nokcl$sd))

n = length(dns$id)
m = 1
dftext = data.frame(label = sprintf("%d/%d", n, m))
dftext$x = mean(dfPlot$x)
dftext$y = max(dfPlot$y) + max(dfPlot$se)

p <- ggplot(dfPlot, aes(x = x, y = y))
p <- p + geom_line(show.legend = F, color = "darkslateblue")
p <- p + geom_ribbon(aes(ymin = y - se, ymax = y + se),
                     color = "transparent", fill = "darkslateblue", alpha = 0.15, show.legend = F)
p <- p + geom_label(aes(label = label), color = "darkslateblue",
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, 
                    show.legend = F) 
p <- p + labs(fill = "", x = "VGluT2 fluorescence [AU]", y = "Proportion of synaptic puncta")
p <- p + theme_classic()

pdf(paste0("figure13e", ".pdf"), width = 2.7, height = 3.5)
p
dev.off()


staining3 = icc2[grepl("RIM_555_Bsn_647", icc2$File),]
staining3 = staining3[staining3$Mean_C3 > 1,]

hist(log(staining3$Mean_C3[staining3$Mean_C3 > 1]), breaks = 100)
hist(staining3$Mean_C3[staining3$Mean_C3 > 1], breaks = 100)

plot(staining3$Mean_C3, staining3$Mean_C4)

#### RIM


breaks = c(-Inf, seq(0, 2000, length.out = 20), Inf)
dns = aggregate(staining3$Mean_C3, by = list(id = staining3$cell), histc)

plot(colMeans(dns$x, na.rm = T), type = "l")

nokcl = get_mean_sd(dns$x)

dfPlot = data.frame(x = rep(breaks[-length(breaks)][-1], 1),
                    y = c(nokcl$mean),
                    se = c(nokcl$sd))

n = length(dns$id)
m = 1
dftext = data.frame(label = sprintf("%d/%d", n, m))
dftext$x = mean(dfPlot$x)
dftext$y = max(dfPlot$y) + max(dfPlot$se)

p <- ggplot(dfPlot, aes(x = x, y = y))
p <- p + geom_line(show.legend = F, color = "darkslateblue")
p <- p + geom_ribbon(aes(ymin = y - se, ymax = y + se),
                     color = "transparent", fill = "darkslateblue", alpha = 0.15, show.legend = F)
p <- p + geom_label(aes(label = label), color = "darkslateblue",
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, 
                    show.legend = F) 
p <- p + labs(fill = "", x = "RIM fluorescence [AU]", y = "Proportion of synaptic puncta")
p <- p + theme_classic()

pdf(paste0("figure13f", ".pdf"), width = 2.7, height = 3.5)
p
dev.off()


hist(log(staining3$Mean_C4[staining3$Mean_C4 > 1]), breaks = 50)
hist(staining3$Mean_C4[staining3$Mean_C4 > 1], breaks = 50)


breaks = c(-Inf, seq(0, 300, length.out = 20), Inf)
dns = aggregate(staining3$Mean_C4, by = list(id = staining3$cell), histc)

plot(colMeans(dns$x, na.rm = T), type = "l")

nokcl = get_mean_sd(dns$x)

dfPlot = data.frame(x = rep(breaks[-length(breaks)][-1], 1),
                    y = c(nokcl$mean),
                    se = c(nokcl$sd))

n = length(dns$id)
m = 1
dftext = data.frame(label = sprintf("%d/%d", n, m))
dftext$x = mean(dfPlot$x)
dftext$y = max(dfPlot$y) + max(dfPlot$se)

p <- ggplot(dfPlot, aes(x = x, y = y))
p <- p + geom_line(show.legend = F, color = "darkslateblue")
p <- p + geom_ribbon(aes(ymin = y - se, ymax = y + se),
                     color = "transparent", fill = "darkslateblue", alpha = 0.15, show.legend = F)
p <- p + geom_label(aes(label = label), color = "darkslateblue",
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, 
                    show.legend = F) 
p <- p + labs(fill = "", x = "Bassoon fluorescence [AU]", y = "Proportion of synaptic puncta")
p <- p + theme_classic()

pdf(paste0("figure13g", ".pdf"), width = 2.7, height = 3.5)
p
dev.off()



