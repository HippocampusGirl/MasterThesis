setwd("~/Dropbox/Univ/masterThesis/analysis/20180918_final/")

fres = res[res$risetime > 0 & res$risetime < 10000 & res$amplitude1 > 0  & !is.na(res$nbqx),]

silent = fres$amplitude1 < 10 
glut = fres$nbqx < 0.1 & !silent

sum(glut & !silent, na.rm = T) / sum(!silent, na.rm = T)

library(ggplot2)
library(ggbeeswarm)
library(dplyr)

# x = na.omit(fres[glut & fres$amplitude1 > 1000 & fres$amplitude1 < 1500,]$CELLID)
# x = sample(x, 1)
# s = data[data$PROTOCOL == "0.2Hz_PP_NBQX.pro" & data$CELLID == x,]
# s_ = as.character(data[which(data$NAME == s$NAME)-1,]$NAME)
# s = as.character(s$NAME)
# x = na.omit(hfres[hfres$nbqx < 0.1 & hfres$amplitude1 > 1000 & hfres$amplitude1 < 1500,]$CELLID)
# x = x[!(x %in% c("20180803_78"))]
# x = sample(x, 1)
# q = hippocampus[hippocampus$PROTOCOL == "0.2Hz_3PP_NBQX.pro" & hippocampus$CELLID == x,]
# q_ = as.character(hippocampus[which(hippocampus$NAME == q$NAME)-2,]$NAME)
# q = as.character(q$NAME)
# write.table(c(s, s_, q, q_), "figure1a.txt", quote = F, row.names = F, col.names= F)

library(freshape2)

dat = melt(fres[c("CELLID", "CULTURE", "amplitude1")])
dat = na.omit(dat)
dat$Type = "DA"

# hdat = melt(hfres[c("CELLID", "CULTURE","amplitude1")])
# hdat = na.omit(hdat)
# hdat$Type = "HC"

# dat = rbind(dat, hdat)

dat$variable = as.character(dat$variable)
dat$variable[dat$variable == "amplitude1"] = "Amplitude"
dat$variable = factor(dat$variable)

dat$Type = factor(dat$Type, levels = c("DA", "HC"))

library(dplyr)

# dat = dat %>% group_by(Type) %>% 
    # mutate(n = n(), m = n_distinct(CULTURE)) %>% 
      # mutate(label = sprintf("%d/%d", m, n))

dftext = (dat %>% group_by(Type) %>% 
  summarise(n = sum(is.finite(value)), m = n_distinct(CULTURE[is.finite(value)])) %>% 
  mutate(label = sprintf("%d/%d", n, m)))[c("Type", "label")]
dftext$value = max(dat$value) + diff(range(dat$value)) / 20

# p <- ggplot(dat, aes(x = value, color = Type)) 
# # p <- p + geom_histogram(bins = 20, aes(y =..count../sum(..count..)))
# p <- p + stat_ecdf(geom = "step", pad = FALSE)
# p <- p + labs(fill = "", x = "Uncorrected amplitude [pA]", y = "%")
# p <- p + scale_y_continuous(labels=scales::percent) + scale_x_log10(breaks = c(10, 100, 1000, 10000))
# p <- p + geom_vline(xintercept = 10, color = "lightblue")
# p <- p + scale_color_manual(values = c("gray50", "darkblue"))
# p <- p + theme_classic() 
# 
# pdf(paste0("figure1b", ".pdf"), width = 4, height = 3)
# p
# dev.off()



p <- ggplot(dat, aes(x = Type, y = value, color = Type)) 
# p <- p + geom_bar(stat = "summary", fun.y = "mean", show.legend = F)
# p <- p + stat_summary(geom = "bar", fun.y = mean, position = "dodge", show.legend = F) 
# p <- p + stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", show.legend = F, width = 0.2)
p <- p + geom_quasirandom(show.legend = F)
p <- p + stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult=1),
                      geom = "pointrange", color = "black",
                      position = position_dodge(0.8), show.legend = F)
p <- p + geom_label(aes(x = Type, label = label), 
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, colour = "black", 
                    show.legend = F) + expand_limits(y = dftext$value[1])
p <- p + labs(fill = "", x = "", y = "Uncorrected amplitude [pA]")
p <- p + geom_hline(yintercept = 10, color = "black")
# p <- p + scale_y_continuous(labels=scales::percent) + scale_x_log10()
p <- p + scale_color_manual(values = c("gray50", "gray50")) 
# p <- p + scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)) + annotation_logticks(sides = "l")
p <- p + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p <- p + facet_grid(. ~ Group)

pdf(paste0("figure1b_", ".pdf"), width = 2, height = 3.5)
p
dev.off()

describe(dat)

# p <- ggplot(fres, aes(x = amplitude1, y = nbqx.amplitude))
# p <- p + geom_point()
# p
# 
# p <- ggplot(fres, aes(x = nbqx.amplitude, y = mix.amplitude))
# p <- p + geom_point()
# p

# plot(na.omit(fres[, c("amplitude1", "nbqx.amplitude")]))
# plot(na.omit(fres[, c("nbqx.amplitude", "mix.amplitude")]))

dat = melt(fres[, c("CELLID", "CULTURE", "nbqx")])
dat = na.omit(dat)
dat$Type = "DA"

# hdat = melt(hfres[hfres$amplitude1 > 10, c("CELLID", "CULTURE", "nbqx")])
# hdat = na.omit(hdat)
# hdat$Type = "HC"
# 
# dat = rbind(dat, hdat)

dat$variable = as.character(dat$variable)
dat$variable[dat$variable == "amplitude1"] = "Amplitude"
dat$variable = factor(dat$variable)

dat$Type = factor(dat$Type, levels = c("DA", "HC"))

dftext = (dat %>% group_by(Type) %>% 
            summarise(n = sum(is.finite(value)), m = n_distinct(CULTURE[is.finite(value)])) %>% 
            mutate(label = sprintf("%d/%d", n, m)))[c("Type", "label")]
dftext$value = max(dat$value) + diff(range(dat$value)) / 20

# p <- ggplot(dat, aes(x = value, color = Type)) 
# # p <- p + geom_histogram(bins = 20, aes(y =..count../sum(..count..))) 
# p <- p + stat_ecdf(geom = "step", pad = FALSE)
# p <- p + labs(fill = "", x = "NBQX Amplitude / control amplitude", y = "%")
# p <- p + scale_y_continuous(labels=scales::percent) + scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1))
# p <- p + geom_vline(xintercept = 1, color = "coral")
# p <- p + geom_vline(xintercept = 0.1, color = "lightblue")
# p <- p + scale_color_manual(values = c("gray50", "darkblue"))
# p <- p + theme_classic() 
# 
# pdf(paste0("figure1c", ".pdf"), width = 4, height = 3)
# p
# dev.off()

options(scipen=10000)
p <- ggplot(dat, aes(x = Type, y = value, color = Type)) 
# p <- p + geom_bar(stat = "summary", fun.y = "mean", show.legend = F)
# p <- p + stat_summary(geom = "bar", fun.y = mean, position = "dodge", show.legend = F) 
# p <- p + stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", show.legend = F, width = 0.2)
# p <- p + geom_dotplot(binaxis = "y", stackdir = "center",
#                       position = position_dodge(0.8), show.legend = F, color = "transparent", binwidth = 0.1)
p <- p + geom_quasirandom(show.legend = F)
p <- p + stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult=1),
                      geom = "pointrange", color = "black",
                      position = position_dodge(0.8), show.legend = F)
p <- p + geom_label(aes(x = Type, label = label), 
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, colour = "black", 
                    show.legend = F) + expand_limits(y = dftext$value[1])
p <- p + geom_hline(yintercept = 0.1, color = "black")
# p <- p + scale_y_continuous(labels=scales::percent) + scale_x_log10()
p <- p + scale_color_manual(values = c("gray50", "gray50")) 
# p <- p + scale_y_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10)) + annotation_logticks(sides = "l")
p <- p + labs(fill = "", x = "", y = "NBQX evoked PSC amplitude / control")
p <- p + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p <- p + facet_grid(. ~ Group)

pdf(paste0("figure1c_", ".pdf"), width = 1.75, height = 3.5)
p
dev.off()


describe(100-dat$value*100)

table(fres$glut )

dftext = (fres %>% group_by(glut) %>%
            summarise(n = sum(is.finite(amplitude1)), m = n_distinct(CULTURE[is.finite(amplitude1)])) %>% 
            mutate(label = sprintf("%s\n%d/%d", glut, n, m)))[c("glut", "label")]
dftext = na.omit(dftext)

dftext$amplitude1 = c(100, 3000)
dftext$nbqx.amplitude = max(fres$nbqx.amplitude, na.rm = T) + diff(range(fres$nbqx.amplitude, na.rm = T)) *2

options(scipen=10000)
p <- ggplot(fres[!silent,], aes(x = amplitude1, y = nbqx.amplitude, color = glut)) 
p <- p + geom_point(show.legend = F)
p <- p + geom_label(aes(label = label),
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, 
                    show.legend = F) + expand_limits(y = dftext$nbqx.amplitude[1]+ diff(range(fres$nbqx.amplitude, na.rm = T)) *2)
p <- p + labs(fill = "", x = "Control amplitude [pC]", y = "NBQX evoked PSC Amplitude [pC]")
p <- p + geom_abline()
p <- p + scale_color_manual(values = c("turquoise", "orange")) 
# p <- p + scale_x_continuous(breaks = c(0, 5000, 10000))
p <- p + scale_y_log10() + scale_x_log10() + annotation_logticks(sides = "lb")
p <- p + theme_classic() # + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p <- p + facet_grid(. ~ Group)

pdf(paste0("figure1c__", ".pdf"), width = 4, height = 3.5)
p
dev.off()

fres$glut[!silent]


fres_ = fres[c("CELLID", "CULTURE")]
fres_$amplitude1 = fres[, c("amplitude1")] * (1 - fres[, c("nbqx")])

dat = melt(fres_)
dat = na.omit(dat)
dat$Type = "DA"

# fres_ = hfres[c("CELLID", "CULTURE")]
# fres_$amplitude1 = hfres[, c("amplitude1")] * (1 - hfres[, c("nbqx")])
# 
# hdat = melt(fres_)
# hdat = na.omit(hdat)
# hdat$Type = "HC"

# dat = rbind(dat, hdat)

dat$variable = as.character(dat$variable)
dat$variable[dat$variable == "amplitude1"] = "Amplitude"
dat$variable = factor(dat$variable)

dat$Type = factor(dat$Type, levels = c("DA", "HC"))

dat$value[dat$value < 0] = 0

dftext = (dat %>% group_by(Type) %>% 
            summarise(n = sum(is.finite(value)), m = n_distinct(CULTURE[is.finite(value)])) %>% 
            mutate(label = sprintf("%d/%d", n, m)))[c("Type", "label")]
dftext$value = max(dat$value) + diff(range(dat$value)) / 20

# p <- ggplot(dat, aes(x = value, color = Type)) 
# # p <- p + geom_histogram(bins = 20, aes(y =..count../sum(..count..))) 
# p <- p + stat_ecdf(geom = "step", pad = FALSE)
# p <- p + labs(fill = "", x = "Amplitude [pA]", y = "%")
# p <- p + scale_y_continuous(labels=scales::percent) + scale_x_log10(breaks = c(10, 100, 1000, 10000))
# p <- p + scale_fill_manual(values = c("gray50", "turquoise"))
# p <- p + theme_classic() 
# 
# pdf(paste0("figure1d", ".pdf"), width = 4, height = 3)
# p
# dev.off()


p <- ggplot(dat, aes(x = Type, y = value, color = Type)) 
# p <- p + geom_bar(stat = "summary", fun.y = "mean", show.legend = F)
# p <- p + stat_summary(geom = "bar", fun.y = mean, position = "dodge", show.legend = F) 
# p <- p + stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", show.legend = F, width = 0.2)
p <- p + geom_quasirandom(show.legend = F)
p <- p + stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult=1),
                      geom = "pointrange", color = "black",
                      position = position_dodge(0.8), show.legend = F)
p <- p + geom_label(aes(x = Type, label = label), 
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, colour = "black", 
                    show.legend = F) + expand_limits(y = dftext$value[1])
p <- p + labs(fill = "", x = "", y = "Amplitude [pA]")
p <- p + geom_hline(yintercept = 10, color = "black")
# p <- p + scale_y_continuous(labels=scales::percent) + scale_x_log10()
p <- p + scale_color_manual(values = c("turquoise", "gray50")) 
# p <- p + scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)) + annotation_logticks(sides = "l")
p <- p + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p <- p + facet_grid(. ~ Group)

pdf(paste0("figure1d_", ".pdf"), width = 1.75, height = 3.5)
p
dev.off()


library(psych)
describe(fres[,])
describe(fres[glut,])
describe(fres[!glut,])




