setwd("~/Dropbox/Univ/masterThesis/analysis/20180918_final/")

silent = res$amplitude1 < 10 
glut = res$nbqx < 0.1 & !silent

hist(res[glut, ]$risetime)

# x = res[glut & res$amplitude1 > 300, ]
# set.seed(1234)
# x = x[sample(nrow(x)),]
# x = x[!duplicated(round(x$risetime/2)), ]
# x = as.character(
#   data[!duplicated(data$CELLID) & data$CELLID %in% x$CELLID,]$NAME)
# 
# # x = sample(x, 3)
# write.table(x, "figure2a.txt", quote = F, row.names = F, col.names= F)

# x = hres[hres$amplitude1 & hres$nbqx < 0.1, c("CELLID", "risetime")]
# x = x[x$risetime < 0.5,]
# x = as.character(
#   hippocampus[!duplicated(hippocampus$CELLID) & hippocampus$CELLID %in% x$CELLID,]$NAME)
# x = sample(x, 1)
# write.table(x, "figure2a_.txt", quote = F, row.names = F, col.names= F)

dat = melt(fres[, c("CELLID", "CULTURE", "risetime")])
dat = na.omit(dat)
dat$Type = "DA"

# hdat = melt(hres[hres$amplitude1 & hres$nbqx < 0.1, c("CELLID", "CULTURE", "risetime")])
# hdat = na.omit(hdat)
# hdat$Type = "HC"
# 
# dat = rbind(dat, hdat)

dat$variable = as.character(dat$variable)
dat$variable[dat$variable == "risetime"] = "Rise time"
dat$variable = factor(dat$variable)

dat$Type = factor(dat$Type, levels = c("DA", "HC"))

dat = dat[dat$value > 0 & dat$value < 10000,]

dftext = (dat %>% group_by(Type) %>% 
            summarise(n = n(), m = n_distinct(CULTURE)) %>% 
            mutate(label = sprintf("%d/%d", n, m)))[c("Type", "label")]
dftext$value = max(dat$value) + diff(range(dat$value)) / 20

# p <- ggplot(dat, aes(x = value, color = Type)) 
# # p <- p + geom_histogram(bins = 20, aes(y =..count../sum(..count..))) 
# p <- p + stat_ecdf(geom = "step", pad = FALSE)
# p <- p + labs(fill = "", x = "Rise time [ms]", y = "%")
# p <- p + scale_y_continuous(labels=scales::percent) + scale_x_log10()
# p <- p + scale_color_manual(values = c("gray50", "darkblue"))
# p <- p + theme_classic() 
# 
# pdf(paste0("figure2b", ".pdf"), width = 4, height = 3)
# p
# dev.off()


p <- ggplot(dat, aes(x = Type, y = value, color = Type)) 
# p <- p + geom_bar(stat = "summary", fun.y = "mean", show.legend = F)
# p <- p + stat_summary(geom = "bar", fun.y = mean, position = "dodge", show.legend = F) 
# p <- p + stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", show.legend = F, width = 0.2)
# p <- p + geom_dotplot(binaxis = "y", stackdir = "center",
#                       position = position_dodge(0.8), show.legend = F, color = "transparent", binwidth = 0.05)
p <- p + geom_quasirandom(show.legend = F)
p <- p + stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult=1),
                      geom = "pointrange", color = "black",
                      position = position_dodge(0.8), show.legend = F)
p <- p + geom_label(aes(x = Type, label = label), 
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, colour = "black", 
                    show.legend = F) + expand_limits(y = dftext$value[1])# p <- p + scale_y_continuous(labels=scales::percent) + scale_x_log10()
p <- p + scale_color_manual(values = c("gray50", "gray50")) 
p <- p + labs(fill = "", x = "", y = "Rise time [ms]")
# p <- p + scale_y_continuous(labels=scales::percent) + scale_x_log10()
p <- p + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p <- p + facet_grid(. ~ Group)

pdf(paste0("figure2b_", ".pdf"), width = 1.75, height = 3.5)
p
dev.off()

wilcox.test(dat$value ~ dat$Type)

describe(dat)

# 
# dat = melt(res[glut, c("CELLID", "npeaks")])
# dat = na.omit(dat)
# 
# dat$variable = as.character(dat$variable)
# dat$variable[dat$variable == "npeaks"] = "Number of peaks"
# dat$variable = factor(dat$variable)
# 
# p <- ggplot(dat, aes(x = value)) 
# p <- p + geom_histogram(bins = 4, aes(y =..count../sum(..count..)))
# # p <- p + stat_ecdf(geom = "step", pad = FALSE)
# p <- p + labs(fill = "", x = "Number of peaks", y = "%")
# p <- p + scale_y_continuous(labels=scales::percent) 
# p <- p + theme_classic() 
# 
# pdf(paste0("figure2c", ".pdf"), width = 3, height = 3)
# p
# dev.off()

describe(fres[,])
