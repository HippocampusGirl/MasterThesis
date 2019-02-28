


minisize = aggregate(base$CHARGE, list(id = base$CELLID), mean, na.omit = T)
suc$minisize = minisize$x[match(suc$CELLID, minisize$id)]

# hminisize = aggregate(hbase$CHARGE, list(id = hbase$CELLID), mean, na.omit = T)
# hsuc$minisize = hminisize$x[match(hsuc$CELLID, hminisize$id)]

suc$rrp = suc$Sucrose / suc$minisize
# hsuc$rrp = hsuc$Sucrose / hsuc$minisize

dat = melt(suc[sglut, c("CELLID", "CULTURE", "rrp")])
dat = na.omit(dat)
dat$Type = "DA"

# hdat = melt(hsuc[c("CELLID", "CULTURE", "rrp")])
# hdat = na.omit(hdat)
# hdat$Type = "HC"
# 
# dat = rbind(dat, hdat)

dat$Type = factor(dat$Type, levels = c("DA", "HC"))

dftext = (dat %>% group_by(Type) %>% 
            summarise(n = n(), m = n_distinct(CULTURE)) %>% 
            mutate(label = sprintf("%d/%d", n, m)))[c("Type", "label")]
dftext$value = max(dat$value) + diff(range(dat$value)) / 20

# p <- ggplot(dat, aes(x = value, color = Type)) 
# # p <- p + geom_histogram(bins = 20, aes(y =..count../sum(..count..))) 
# p <- p + stat_ecdf(geom = "step", pad = FALSE)
# p <- p + labs(fill = "", x = "# SV/RRP", y = "%")
# p <- p + scale_y_continuous(labels=scales::percent) 
# p <- p + scale_color_manual(values = c("turquoise", "gray50"))
# p <- p + theme_classic() 
# 
# pdf(paste0("figure3d", ".pdf"), width = 3, height = 3)
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
p <- p + labs(fill = "", x = "", y = "# SV/RRP")
# p <- p + scale_y_log10() + annotation_logticks(sides = "l")
# p <- p + scale_y_continuous(labels=scales::percent) + scale_x_log10()
p <- p + scale_color_manual(values = c("turquoise", "gray50"))
p <- p + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p <- p + facet_grid(. ~ Group)

pdf(paste0("figure3d_", ".pdf"), width = 1.75, height = 3.5)
p
dev.off()

describe(dat)

describe(suc[!sglut, c("CELLID", "CULTURE", "rrp")])
