
silent = res$amplitude1 < 10
glut = res$nbqx < 0.1 & !silent

# rres = res[glut,]
# rres$sulpiride_ = NA
# rres$sch23390_ = NA
# rres$quinpirole_ = NA
# for (i in 1:nrow(rres)) {
#   s = data[data$PROTOCOL == "0.2Hz_PP_Sulpiride.pro" & data$CELLID == rres$CELLID[i],]
#   j = which(data$NAME == s$NAME)
#   if (length(j) == 1) {  
#     rres$sulpiride_[i] = data[j-1,]$AMPLITUDE_REGION_00_SERIES_00_1.0.pA
#   }
# }
# for (i in 1:nrow(rres)) {
#   s = data[data$PROTOCOL == "0.2Hz_PP_Quinpirole.pro" & data$CELLID == rres$CELLID[i],]
#   j = which(data$NAME == s$NAME)
#   if (length(j) == 1) {  
#     rres$quinpirole_[i] = data[j-1,]$AMPLITUDE_REGION_00_SERIES_00_1.0.pA
#   }
# }
# for (i in 1:nrow(rres)) {
#   s = data[data$PROTOCOL == "0.2Hz_PP_SCH23390.pro" & data$CELLID == rres$CELLID[i],]
#   j = which(data$NAME == s$NAME)
#   if (length(j) == 1) {  
#     rres$sch23390_[i] = data[j-1,]$AMPLITUDE_REGION_00_SERIES_00_1.0.pA
#   }
# }
# 
# x = na.omit(rres[rres$sulpiride_ > 1000 & rres$sulpiride_ < 1500 & rres$sulpiride > 1.1,]$CELLID)
# s = data[data$PROTOCOL == "0.2Hz_PP_Sulpiride.pro" & data$CELLID == x,]
# s_ = as.character(data[which(data$NAME == s$NAME)-1,]$NAME)
# s = as.character(s$NAME)
# x = na.omit(rres[rres$quinpirole_ > 1000 & rres$quinpirole_ < 1500 & rres$quinpirole < 0.5,]$CELLID)
# x = sample(x, 1)
# q = data[data$PROTOCOL == "0.2Hz_PP_Quinpirole.pro" & data$CELLID == x,]
# q_ = as.character(data[which(data$NAME == q$NAME)-1,]$NAME)
# q = as.character(q$NAME)
# x = na.omit(rres[rres$sch23390_ > 1000 & rres$sch23390_ < 1500,]$CELLID)
# x = sample(x, 1)
# d = data[data$PROTOCOL == "0.2Hz_PP_SCH23390.pro" & data$CELLID == x,]
# d_ = as.character(data[which(data$NAME == d$NAME)-1,]$NAME)
# d = as.character(d$NAME)
# write.table(c(s, s_, q, q_, d, d_), "figure6a.txt", quote = F, row.names = F, col.names= F)

library(reshape2)


dat = melt(res[glut, c("CELLID", "CULTURE", "mix")])
dat = na.omit(dat)

dat$variable = as.character(dat$variable)
dat$variable[dat$variable == "mix"] = "NBQX\nSCH23390\nSulpiride"
dat$variable[dat$variable == "sulpiride"] = "Sulpiride"
dat$variable[dat$variable == "sch23390"] = "SCH23390"
dat$variable[dat$variable == "quinpirole"] = "Quinpirole"
dat$variable[dat$variable == "baclofen"] = "Baclofen"
dat$variable = factor(dat$variable, levels = c("SCH23390", "Sulpiride", "NBQX\nSCH23390\nSulpiride", "Quinpirole", "Baclofen"))

dftext = (dat %>% group_by(variable) %>% 
            summarise(n = n(), m = n_distinct(CULTURE)) %>% 
            mutate(label = sprintf("%d/%d", n, m)))[c("variable", "label")]
dftext$value = max(dat$value) + diff(range(dat$value)) / 20

p <- ggplot(dat, aes(x = variable, y = value, color = variable))
p <- p + geom_quasirandom(show.legend = F)
p <- p + stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult=1),
                      geom = "pointrange", color = "black",
                      position = position_dodge(0.8), show.legend = F)
p <- p + geom_label(aes(label = label), 
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, 
                    show.legend = F) + expand_limits(y = dftext$value[1])
p <- p + geom_hline(yintercept = 1, color = "black")
p <- p + scale_color_manual(values = c( "deeppink", "lightseagreen", "limegreen"))
p <- p + labs(fill = "", x = "", y = "Evoked PSC amplitude / control")
p <- p + theme_classic() #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste0("figure6b_i", ".pdf"), width = 1.75, height = 3.5)
p
dev.off()


dat = melt(res[glut, c("CELLID", "CULTURE", "sch23390")])
dat = na.omit(dat)

dat$variable = as.character(dat$variable)
dat$variable[dat$variable == "mix"] = "NBQX\nSCH23390\nSulpiride"
dat$variable[dat$variable == "sulpiride"] = "Sulpiride"
dat$variable[dat$variable == "sch23390"] = "SCH23390"
dat$variable[dat$variable == "quinpirole"] = "Quinpirole"
dat$variable[dat$variable == "baclofen"] = "Baclofen"
dat$variable = factor(dat$variable, levels = c("SCH23390", "Sulpiride", "NBQX\nSCH23390\nSulpiride", "Quinpirole", "Baclofen"))

dftext = (dat %>% group_by(variable) %>% 
            summarise(n = n(), m = n_distinct(CULTURE)) %>% 
            mutate(label = sprintf("%d/%d", n, m)))[c("variable", "label")]
#dftext$value = max(dat$value) + diff(range(dat$value)) / 20

dftext$value = 1.2

p <- ggplot(dat, aes(x = variable, y = value, color = variable))
p <- p + geom_quasirandom(show.legend = F)
p <- p + stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult=1),
                      geom = "pointrange", color = "black",
                      position = position_dodge(0.8), show.legend = F)
p <- p + geom_label(aes(label = label), 
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, 
                    show.legend = F) + expand_limits(y = dftext$value[1])
p <- p + geom_hline(yintercept = 1, color = "black")
p <- p + scale_color_manual(values = c("darkslateblue", "lightblue", "deeppink", "lightseagreen", "limegreen"))
p <- p + labs(fill = "", x = "", y = "Evoked PSC amplitude / control")
p <- p + theme_classic() #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste0("figure6b_ii", ".pdf"), width = 1.75, height = 3.5)
p
dev.off()


dat = melt(res[glut, c("CELLID", "CULTURE", "quinpirole", "sulpiride", "baclofen")])
dat = na.omit(dat)

dat$variable = as.character(dat$variable)
dat$variable[dat$variable == "mix"] = "NBQX\nSCH23390\nSulpiride"
dat$variable[dat$variable == "sulpiride"] = "Sulpiride"
dat$variable[dat$variable == "sch23390"] = "SCH23390"
dat$variable[dat$variable == "quinpirole"] = "Quinpirole"
dat$variable[dat$variable == "baclofen"] = "Baclofen"
dat$variable = factor(dat$variable, levels = c("SCH23390", "Sulpiride", "NBQX\nSCH23390\nSulpiride", "Quinpirole", "Baclofen"))

dftext = (dat %>% group_by(variable) %>% 
            summarise(n = n(), m = n_distinct(CULTURE)) %>% 
            mutate(label = sprintf("%d/%d", n, m)))[c("variable", "label")]
dftext$value = max(dat$value) + diff(range(dat$value)) / 20

p <- ggplot(dat, aes(x = variable, y = value, color = variable))
p <- p + geom_quasirandom(show.legend = F)
p <- p + stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult=1),
                      geom = "pointrange", color = "black",
                      position = position_dodge(0.8), show.legend = F)
p <- p + geom_label(aes(label = label), 
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, 
                    show.legend = F) + expand_limits(y = dftext$value[1])
p <- p + geom_hline(yintercept = 1, color = "black")
p <- p + scale_color_manual(values = c("lightblue", "lightseagreen", "limegreen"))
p <- p + labs(fill = "", x = "", y = "Evoked PSC amplitude / control")
p <- p + theme_classic() #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste0("figure6b_iii", ".pdf"), width = 2.5, height = 3.5)
p
dev.off()









# Mix

dat = melt(res[glut, c("CELLID", "CULTURE", "cmix")])
dat$value[dat$value < 0] = 0
dat = na.omit(dat)

dat$variable = as.character(dat$variable)
dat$variable[dat$variable == "cmix"] = "NBQX\nSCH23390\nSulpiride"
dat$variable[dat$variable == "csulpiride"] = "Sulpiride"
dat$variable[dat$variable == "csch23390"] = "SCH23390"
dat$variable[dat$variable == "cquinpirole"] = "Quinpirole"
dat$variable[dat$variable == "cbaclofen"] = "Baclofen"
dat$variable = factor(dat$variable, levels = c("SCH23390", "Sulpiride", "NBQX\nSCH23390\nSulpiride", "Quinpirole", "Baclofen"))

dftext = (dat %>% group_by(variable) %>% 
            summarise(n = n(), m = n_distinct(CULTURE)) %>% 
            mutate(label = sprintf("%d/%d", n, m)))[c("variable", "label")]
dftext$value = max(dat$value) + diff(range(dat$value)) / 20

p <- ggplot(dat, aes(x = variable, y = value, color = variable))
p <- p + geom_quasirandom(show.legend = F)
p <- p + stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult=1),
                      geom = "pointrange", color = "black",
                      position = position_dodge(0.8), show.legend = F)
p <- p + geom_label(aes(label = label), 
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, 
                    show.legend = F) + expand_limits(y = dftext$value[1])
p <- p + geom_hline(yintercept = 1, color = "black")
p <- p + scale_color_manual(values = c("deeppink", "lightseagreen", "limegreen"))
p <- p + labs(fill = "", x = "", y = "Evoked PSC charge / control")
p <- p + theme_classic() #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste0("figure6c_i", ".pdf"), width = 1.75, height = 3.5)
p
dev.off()

wilcox.test(dat$value, mu = 1.0)

# D1

dat = melt(res[glut, c("CELLID", "CULTURE",  "csch23390")])
dat$value[dat$value < 0] = 0
dat = na.omit(dat)

dat$variable = as.character(dat$variable)
dat$variable[dat$variable == "cmix"] = "NBQX\nSCH23390\nSulpiride"
dat$variable[dat$variable == "csulpiride"] = "Sulpiride"
dat$variable[dat$variable == "csch23390"] = "SCH23390"
dat$variable[dat$variable == "cquinpirole"] = "Quinpirole"
dat$variable[dat$variable == "cbaclofen"] = "Baclofen"
dat$variable = factor(dat$variable, levels = c("SCH23390", "Sulpiride", "NBQX\nSCH23390\nSulpiride", "Quinpirole", "Baclofen"))

dftext = (dat %>% group_by(variable) %>% 
            summarise(n = n(), m = n_distinct(CULTURE)) %>% 
            mutate(label = sprintf("%d/%d", n, m)))[c("variable", "label")]
dftext$value = max(dat$value) + diff(range(dat$value)) / 20

p <- ggplot(dat, aes(x = variable, y = value, color = variable))
p <- p + geom_quasirandom(show.legend = F)
p <- p + stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult=1),
                      geom = "pointrange", color = "black",
                      position = position_dodge(0.8), show.legend = F)
p <- p + geom_label(aes(label = label), 
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, 
                    show.legend = F) + expand_limits(y = dftext$value[1])
p <- p + geom_hline(yintercept = 1, color = "black")
p <- p + scale_color_manual(values = c("darkslateblue", "lightblue", "deeppink", "lightseagreen", "limegreen"))
p <- p + labs(fill = "", x = "", y = "Evoked PSC charge / control")
p <- p + theme_classic() #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste0("figure6c_ii", ".pdf"), width = 1.75, height = 3.5)
p
dev.off()


dat = melt(res[glut, c("CELLID", "CULTURE", "cquinpirole", "csulpiride", "cbaclofen")])
dat$value[dat$value < 0] = 0
dat = na.omit(dat)

dat$variable = as.character(dat$variable)
dat$variable[dat$variable == "cmix"] = "NBQX\nSCH23390\nSulpiride"
dat$variable[dat$variable == "csulpiride"] = "Sulpiride"
dat$variable[dat$variable == "csch23390"] = "SCH23390"
dat$variable[dat$variable == "cquinpirole"] = "Quinpirole"
dat$variable[dat$variable == "cbaclofen"] = "Baclofen"
dat$variable = factor(dat$variable, levels = c("SCH23390", "Sulpiride", "NBQX\nSCH23390\nSulpiride", "Quinpirole", "Baclofen"))

dftext = (dat %>% group_by(variable) %>% 
            summarise(n = n(), m = n_distinct(CULTURE)) %>% 
            mutate(label = sprintf("%d/%d", n, m)))[c("variable", "label")]
dftext$value = max(dat$value) + diff(range(dat$value)) / 20

p <- ggplot(dat, aes(x = variable, y = value, color = variable))
p <- p + geom_quasirandom(show.legend = F)
p <- p + stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult=1),
                      geom = "pointrange", color = "black",
                      position = position_dodge(0.8), show.legend = F)
p <- p + geom_label(aes(label = label), 
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, 
                    show.legend = F) + expand_limits(y = dftext$value[1])
p <- p + geom_hline(yintercept = 1, color = "black")
p <- p + scale_color_manual(values = c("lightblue", "lightseagreen", "limegreen"))
p <- p + labs(fill = "", x = "", y = "Evoked PSC charge / control")
p <- p + theme_classic() #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste0("figure6c_iii", ".pdf"), width = 2.5, height = 3.5)
p
dev.off()

wilcox.test(dat$value[dat$variable == "Sulpiride"], mu = 1.0)
wilcox.test(dat$value[dat$variable == "Quinpirole"], mu = 1.0)
wilcox.test(dat$value[dat$variable == "Baclofen"], mu = 1.0)

# 
# library(reshape2)
# dat = melt(res[glut, c("CELLID", "baclofen")])
# dat = na.omit(dat)
# 
# dat$variable = as.character(dat$variable)
# dat$variable[dat$variable == "sulpiride"] = "Sulpiride"
# dat$variable[dat$variable == "quinpirole"] = "Quinpirole"
# dat$variable[dat$variable == "baclofen"] = "Baclofen"
# dat$variable = factor(dat$variable)
# 
# length(unique(dat$CELLID))
# length(unique(dat$CELLID))
# 
# p <- ggplot(dat, aes(x=variable, y=value, fill = variable))
# p <- p + geom_dotplot(binaxis = "y", stackdir = "center",
#                       position = position_dodge(0.8))
# p <- p + labs(fill = "", x = "", y = "Amplitude / control amplitude")
# p <- p + stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult=1),
#                       geom = "pointrange", color = "black",
#                       position = position_dodge(0.8), show.legend = F)
# p <- p + geom_hline(yintercept = 1, color = "coral")
# p <- p + theme_classic()
# 
# pdf(paste0("figure6c", ".pdf"), width = 3.5, height = 3)
# p
# dev.off()

library(psych)

wilcox.test(res[glut, "nbqx"], res[glut, "mix"])

wilcox.test(res[glut, "sch23390"], mu = 1.0)
wilcox.test(res[glut, "csch23390"], mu = 1.0)

wilcox.test(res[glut, "sulpiride"], mu = 1.0)
wilcox.test(res[glut, "csulpiride"], mu = 1.0)

wilcox.test(res[glut, "baclofen"], mu = 1.0)
wilcox.test(res[glut, "cbaclofen"], mu = 1.0)

wilcox.test(res[glut, "mix"], mu = 1.0)
wilcox.test(res[glut, "cmix"], mu = 1.0)


wilcox.test(res[glut, "quinpirole"], mu = 1.0)
wilcox.test(res[glut, "cquinpirole"], mu = 1.0)

wilcox.test(res[!silent, "sulpiride"], mu = 1.0)

corr.test(res[glut, "quinpirole"], res[glut, "sulpiride"])
plot(res[glut, "quinpirole"], res[glut, "sulpiride"])

corr.test(res[glut, "quinpirole"], res[glut, "baclofen"])
plot(res[glut, "quinpirole"], res[glut, "baclofen"])
# 
# 
# dat = melt(res[c("CELLID", "ppf")])
# dat = na.omit(dat)
# 
# p <- ggplot(dat, aes(x = variable, y = value, fill = variable)) 
# p <- p + geom_dotplot(binaxis = "y", stackdir = "center", 
#                       position = position_dodge(0.8), show.legend = F) 
# p <- p + stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult=1), 
#                       geom = "pointrange", color = "black", 
#                       position = position_dodge(0.8), show.legend = F)
# p <- p + labs(fill = "", x = "", y = "Paired pulse ratio")
# p <- p + geom_hline(yintercept = 1, color = "coral")
# p <- p + theme_classic() + theme(axis.text.x=element_blank())
# 
# pdf(paste0("ppf", ".pdf"), width = 1.8, height = 4)
# p
# dev.off()
# 
# 
# dat = melt(res_[c("CELLID", "nbqx")])
# dat = na.omit(dat)
# 
# p <- ggplot(dat, aes(x = variable, y = value, fill = variable)) 
# p <- p + geom_dotplot(binaxis = "y", stackdir = "center", 
#                       position = position_dodge(0.8), show.legend = F, binwidth = 0.01) 
# p <- p + labs(fill = "", x = "", y = "Amplitude / control amplitude")
# p <- p + stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult=1), 
#                       geom = "pointrange", color = "black", 
#                       position = position_dodge(0.8), show.legend = F)
# p <- p + geom_hline(yintercept = 1, color = "coral")
# p <- p + theme_classic() + theme(axis.text.x=element_blank())
# 
# pdf(paste0("nbqx", ".pdf"), width = 1.8, height = 4)
# p
# dev.off()
# 
# suc <- read.xlsx("sucrose.xlsx", 1)
# 
# cellid = gsub("(.*)/cell_", "\\1/", as.character(suc$NAME))
# suc$DATE = gsub("(\\d+)[AB]?/(\\d+)_\\d{4}.abf$", "\\1", cellid)
# suc$CELLNO = as.numeric(gsub("(.+)/(\\d+)_\\d{4}.abf$", "\\2", cellid))
# suc$CELLID = paste0(suc$DATE, "_", suc$CELLNO)
# 
# suc = suc[suc$CELLID %in% res$CELLID,]
# suc["Sucrose"] = suc$CHARGE_SERIES_00_1.0.s.pA / 1000
# 
# dat = melt(suc[c("CELLID", "Sucrose")])
# dat = na.omit(dat)
# 
# p <- ggplot(dat, aes(x = variable, y = value, fill = variable)) 
# p <- p + geom_dotplot(binaxis = "y", stackdir = "center", 
#                       position = position_dodge(0.8), show.legend = F) 
# p <- p + labs(fill = "", x = "", y = "Charge [nC]")
# p <- p + stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult=1), 
#                       geom = "pointrange", color = "black", 
#                       position = position_dodge(0.8), show.legend = F)
# p <- p + theme_classic() + theme(axis.text.x=element_blank())
# 
# pdf(paste0("suc", ".pdf"), width = 1.8, height = 4)
# p
# dev.off()
# 
# ##
# 
# dat = melt(res[c("CELLID", "baclofen")])
# dat = na.omit(dat)
# 
# dat$variable = as.character(dat$variable)
# dat$variable[dat$variable == "baclofen"] = "Baclofen"
# dat$variable = factor(dat$variable)
# 
# dat = dat[dat$value < 5, ]
# length(unique(dat$CELLID))
# 
# p <- ggplot(dat, aes(x=variable, y=value, fill = variable)) 
# p <- p + geom_dotplot(binaxis = "y", stackdir = "center", 
#                       position = position_dodge(0.8), show.legend = F) 
# p <- p + labs(fill = "", x = "", y = "Amplitude / control amplitude")
# p <- p + stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult=1), 
#                       geom = "pointrange", color = "black", 
#                       position = position_dodge(0.8), show.legend = F)
# p <- p + geom_hline(yintercept = 1, color = "coral")
# p <- p + theme_classic() 
# 
# pdf(paste0("baclofen", ".pdf"), width = 3, height = 4)
# p
# dev.off()
# 
# 
# 

ppfres = res[glut, c("CELLID", "CULTURE", "quinpirole.relppf", "sulpiride.relppf", "baclofen.relppf", "sch23390.relppf", "nbqx.relppf")]

# ppfres[c("quinpirole.relppf", "sulpiride.relppf", "baclofen.relppf", "sch23390.relppf", "mix.relppf")] = 
# sweep(ppfres[c("quinpirole.relppf", "sulpiride.relppf", "baclofen.relppf", "sch23390.relppf", "mix.relppf")], 1,  res[glut, "ppf"], `/`)


# dat = melt(ppfres[c("CELLID", "CULTURE",  "mix.relppf")])
# dat = na.omit(dat)
# 
# dat$variable = as.character(dat$variable)
# dat$variable[dat$variable == "mix.relppf"] = "NBQX\nSCH23390\nSulpiride"
# dat$variable[dat$variable == "sulpiride.relppf"] = "Sulpiride"
# dat$variable[dat$variable == "sch23390.relppf"] = "SCH23390"
# dat$variable[dat$variable == "quinpirole.relppf"] = "Quinpirole"
# dat$variable[dat$variable == "baclofen.relppf"] = "Baclofen"
# dat$variable = factor(dat$variable, levels = c("SCH23390", "Sulpiride", "NBQX\nSCH23390\nSulpiride", "Quinpirole", "Baclofen"))
# 
# length(unique(dat$CELLID))
# length(unique(dat$CELLID))
# 
# table(dat$variable)
# 
# dftext = (dat %>% group_by(variable) %>% 
#             summarise(n = n(), m = n_distinct(CULTURE)) %>% 
#             mutate(label = sprintf("%d/%d", n, m)))[c("variable", "label")]
# dftext$value = max(dat$value) + diff(range(dat$value)) / 20
# 
# # 
# p <- ggplot(dat, aes(x = variable, y = value, color = variable))
# # p <- p + geom_bar(stat = "summary", fun.y = "mean", show.legend = F)
# # p <- p + stat_summary(geom = "bar", fun.y = mean, position = "dodge", show.legend = F)
# # p <- p + stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", show.legend = F, width = 0.2)
# p <- p + geom_quasirandom(show.legend = F)
# p <- p + stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult=1),
#                       geom = "pointrange", color = "black",
#                       position = position_dodge(0.8), show.legend = F)
# p <- p + geom_label(aes(label = label), 
#                     data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, 
#                     show.legend = F) + expand_limits(y = dftext$value[1])
# p <- p + geom_hline(yintercept = 1, color = "black")
# p <- p + scale_color_manual(values = c( "deeppink", "lightseagreen", "limegreen"))
# p <- p + labs(fill = "", x = "", y = "Paired pulse ratio / control")
# # p <- p + scale_y_continuous(labels=scales::percent) + scale_x_log10()
# p <- p + theme_classic() ##+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
# # p <- p + facet_grid(. ~ Group)
# 
# pdf(paste0("figure6d_i", ".pdf"), width = 1.75, height = 3.5)
# p
# dev.off()

dat = melt(ppfres[c("CELLID", "CULTURE",  "sch23390.relppf")])
dat = na.omit(dat)

dat$variable = as.character(dat$variable)
dat$variable[dat$variable == "mix.relppf"] = "NBQX\nSCH23390\nSulpiride"
dat$variable[dat$variable == "sulpiride.relppf"] = "Sulpiride"
dat$variable[dat$variable == "sch23390.relppf"] = "SCH23390"
dat$variable[dat$variable == "quinpirole.relppf"] = "Quinpirole"
dat$variable[dat$variable == "baclofen.relppf"] = "Baclofen"
dat$variable = factor(dat$variable, levels = c("SCH23390", "Sulpiride", "NBQX\nSCH23390\nSulpiride", "Quinpirole", "Baclofen"))

length(unique(dat$CELLID))
length(unique(dat$CELLID))

table(dat$variable)

dftext = (dat %>% group_by(variable) %>% 
            summarise(n = n(), m = n_distinct(CULTURE)) %>% 
            mutate(label = sprintf("%d/%d", n, m)))[c("variable", "label")]
dftext$value = max(dat$value) + diff(range(dat$value)) / 20

# 
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
p <- p + geom_hline(yintercept = 1, color = "black")
p <- p + scale_color_manual(values = c("darkslateblue", "lightblue", "deeppink", "lightseagreen", "limegreen"))
p <- p + labs(fill = "", x = "", y = "Paired pulse ratio / control")
# p <- p + scale_y_continuous(labels=scales::percent) + scale_x_log10()
p <- p + theme_classic() #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p <- p + facet_grid(. ~ Group)

pdf(paste0("figure6d_ii", ".pdf"), width = 1.75, height = 3.5)
p
dev.off()

dat = melt(ppfres[c("CELLID", "CULTURE", "quinpirole.relppf", "sulpiride.relppf", "baclofen.relppf")])
dat = na.omit(dat)

dat$variable = as.character(dat$variable)
dat$variable[dat$variable == "mix.relppf"] = "NBQX\nSCH23390\nSulpiride"
dat$variable[dat$variable == "sulpiride.relppf"] = "Sulpiride"
dat$variable[dat$variable == "sch23390.relppf"] = "SCH23390"
dat$variable[dat$variable == "quinpirole.relppf"] = "Quinpirole"
dat$variable[dat$variable == "baclofen.relppf"] = "Baclofen"
dat$variable = factor(dat$variable, levels = c("SCH23390", "Sulpiride", "NBQX\nSCH23390\nSulpiride", "Quinpirole", "Baclofen"))

length(unique(dat$CELLID))
length(unique(dat$CELLID))

table(dat$variable)

dftext = (dat %>% group_by(variable) %>% 
            summarise(n = n(), m = n_distinct(CULTURE)) %>% 
            mutate(label = sprintf("%d/%d", n, m)))[c("variable", "label")]
dftext$value = max(dat$value) + diff(range(dat$value)) / 20

# 
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
p <- p + geom_hline(yintercept = 1, color = "black")
p <- p + scale_color_manual(values = c("lightblue","lightseagreen", "limegreen"))
p <- p + labs(fill = "", x = "", y = "Paired pulse ratio / control")
# p <- p + scale_y_continuous(labels=scales::percent) + scale_x_log10()
p <- p + theme_classic() #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p <- p + facet_grid(. ~ Group)

pdf(paste0("figure6d_iii", ".pdf"), width = 2.5, height = 3.5)
p
dev.off()

plot(res$amplitude1[glut], ppfres$nbqx.relppf)

wilcox.test(ppfres$nbqx.relppf, mu = 1.0)

t.test(ppfres$sch23390.relppf, mu = 1.0)

wilcox.test(ppfres$sch23390.relppf, mu = 1.0)
wilcox.test(ppfres$sulpiride.relppf, mu = 1.0)
wilcox.test(ppfres$mix.relppf, mu = 1.0)
wilcox.test(ppfres$quinpirole.relppf, mu = 1.0)
wilcox.test(ppfres$baclofen.relppf, mu = 1.0)


wilcox.test(res[glut, "sch23390"], mu = 1.0)
wilcox.test(res[glut, "csch23390"], mu = 1.0)

wilcox.test(res[glut, "sulpiride"], mu = 1.0)
wilcox.test(res[glut, "csulpiride"], mu = 1.0)

wilcox.test(res[glut, "quinpirole"], mu = 1.0)
wilcox.test(res[glut, "cquinpirole"], mu = 1.0)

wilcox.test(res[glut, "baclofen"], mu = 1.0)
wilcox.test(res[glut, "cbaclofen"], mu = 1.0)
