setwd("~/Dropbox/Univ/masterThesis/analysis/20180918_final/")

library(ggplot2)
library(ggbeeswarm)
library(dplyr)

silent = res$amplitude1 < 10 
glut = res$nbqx < 0.1 & !silent

events = read.csv("events.csv")
events = events[events$AMPLITUDE > 5,]
events = events[events$CHARGE > 0,]

cellid = gsub("(.*)/cell_", "\\1/", as.character(events$FILENAME))
events$DATE = gsub("(\\d+)[AB]?/(\\d+)_\\d{4}.abf$", "\\1", cellid)
events$CELLNO = as.numeric(gsub("(.+)/(\\d+)_\\d{4}.abf$", "\\2", cellid))
events$RECNO = as.numeric(gsub("(.+)/(\\d+)_(\\d{4}).abf$", "\\3", cellid))
events$CELLID = paste0(events$DATE, "_", events$CELLNO)

events$PROTOCOL_RECNO = paste0(events$PROTOCOL, "_", events$RECNO)

events = events[events$CELLID %in% res$CELLID,]

base = events[events$PROTOCOL == "0.2Hz_PP_Base.pro" & events$RECNO < 2,]
nbqx = events[grepl("NBQX", events$PROTOCOL),]
sch = events[grepl("SCH23390", events$PROTOCOL),]
sulp = events[grepl("Sulpiride", events$PROTOCOL),]
quin = events[grepl("Quinpirole", events$PROTOCOL),]
mix = events[grepl("5AP_Mix", events$PROTOCOL),]

silent = res$amplitude1 < 10 & res$amplitude2 < 10
glut = res$nbqx < 0.1 & !silent





ulen <- function(x) {
  return(length(unique(x)))
}

freq = data.frame(CELLID = unique(events$CELLID))
freq$CULTURE = res$CULTURE[match(freq$CELLID, res$CELLID)]
freq$glut = res$glut[match(freq$CELLID, res$CELLID)]
freq = freq[!is.na(freq$glut),]

basefreq = aggregate(base$AMPLITUDE, list(id = base$CELLID), length) 
basefreq$x = basefreq$x / (aggregate(base$PROTOCOL_RECNO, list(id = base$CELLID), ulen)$x * 6 * 7)
mf = merge(freq, basefreq, by.x = "CELLID", by.y = "id", all.x = T)
freq$base = mf$x[match(as.character(freq$CELLID), mf$CELLID)]

nbqxfreq = aggregate(nbqx$AMPLITUDE, list(id = nbqx$CELLID), length) 
nbqxfreq$x = nbqxfreq$x / (aggregate(nbqx$PROTOCOL_RECNO, list(id = nbqx$CELLID), ulen)$x * 6 * 7)
mf = merge(freq, nbqxfreq, by.x = "CELLID", by.y = "id", all.x = T)
freq$nbqx = mf$x[match(as.character(freq$CELLID), mf$CELLID)]

sulpfreq = aggregate(sulp$AMPLITUDE, list(id = sulp$CELLID), length) 
sulpfreq$x = sulpfreq$x / (aggregate(sulp$PROTOCOL_RECNO, list(id = sulp$CELLID), ulen)$x * 6 * 7)
mf = merge(freq, sulpfreq, by.x = "CELLID", by.y = "id", all.x = T)
freq$sulp = mf$x[match(as.character(freq$CELLID), mf$CELLID)]

schfreq = aggregate(sch$AMPLITUDE, list(id = sch$CELLID), length) 
schfreq$x = schfreq$x / (aggregate(sch$PROTOCOL_RECNO, list(id = sch$CELLID), ulen)$x * 6 * 7)
mf = merge(freq, schfreq, by.x = "CELLID", by.y = "id", all.x = T)
freq$sch = mf$x[match(as.character(freq$CELLID), mf$CELLID)]

mixfreq = aggregate(mix$AMPLITUDE, list(id = mix$CELLID), length) 
mixfreq$x = mixfreq$x / (aggregate(mix$PROTOCOL_RECNO, list(id = mix$CELLID), ulen)$x * 6 * 7)
mf = merge(freq, mixfreq, by.x = "CELLID", by.y = "id", all.x = T)
freq$mix = mf$x[match(as.character(freq$CELLID), mf$CELLID)]

freq$corr = freq$base - freq$nbqx
freq$corr[freq$corr < 0] = 0

describe(freq)

filenames = sort(unique(events$FILENAME))



allfreq = aggregate(events$AMPLITUDE, list(id = events$CELLID, filename = events$FILENAME, protocol = events$PROTOCOL, recno = events$RECNO), length)
allfreq$x = allfreq$x / (aggregate(allfreq$recno, 
        list(id = allfreq$id, filename = allfreq$filename, protocol = allfreq$protocol, recno = allfreq$recno), ulen)$x * 6 * 7)

relfreq = data.frame(CELLID = unique(events$CELLID))
relfreq$CULTURE = res$CULTURE[match(relfreq$CELLID, res$CELLID)]
fields = c("nbqx", "sch", "sulp", "quin", "mix", "baclofen")
for (field in fields) {
  relfreq[field] = NA
  fieldfreq = allfreq[grepl(field, allfreq$protocol, ignore.case = T),]
  
  for (i in 1:nrow(fieldfreq)) {
    id = fieldfreq[i, "id"]
    recno = fieldfreq[i, "recno"]
    index0 = which(allfreq$id == id & allfreq$recno == (recno - 1) & grepl("Base", allfreq$protocol))
    index1 = which(allfreq$id == id & allfreq$recno == (recno + 1) & grepl("Base", allfreq$protocol))
    
    if (length(index0) == 1 && length(index1) == 1) {
      reference = (allfreq[index0, "x"] + allfreq[index1, "x"]) / 2
    } else if (length(index0) == 1) {
      reference = allfreq[index0, "x"] 
    } else if (length(index1) == 1) {
      reference = allfreq[index1, "x"]
    }
    
    relfreq[relfreq$CELLID == id,field] = fieldfreq[i, "x"] / reference
  }
}


wilcox.test(relfreq$mix, mu = 1.0)
wilcox.test(relfreq$sulp, mu = 1.0)
wilcox.test(relfreq$sch, mu = 1.0)
wilcox.test(relfreq$quin, mu = 1.0)

######

# 
# colMeans(freq[,3:7], na.rm =T)
# relfreq = cbind(freq[c("CELLID", "CULTURE")], sweep(freq[c("mix", "sulp", "sch")], 1, freq$base, `/`))
# colMeans(relfreq[, 3:5], na.rm = T)



#####

dat = melt(freq[c("CELLID", "CULTURE", "glut", "corr")])

# dat$variable = as.character(dat$variable)
# dat$variable[dat$variable == "corr"] = "DA"
# dat$variable = factor(dat$variable)

dat$variable = dat$glut

dftext = (dat %>% group_by(variable) %>% 
            summarise(n = sum(is.finite(value)), m = n_distinct(CULTURE[is.finite(value)])) %>% 
            mutate(label = sprintf("%d/%d", n, m)))[c("variable", "label")]
dftext$value = max(dat$value, na.rm = T) + diff(range(dat$value, na.rm = T)) / 20

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
p <- p + labs(fill = "", x = "", y = "Spontaneous PSC frequency [Hz]")
# p <- p + scale_y_continuous(labels=scales::percent) + scale_x_log10()
p <- p + scale_color_manual(values = c("turquoise", "orange")) 
# p <- p + scale_y_log10() + annotation_logticks(sides = "l")
p <- p + theme_classic() #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p <- p + facet_grid(. ~ Group)

pdf(paste0("figure5c__", ".pdf"), width = 2, height = 3.5)
p
dev.off()


describe(dat)

describe(freq[glut,"corr"])
describe(freq[!glut,"corr"])




wilcox.test(unlist(freq$nbqx), unlist(freq$mix))

#
amp = data.frame(CELLID = unique(events$CELLID))
amp$CULTURE = res$CULTURE[match(amp$CELLID, res$CELLID)]
amp$glut = res$glut[match(amp$CELLID, res$CELLID)]
amp = amp[!is.na(amp$glut),]

baseamp = aggregate(base$AMPLITUDE, list(id = base$CELLID), mean, na.omit = T)
baseamp$x = baseamp$x
mf = merge(amp, baseamp, by.x = "CELLID", by.y = "id", all.x = T)
amp$base = mf$x[match(as.character(amp$CELLID), mf$CELLID)]

nbqxamp = aggregate(nbqx$AMPLITUDE, list(id = nbqx$CELLID), mean, na.omit = T)
nbqxamp$x = nbqxamp$x
mf = merge(amp, nbqxamp, by.x = "CELLID", by.y = "id", all.x = T)
amp$nbqx = mf$x[match(as.character(amp$CELLID), mf$CELLID)]

sulpamp = aggregate(sulp$AMPLITUDE, list(id = sulp$CELLID), mean, na.omit = T)
sulpamp$x = sulpamp$x
mf = merge(amp, sulpamp, by.x = "CELLID", by.y = "id", all.x = T)
amp$sulp = mf$x[match(as.character(amp$CELLID), mf$CELLID)]

quinamp = aggregate(quin$AMPLITUDE, list(id = quin$CELLID), mean, na.omit = T)
quinamp$x = quinamp$x
mf = merge(amp, quinamp, by.x = "CELLID", by.y = "id", all.x = T)
amp$quin = mf$x[match(as.character(amp$CELLID), mf$CELLID)]

schamp = aggregate(sch$AMPLITUDE, list(id = sch$CELLID), mean, na.omit = T)
schamp$x = schamp$x 
mf = merge(amp, schamp, by.x = "CELLID", by.y = "id", all.x = T)
amp$sch = mf$x[match(as.character(amp$CELLID), mf$CELLID)]

mixamp = aggregate(mix$AMPLITUDE, list(id = mix$CELLID), mean, na.omit = T)
mixamp$x = mixamp$x
mf = merge(amp, mixamp, by.x = "CELLID", by.y = "id", all.x = T)
amp$mix = mf$x[match(as.character(amp$CELLID), mf$CELLID)]

amp$corr = amp$base - amp$nbqx
amp$corr[amp$corr < 0] = 0




allamp = aggregate(events$AMPLITUDE, 
                   list(id = events$CELLID, filename = events$FILENAME, protocol = events$PROTOCOL, recno = events$RECNO), mean, na.omit = T)
allamp$x = allamp$x 

relamp = data.frame(CELLID = unique(events$CELLID))
relamp$CULTURE = res$CULTURE[match(relamp$CELLID, res$CELLID)]
fields = c("nbqx", "sch", "sulp", "quin", "mix", "baclofen")
for (field in fields) {
  allamp[field] = NA
  fieldamp = allamp[grepl(field, allamp$protocol, ignore.case = T),]
  
  for (i in 1:nrow(fieldamp)) {
    id = fieldamp[i, "id"]
    recno = fieldamp[i, "recno"]
    index0 = which(allamp$id == id & allamp$recno == (recno - 1) & grepl("Base", allamp$protocol))
    index1 = which(allamp$id == id & allamp$recno == (recno + 1) & grepl("Base", allamp$protocol))
    
    if (length(index0) == 1 && length(index1) == 1) {
      reference = (allamp[index0, "x"] + allamp[index1, "x"]) / 2
    } else if (length(index0) == 1) {
      reference = allamp[index0, "x"] 
    } else if (length(index1) == 1) {
      reference = allamp[index1, "x"]
    }
    
    relamp[relamp$CELLID == id,field] = fieldamp[i, "x"] / reference
  }
}
describe(relamp)


# colMeans(amp[,3:7], na.rm =T)
# relamp = cbind(amp[c("CELLID", "CULTURE")], sweep(amp[c("mix", "sulp", "sch", "quin")], 1, amp$base, `/`))
# colMeans(relamp[, 3:6], na.rm = T)

dat = na.omit(melt(relamp[c("CELLID", "CULTURE", "mix", "sulp", "sch")]))

wilcox.test(relamp$mix, mu = 1.0)
wilcox.test(relamp$sulp, mu = 1.0)
wilcox.test(relamp$sch, mu = 1.0)
wilcox.test(relamp$quin, mu = 1.0)

describe(amp)

dat$variable = as.character(dat$variable)
# dat$variable[dat$variable == "base"] = "DA"
# dat$variable[dat$variable == "nbqx"] = "DA NBQX"
dat$variable[dat$variable == "mix"] = "DA NBQX\nSCH23390\nSulpiride"
dat$variable[dat$variable == "sulp"] = "DA Sulpiride"
dat$variable[dat$variable == "sch"] = "DA SCH23390"
dat$variable = factor(dat$variable, levels = c("DA SCH23390", "DA Sulpiride", "DA NBQX\nSCH23390\nSulpiride"))

dftext = (dat %>% group_by(variable) %>% 
            summarise(n = sum(is.finite(value)), m = n_distinct(CULTURE[is.finite(value)])) %>% 
            mutate(label = sprintf("%d/%d", n, m)))[c("variable", "label")]
dftext$value = max(dat$value, na.rm = T) + diff(range(dat$value, na.rm = T)) / 20

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
p <- p + scale_color_manual(values = c("darkslateblue", "lightblue", "deeppink"))
p <- p + geom_hline(yintercept = 1, color = "black")
p <- p + labs(fill = "", x = "", y = "Spontaneous PSC amplitude / control")
# p <- p + scale_y_continuous(labels=scales::percent) + scale_x_log10()
p <- p + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p <- p + facet_grid(. ~ Group)

pdf(paste0("figure5b_", ".pdf"), width = 2.5, height = 3.5)
p
dev.off()
# 
wilcox.test(amp$nbqx, amp$mix)
wilcox.test(amp$base, amp$sch)


wilcox.test(freq$corr)
mean(freq$corr, na.rm = T)

dat = melt(amp[c("CELLID", "CULTURE", "glut", "corr")])

# dat$variable = as.character(dat$variable)
# dat$variable[dat$variable == "corr"] = "DA"
# dat$variable = factor(dat$variable)

dat$variable = dat$glut

dftext = (dat %>% group_by(variable) %>% 
            summarise(n = sum(is.finite(value)), m = n_distinct(CULTURE[is.finite(value)])) %>% 
            mutate(label = sprintf("%d/%d", n, m)))[c("variable", "label")]
dftext$value = max(dat$value, na.rm = T) + diff(range(dat$value, na.rm = T)) / 20

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
p <- p + labs(fill = "", x = "", y = "Spontaneous PSC amplitude [pA]")
# p <- p + scale_y_continuous(labels=scales::percent) + scale_x_log10()
p <- p + scale_color_manual(values = c("turquoise", "orange")) 
# p <- p + scale_y_log10() + annotation_logticks(sides = "l")
p <- p + theme_classic() #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p <- p + facet_grid(. ~ Group)

pdf(paste0("figure5d__", ".pdf"), width = 2, height = 3.5)
p
dev.off()

# #
# # pdf(paste0("figure5c", ".pdf"), width = 4, height = 3)
# # p
# # dev.off()
# #
# #
# #




rise = data.frame(CELLID = unique(events$CELLID))
rise$CULTURE = res$CULTURE[match(rise$CELLID, res$CELLID)]
rise$glut = res$glut[match(rise$CELLID, res$CELLID)]
rise = rise[!is.na(rise$glut),]

baserise = aggregate(base$RISE_TIME, list(id = base$CELLID), mean, na.omit = T)
baserise$x = baserise$x
mf = merge(rise, baserise, by.x = "CELLID", by.y = "id", all.x = T)
rise$base = mf$x[match(as.character(rise$CELLID), mf$CELLID)]

nbqxrise = aggregate(nbqx$RISE_TIME, list(id = nbqx$CELLID), mean, na.omit = T)
nbqxrise$x = nbqxrise$x
mf = merge(rise, nbqxrise, by.x = "CELLID", by.y = "id", all.x = T)
rise$nbqx = mf$x[match(as.character(rise$CELLID), mf$CELLID)]

sulprise = aggregate(sulp$RISE_TIME, list(id = sulp$CELLID), mean, na.omit = T)
sulprise$x = sulprise$x
mf = merge(rise, sulprise, by.x = "CELLID", by.y = "id", all.x = T)
rise$sulp = mf$x[match(as.character(rise$CELLID), mf$CELLID)]

schrise = aggregate(sch$RISE_TIME, list(id = sch$CELLID), mean, na.omit = T)
schrise$x = schrise$x 
mf = merge(rise, schrise, by.x = "CELLID", by.y = "id", all.x = T)
rise$sch = mf$x[match(as.character(rise$CELLID), mf$CELLID)]

mixrise = aggregate(mix$RISE_TIME, list(id = mix$CELLID), mean, na.omit = T)
mixrise$x = mixrise$x
mf = merge(rise, mixrise, by.x = "CELLID", by.y = "id", all.x = T)
rise$mix = mf$x[match(as.character(rise$CELLID), mf$CELLID)]


describe(rise)



allrise = aggregate(events$RISE_TIME, 
                   list(id = events$CELLID, filename = events$FILENAME, protocol = events$PROTOCOL, recno = events$RECNO), mean, na.omit = T)
allrise$x = allrise$x 

relrise = data.frame(CELLID = unique(events$CELLID))
relrise$CULTURE = res$CULTURE[match(relrise$CELLID, res$CELLID)]
fields = c("nbqx", "sch", "sulp", "quin", "mix", "baclofen")
for (field in fields) {
  allrise[field] = NA
  fieldrise = allrise[grepl(field, allrise$protocol, ignore.case = T),]
  
  for (i in 1:nrow(fieldrise)) {
    id = fieldrise[i, "id"]
    recno = fieldrise[i, "recno"]
    index0 = which(allrise$id == id & allrise$recno == (recno - 1) & grepl("Base", allrise$protocol))
    index1 = which(allrise$id == id & allrise$recno == (recno + 1) & grepl("Base", allrise$protocol))
    
    if (length(index0) == 1 && length(index1) == 1) {
      reference = (allrise[index0, "x"] + allrise[index1, "x"]) / 2
    } else if (length(index0) == 1) {
      reference = allrise[index0, "x"] 
    } else if (length(index1) == 1) {
      reference = allrise[index1, "x"]
    }
    
    relrise[relrise$CELLID == id,field] = fieldrise[i, "x"] / reference
  }
}
describe(relrise)



# colMeans(rise[,3:7], na.rm =T)
# relrise = cbind(rise[c("CELLID", "CULTURE")], sweep(rise[c("mix", "sulp", "sch")], 1, rise$base, `/`))
# colMeans(relrise[, 3:5], na.rm = T)
dat = na.omit(melt(relrise[c("CELLID", "CULTURE", "mix", "sulp", "sch")]))

wilcox.test(relrise$mix, mu = 1.0)
wilcox.test(relrise$sulp, mu = 1.0)
  wilcox.test(relrise$sch, mu = 1.0)

dat$variable = as.character(dat$variable)
# dat$variable[dat$variable == "base"] = "DA"
# dat$variable[dat$variable == "nbqx"] = "DA NBQX"
dat$variable[dat$variable == "mix"] = "DA NBQX\nSCH23390\nSulpiride"
dat$variable[dat$variable == "sulp"] = "DA Sulpiride"
dat$variable[dat$variable == "sch"] = "DA SCH23390"
dat$variable = factor(dat$variable, levels = c("DA SCH23390", "DA Sulpiride", "DA NBQX\nSCH23390\nSulpiride"))

dftext = (dat %>% group_by(variable) %>% 
            summarise(n = sum(is.finite(value)), m = n_distinct(CULTURE[is.finite(value)])) %>% 
            mutate(label = sprintf("%d/%d", n, m)))[c("variable", "label")]
dftext$value = max(dat$value, na.rm = T) + diff(range(dat$value, na.rm = T)) / 20

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
p <- p + scale_color_manual(values = c("darkslateblue", "lightblue", "deeppink"))
p <- p + labs(fill = "", x = "", y = "Spontaneous PSC rise time / control")
# p <- p + scale_y_continuous(labels=scales::percent) + scale_x_log10()
p <- p + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p <- p + facet_grid(. ~ Group)

pdf(paste0("figure5e_", ".pdf"), width = 2.5, height = 3.5)
p
dev.off()



dat = melt(rise[c("CELLID", "CULTURE", "glut", "base")])
# 
# dat$variable = as.character(dat$variable)
# dat$variable[dat$variable == "base"] = "DA"
# dat$variable = factor(dat$variable)

dat$variable = dat$glut

dftext = (dat %>% group_by(variable) %>% 
            summarise(n = sum(is.finite(value)), m = n_distinct(CULTURE[is.finite(value)])) %>% 
            mutate(label = sprintf("%d/%d", n, m)))[c("variable", "label")]
dftext$value = max(dat$value, na.rm = T) + diff(range(dat$value, na.rm = T)) / 20

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
p <- p + labs(fill = "", x = "", y = "Spontaneous PSC rise time [ms]")
# p <- p + scale_y_continuous(labels=scales::percent) + scale_x_log10()
p <- p + scale_color_manual(values = c("turquoise", "orange")) 
# p <- p + scale_y_log10() + annotation_logticks(sides = "l")
p <- p + theme_classic() #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p <- p + facet_grid(. ~ Group)

pdf(paste0("figure5e__", ".pdf"), width = 2, height = 3.5)
p
dev.off()





chg = data.frame(CELLID = unique(events$CELLID))
chg$CULTURE = res$CULTURE[match(chg$CELLID, res$CELLID)]
chg$glut = res$glut[match(chg$CELLID, res$CELLID)]
chg = chg[!is.na(chg$glut),]

basechg = aggregate(base$CHARGE, list(id = base$CELLID), mean, na.omit = T)
basechg$x = basechg$x
mf = merge(chg, basechg, by.x = "CELLID", by.y = "id", all.x = T)
chg$base = mf$x[match(as.character(chg$CELLID), mf$CELLID)]* 1000

nbqxchg = aggregate(nbqx$CHARGE, list(id = nbqx$CELLID), mean, na.omit = T)
nbqxchg$x = nbqxchg$x
mf = merge(chg, nbqxchg, by.x = "CELLID", by.y = "id", all.x = T)
chg$nbqx = mf$x[match(as.character(chg$CELLID), mf$CELLID)]* 1000

sulpchg = aggregate(sulp$CHARGE, list(id = sulp$CELLID), mean, na.omit = T)
sulpchg$x = sulpchg$x
mf = merge(chg, sulpchg, by.x = "CELLID", by.y = "id", all.x = T)
chg$sulp = mf$x[match(as.character(chg$CELLID), mf$CELLID)]* 1000

schchg = aggregate(sch$CHARGE, list(id = sch$CELLID), mean, na.omit = T)
schchg$x = schchg$x 
mf = merge(chg, schchg, by.x = "CELLID", by.y = "id", all.x = T)
chg$sch = mf$x[match(as.character(chg$CELLID), mf$CELLID)]* 1000

mixchg = aggregate(mix$CHARGE, list(id = mix$CELLID), mean, na.omit = T)
mixchg$x = mixchg$x
mf = merge(chg, mixchg, by.x = "CELLID", by.y = "id", all.x = T)
chg$mix = mf$x[match(as.character(chg$CELLID), mf$CELLID)] * 1000

describe(chg)



allchg = aggregate(events$CHARGE, 
                    list(id = events$CELLID, filename = events$FILENAME, protocol = events$PROTOCOL, recno = events$RECNO), mean, na.omit = T)
allchg$x = allchg$x 

relchg = data.frame(CELLID = unique(events$CELLID))
relchg$CULTURE = res$CULTURE[match(relchg$CELLID, res$CELLID)]
fields = c("nbqx", "sch", "sulp", "quin", "mix", "baclofen")
for (field in fields) {
  allchg[field] = NA
  fieldchg = allchg[grepl(field, allchg$protocol, ignore.case = T),]
  
  for (i in 1:nrow(fieldchg)) {
    id = fieldchg[i, "id"]
    recno = fieldchg[i, "recno"]
    index0 = which(allchg$id == id & allchg$recno == (recno - 1) & grepl("Base", allchg$protocol))
    index1 = which(allchg$id == id & allchg$recno == (recno + 1) & grepl("Base", allchg$protocol))
    
    if (length(index0) == 1 && length(index1) == 1) {
      reference = (allchg[index0, "x"] + allchg[index1, "x"]) / 2
    } else if (length(index0) == 1) {
      reference = allchg[index0, "x"] 
    } else if (length(index1) == 1) {
      reference = allchg[index1, "x"]
    }
    
    relchg[relchg$CELLID == id,field] = fieldchg[i, "x"] / reference
  }
}
describe(relchg)

# colMeans(chg[,3:7], na.rm =T)
# relchg = cbind(chg[c("CELLID", "CULTURE")], sweep(chg[c("mix", "sulp", "sch")], 1, chg$base, `/`))
# colMeans(relchg[, 3:5], na.rm = T)
dat = na.omit(melt(relchg[c("CELLID", "CULTURE", "mix", "sulp", "sch")]))

wilcox.test(relchg$mix, mu = 1.0)
wilcox.test(relchg$sulp, mu = 1.0)
wilcox.test(relchg$sch, mu = 1.0)

dat$variable = as.character(dat$variable)
# dat$variable[dat$variable == "base"] = "DA"
# dat$variable[dat$variable == "nbqx"] = "DA NBQX"
dat$variable[dat$variable == "mix"] = "DA NBQX\nSCH23390\nSulpiride"
dat$variable[dat$variable == "sulp"] = "DA Sulpiride"
dat$variable[dat$variable == "sch"] = "DA SCH23390"
dat$variable = factor(dat$variable, levels = c("DA SCH23390", "DA Sulpiride", "DA NBQX\nSCH23390\nSulpiride"))

dftext = (dat %>% group_by(variable) %>% 
            summarise(n = sum(is.finite(value)), m = n_distinct(CULTURE[is.finite(value)])) %>% 
            mutate(label = sprintf("%d/%d", n, m)))[c("variable", "label")]
dftext$value = max(dat$value, na.rm = T) + diff(range(dat$value, na.rm = T)) / 20

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
p <- p + scale_color_manual(values = c("darkslateblue", "lightblue", "deeppink"))
p <- p + labs(fill = "", x = "", y = "Spontaneous PSC charge / control")
# p <- p + scale_y_continuous(labels=scales::percent) + scale_x_log10()
p <- p + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p <- p + facet_grid(. ~ Group)

pdf(paste0("figure5f_", ".pdf"), width = 2.5, height = 3.5)
p
dev.off()





chg$corr = chg$base - chg$nbqx
chg$corr[chg$corr < 0] = 0

mean(chg$base, na.rm = T)
library(psych)
describe(freq$corr)
mean(amp$corr, na.rm = T)
mean(rise$base, na.rm = T)

dat = melt(chg[c("CELLID", "CULTURE", "glut", "corr")])

# dat$variable = as.character(dat$variable)
# dat$variable[dat$variable == "corr"] = "DA"
# dat$variable = factor(dat$variable)

dat$variable = dat$glut

dftext = (dat %>% group_by(variable) %>% 
            summarise(n = sum(is.finite(value)), m = n_distinct(CULTURE[is.finite(value)])) %>% 
            mutate(label = sprintf("%d/%d", n, m)))[c("variable", "label")]
dftext$value = max(dat$value, na.rm = T) + diff(range(dat$value, na.rm = T)) / 20

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
p <- p + labs(fill = "", x = "", y = "Spontaneous PSC charge [fC]")
# p <- p + scale_y_continuous(labels=scales::percent) + scale_x_log10()
p <- p + scale_color_manual(values = c("turquoise", "orange")) 
# p <- p + scale_y_log10() + annotation_logticks(sides = "l")
p <- p + theme_classic() #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p <- p + facet_grid(. ~ Group)

pdf(paste0("figure5f__", ".pdf"), width = 2, height = 3.5)
p
dev.off()

relfreq[!is.na(relfreq$sch),]

rels = list(relfreq, relamp, relrise, relchg)
texts = c("frequency", "amplitude", "rise time", "charge")
letters = c("c", "b", "e", "f")
for (i in 1:length(texts)) {
  rel = rels[[i]]
  text = texts[i]
  letter = letters[i]
  
  dat = na.omit(melt(rel[c("CELLID", "CULTURE", "mix")]))
  dat = dat[dat$CELLID %in% res$CELLID[glut],]
  dat <- fix.names(dat)
  
  dftext = (dat %>% group_by(variable) %>% 
              summarise(n = sum(is.finite(value)), m = n_distinct(CULTURE[is.finite(value)])) %>% 
              mutate(label = sprintf("%d/%d", n, m)))[c("variable", "label")]
  dftext$value = max(dat$value, na.rm = T) + diff(range(dat$value, na.rm = T)) / 20
  
  p <- ggplot(dat, aes(x = variable, y = value, color = variable))
  p <- p + geom_quasirandom(show.legend = F)
  p <- p + stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult=1),
                        geom = "pointrange", color = "black",
                        position = position_dodge(0.8), show.legend = F)
  p <- p + geom_label(aes(label = label), 
                      data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, 
                      show.legend = F) + expand_limits(y = dftext$value[1])
  p <- p + geom_hline(yintercept = 1, color = "black")
  p <- p + scale_color_manual(values = c("deeppink", "green"))
  p <- p + labs(fill = "", x = "", y = paste0("Spontaneous PSC ", text, " / control"))
  p <- p + theme_classic() #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  pdf(paste0("figure5", letter, "_i", ".pdf"), width = 1.75, height = 3.5)
  print(p)
  dev.off()
  
  dat = na.omit(melt(rel[c("CELLID", "CULTURE", "sch")]))
  dat = dat[dat$CELLID %in% res$CELLID[glut],]
  dat <- fix.names(dat)
  
  dftext = (dat %>% group_by(variable) %>% 
              summarise(n = sum(is.finite(value)), m = n_distinct(CULTURE[is.finite(value)])) %>% 
              mutate(label = sprintf("%d/%d", n, m)))[c("variable", "label")]
  dftext$value = max(dat$value, na.rm = T) + diff(range(dat$value, na.rm = T)) / 20
  
  p <- ggplot(dat, aes(x = variable, y = value, color = variable))
  p <- p + geom_quasirandom(show.legend = F)
  p <- p + stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult=1),
                        geom = "pointrange", color = "black",
                        position = position_dodge(0.8), show.legend = F)
  p <- p + geom_label(aes(label = label), 
                      data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, 
                      show.legend = F) + expand_limits(y = dftext$value[1])
  p <- p + geom_hline(yintercept = 1, color = "black")
  p <- p + scale_color_manual(values = c("darkslateblue", "lightblue", "deeppink", "green"))
  p <- p + labs(fill = "", x = "", y = paste0("Spontaneous PSC ", text, " / control"))
  p <- p + theme_classic() #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  pdf(paste0("figure5", letter, "_ii", ".pdf"), width = 1.75, height = 3.5)
  print(p)
  dev.off()
  
  dat = na.omit(melt(rel[c("CELLID", "CULTURE", "sulp", "quin", "baclofen")]))
  dat = dat[dat$CELLID %in% res$CELLID[glut],]
  dat <- fix.names(dat)
  
  dftext = (dat %>% group_by(variable) %>% 
              summarise(n = sum(is.finite(value)), m = n_distinct(CULTURE[is.finite(value)])) %>% 
              mutate(label = sprintf("%d/%d", n, m)))[c("variable", "label")]
  dftext$value = max(dat$value, na.rm = T) + diff(range(dat$value, na.rm = T)) / 20
  
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
  p <- p + labs(fill = "", x = "", y = paste0("Spontaneous PSC ", text, " / control"))
  p <- p + theme_classic() #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  pdf(paste0("figure5", letter, "_iii", ".pdf"), width = 2.5, height = 3.5)
  print(p)
  dev.off()
  
}



wilcox.test(relamp$sch[relamp$CELLID %in% res$CELLID[glut]], mu = 1)
wilcox.test(relamp$sulp[relamp$CELLID %in% res$CELLID[glut]], mu = 1)
wilcox.test(relamp$quin[relamp$CELLID %in% res$CELLID[glut]], mu = 1)
wilcox.test(relamp$baclofen[relamp$CELLID %in% res$CELLID[glut]], mu = 1)
wilcox.test(relamp$mix[relamp$CELLID %in% res$CELLID[glut]], mu = 1)

wilcox.test(relchg$sch[relamp$CELLID %in% res$CELLID[glut]], mu = 1)
wilcox.test(relchg$sulp[relamp$CELLID %in% res$CELLID[glut]], mu = 1)
wilcox.test(relchg$quin[relamp$CELLID %in% res$CELLID[glut]], mu = 1)
wilcox.test(relchg$baclofen[relamp$CELLID %in% res$CELLID[glut]], mu = 1)
wilcox.test(relchg$mix[relamp$CELLID %in% res$CELLID[glut]], mu = 1)

wilcox.test(relrise$sch[relamp$CELLID %in% res$CELLID[glut]], mu = 1)
wilcox.test(relrise$sulp[relamp$CELLID %in% res$CELLID[glut]], mu = 1)
wilcox.test(relrise$quin[relamp$CELLID %in% res$CELLID[glut]], mu = 1)
wilcox.test(relrise$baclofen[relamp$CELLID %in% res$CELLID[glut]], mu = 1)
wilcox.test(relrise$mix[relamp$CELLID %in% res$CELLID[glut]], mu = 1)

wilcox.test(relfreq$sch[relamp$CELLID %in% res$CELLID[glut]], mu = 1)
wilcox.test(relfreq$sulp[relfreq$CELLID %in% res$CELLID[glut]], mu = 1)
wilcox.test(relfreq$quin[relamp$CELLID %in% res$CELLID[glut]], mu = 1)
wilcox.test(relfreq$baclofen[relamp$CELLID %in% res$CELLID[glut]], mu = 1)
wilcox.test(relfreq$mix[relamp$CELLID %in% res$CELLID[glut]], mu = 1)


describe(freq[freq$glut == "Glu","corr"])
describe(freq[freq$glut != "Glu","corr"])
wilcox.test(freq$corr ~ freq$glut)

describe(amp[amp$glut == "Glu","corr"])
describe(amp[amp$glut != "Glu","corr"])
wilcox.test(amp$corr ~ amp$glut)

describe(chg[chg$glut == "Glu","corr"])
describe(chg[chg$glut != "Glu","corr"])
wilcox.test(chg$corr ~ chg$glut)

describe(rise[rise$glut == "Glu","base"])
describe(rise[rise$glut != "Glu","base"])
wilcox.test(rise$base ~ rise$glut)

