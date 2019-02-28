setwd("~/Dropbox/Univ/masterThesis/analysis/20181107/")

library(ggplot2)
library(dplyr)

gfp = read.table("ResultsParticles_C=0.xls", header = T)
th = read.table("ResultsParticles_C=2.xls", header = T)

merge = merge(x = gfp, y = th, by = c("File", "Particle"), suffixes = c(".gfp", ".th"))

merge$Culture = factor(gsub("mass(\\d).*_.*", "\\1", merge$File))
merge$GDNF = grepl("GDNF", merge$File)

merge = merge[merge$Area.gfp < 250,]
hist(merge$Area.gfp)

nogdnf = merge[merge$GDNF == F,]

# library("mixtools")

# hist(nogdnf$Mean.th)
# mix <- normalmixEM(nogdnf$Mean.th, k = 2)
# lik1 <- function(x) {
#   lam = mix$lambda[1]
#   mu = mix$mu[1]
#   sigma = mix$sigma[1]
#   dnorm(x, mu, sigma)
# }
# lik2 <- function(x) {
#   lam = mix$lambda[2]
#   mu = mix$mu[2]
#   sigma = mix$sigma[2]
#   dnorm(x, mu, sigma)
# }
# dl <- function (x) {
#   lik1(x)-lik2(x)
# }
# uniroot(dl, c(0, 4000))



dftext = (nogdnf %>% 
            summarise(n = n(), m = n_distinct(Culture)) %>% 
            mutate(label = sprintf("%d/%d", n, m)))[c("label")]
dftext$Mean.gfp = 40
dftext$Mean.th = max(nogdnf$Mean.th)

boundary = 200

mean(nogdnf$Mean.th > boundary)

p <- ggplot(nogdnf, aes(x = Mean.gfp, y = Mean.th, color = "turquoise"))
p <- p + geom_point(show.legend = F)
p <- p + geom_hline(yintercept = boundary, color = "black")
p <- p + scale_x_log10() + scale_y_log10() + annotation_logticks(sides = "bl")
p <- p + geom_label(aes(label = label), 
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, color = "black",
                    show.legend = F) + expand_limits(y = dftext$value[1])
p <- p + theme_classic() 
p <- p + theme(legend.margin = margin(t = 0, unit = "cm"))
p <- p + scale_color_manual(values = c("turquoise", "turquoise4", "slateblue1"))
p <- p + labs(fill = "", x = "GFP fluorescence [AU]", y = "TH fluorescence [AU]")

pdf(paste0("figure9a", ".pdf"), width = 2.5, height = 3.5)
p
dev.off()

# hist(nogdnf$Mean.th, xlim = c(0, 1500), ylim = c(0, 100))
# 
# labels = nogdnf$Mean.th > 1000
# scores = nogdnf$Mean.gfp
# 
# ord = order(scores, decreasing = T)
# scores = scores[ord]
# labels = labels[ord]
# tpr = cumsum(labels)/sum(labels)
# fpr = cumsum(!labels)/sum(!labels)
# 
# cutoff = scores[which(fpr < 0.05)]
# cutoff = cutoff[length(cutoff)]
# mean(labels[scores >= cutoff])
# mean(labels[scores < cutoff])

# p <- ggplot(nogdnf, aes(x = fpr, y = tpr))
# p <- p + geom_line()
# p <- p + theme_classic() 
# p <- p + labs(fill = "", x = "1 - FPR", y = "TPR")
# # p <- p + scale_x_continuous(labels=scales::percent) + scale_y_continuous(labels=scales::percent)
# p

gdnf = merge[merge$Culture == 4,]

dftext = (gdnf %>% 
            summarise(n = n(), m = n_distinct(Culture)) %>% 
            mutate(label = sprintf("%d/%d", n, m)))[c("label")]
dftext$Mean.gfp = 100
dftext$Mean.th = max(gdnf$Mean.th)

boundary = 300

gdnf$GDNF = factor(gdnf$GDNF, labels = c("Control", "10 ng/ml GDNF"))
p <- ggplot(gdnf, aes(x = Mean.gfp, y = Mean.th, color = GDNF))
p <- p + geom_point(show.legend = F)
p <- p + geom_hline(yintercept = boundary, color = "black")
# p <- p + geom_hex(bins = 64)
p <- p + scale_x_log10() + scale_y_log10() + annotation_logticks(sides = "bl")
p <- p + geom_label(aes(label = label), 
                    data = dftext, stat = "identity", position = position_dodge(0.9), alpha = 0.5, color = "black",
                    show.legend = F) + expand_limits(y = dftext$value[1])
p <- p + geom_label(data = data.frame(GDNF = c("Control", "10 ng/ml GDNF"), Mean.gfp = c(300, 55), Mean.th = c(1000, 1800)), 
                   aes(label = GDNF), show.legend = F) 
p <- p + scale_color_manual(values = c("turquoise", "orchid"))
p <- p + theme_classic() 
p <- p + theme(legend.margin = margin(t = 0, unit = "cm"))
p <- p + labs(fill = "", x = "GFP fluorescence [AU]", y = "TH fluorescence [AU]")

pdf(paste0("figure9b", ".pdf"), width = 3, height = 3.5)
p
dev.off()

library(psych)

corr.test(merge$Mean.gfp, merge$Mean.th)
corr.test(gdnf$Mean.gfp, gdnf$Mean.th)

# p

