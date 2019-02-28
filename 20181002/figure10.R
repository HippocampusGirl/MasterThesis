
# setwd("~/Dropbox/Univ/masterThesis/analysis/20181002/")
# setwd("/media/lw/leadisk/masterThesis/data/FACS/")
setwd("/media/lw/leadisk/masterThesis/analysis/20181002")

library(flowCore)
library(ggcyto)

# x <- read.FCS("../../data/FACS/Specimen_001_gfp.fcs")
x <- read.flowSet(files = c("../../data/FACS/Specimen_001_gfp.fcs", 
                            "../../data/FACS/Specimen_001_gfp_001.fcs"),
                  column.pattern = ".*(SC|PE-A|FITC-A).*")

for (i in 1:10) gc()

boundaries <- t(matrix(c(50000, 1000,
                       50000, 100000,
                       150000, 250000,
                       200000, 200000,
                       200000, 1000),
                     ncol = 5,nrow = 2))
colnames(boundaries) <- c("FSC-A","SSC-A")
gate1 <- polygonGate(filterId = "SC-A", boundaries)
gate2 <- rectangleGate(filterId = "SSC", "SSC-H" = c(0, 120000), "SSC-W" = c(30000, 170000))
gate3 <- rectangleGate(filterId = "FSC", "FSC-H" = c(0, 150000), "FSC-W" = c(50000, 160000))

fr <- Subset(x, gate1 %&% gate2 %&% gate3, c("PE-A", "FITC-A"))

for (i in 1:10) gc()

m = length(x)
n = nrow(x@frames$Specimen_001_gfp_001.fcs) +
  nrow(x@frames$Specimen_001_gfp.fcs)
dftext = data.frame("FSC-A" = c(1.2e5), "SSC-A" = c(2.5e5), label = sprintf("%d/%d", n, m), check.names = F)

p <- ggcyto(x, aes(x = "FSC-A", y = "SSC-A")) 
p <- p + geom_hex(bins = 256)
p <- p + geom_gate(gate1, colour = "coral") + geom_stats()
p <- p + geom_label(data = dftext, aes(label = label), show.legend = F)
# p <- p + scale_x_log10() + scale_y_log10()
p <- p + theme_classic() 
p <- p + theme(strip.background = element_blank(), strip.text = element_blank())
p <- p + facet_null()

pdf("figure10a.pdf", width = 3, height = 3.5)
p
dev.off()

kf <- kmeansFilter("FITC-A" = c("Pop1", "Pop2"), filterId = "kmFilter")
kfr <- filter(fr, kf)
summary(kfr)

pop1 = split(fr, kfr, population = "Pop1")[[1]]
n1 = fsApply(pop1, nrow)
pop2 = split(fr, kfr, population = "Pop2")[[1]]
n2 = fsApply(pop2, nrow)

pop = pop2
if (sum(n1) < sum(n2)) pop = pop1
fsApply(pop, nrow)

popmean = apply(fsApply(pop, exprs), 2, mean) * 1.6
popcov = as.matrix(cov(fsApply(pop, exprs))) * 10

gate4 = ellipsoidGate(filterId = "Pos", mean = popmean, cov = popcov)

m = length(fr)
n = nrow(fr@frames$Specimen_001_gfp_001.fcs) +
      nrow(fr@frames$Specimen_001_gfp.fcs)
dftext = data.frame("PE-A" = c(100), "FITC-A" = c(30000), label = sprintf("%d/%d", n, m), check.names = F)

sum(is.na(fsApply(fr, exprs)))

p <- ggcyto(fr, aes(x = "PE-A", y = "FITC-A")) 
p <- p + geom_hex(bins = 256) 
p <- p + geom_gate(gate4, colour = "coral") + geom_stats(hjust = 0, vjust = 1, digits = 2)
p <- p + geom_label(data = dftext, aes(label = label), show.legend = F)
p <- p + scale_x_logicle() + scale_y_logicle()
p <- p + theme_classic() 
p <- p + theme(strip.background = element_blank(), strip.text = element_blank())
p <- p + facet_null()
p <- p + ggtitle("")

pdf("figure10b.pdf", width = 3, height = 3.5)
p
dev.off()
