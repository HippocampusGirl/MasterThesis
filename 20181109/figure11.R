setwd("~/Dropbox/Univ/masterThesis/analysis/20181109/")

library(ggplot2)
library(ggrepel)
library(reshape2)
library(MASS)
library(scales)

egfp = read.table("egfpem.txt", header = F)
names(egfp) = c("Wavelength", "Density")

egfp = rbind(egfp, data.frame(Wavelength = c(350, 700, egfp$Wavelength[1]-(1:50)), Density = c(0, 0, exp(-(1:50) / 5) * egfp$Density[1])))

filters = data.frame(
  label = factor(c("SC filter", "FITC filter", "PE filter", "Laser", "eGFP emission"), levels = c("SC filter", "FITC filter", "PE filter", "Laser", "eGFP emission")),
  xmin = c(488 - 5, 530 - 15, 575 - 13, NA, NA),
  xmax = c(488 + 5, 530 + 15, 575 + 13, NA, NA),
  Wavelength = c(488, 530, 575, 488, 518), 
  Density = c(1, 0.9, 0.8, 0.7, 0.6),
  ymin = rep(0, 5), ymax = rep(1, 5)
)

p <- ggplot(egfp, aes(x = Wavelength, y = Density)) 
p <- p + geom_line(color = "green") + geom_area(fill = "green", position = position_identity(), alpha = 0.5)
p <- p + geom_rect(data = filters, 
          aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf, fill = label), alpha = 0.5,
          show.legend = F)
p <- p + geom_vline(xintercept = 488, color = "blue")
p <- p + geom_label_repel(data = filters, aes(label = label, color = label), show.legend = F, seed = 1001)
p <- p + labs(fill = "", x = "Wavelength [nm]", y = "Relative intensity")
p <- p + scale_color_manual(values = c("gray50", "gray50", "gray50", "blue", "green"))
p <- p + scale_fill_manual(values = c("gray50", "gray50", "gray50", "blue", "green"))
p <- p + theme_classic() 
p <- p + theme(strip.background = element_blank(), strip.text = element_blank())
p <- p + facet_null()

pdf("figure11a.pdf", width = 3, height = 3.5)
p
dev.off()

ex = data.frame(Wavelength = 200:800)
for (mol in c("alexa350", "alexa405", "egfp", "alexa555", "alexa647", "fm464")) {
  data = read.table(sprintf("%sex.txt", mol), header = F)
  names(data) = c("Wavelength", "Density")
  data = rbind(data, data.frame(Wavelength = c(200, 800), Density = c(0, 0)))
  ex[mol] = approx(x = data$Wavelength, y = data$Density, xout = ex$Wavelength)$y
}

levels = c("350 nm filter", "403 nm filter", "490 nm filter", "555 nm filter", "635 nm filter", 
           "DyLight 350", "Alexa 405", "eGFP", "Alexa 555", "Alexa 647", "FM4-64FX")

dat = melt(ex, id = "Wavelength")

dat$variable = as.character(dat$variable)
dat$variable[dat$variable == "alexa350"] = "DyLight 350"
dat$variable[dat$variable == "alexa405"] = "Alexa 405"
dat$variable[dat$variable == "egfp"] = "eGFP"
dat$variable[dat$variable == "alexa555"] = "Alexa 555"
dat$variable[dat$variable == "alexa647"] = "Alexa 647"
dat$variable[dat$variable == "fm464"] = "FM4-64FX"
dat$variable = factor(dat$variable, 
                      levels = levels)

exfilters = data.frame(
  variable = factor(c("350 nm filter", "403 nm filter", "490 nm filter", "555 nm filter", "635 nm filter", 
                   "DyLight 350", "Alexa 405", "eGFP", "Alexa 555", "Alexa 647", "FM4-64FX"), 
                 levels = levels),
  xmin = c(350 - 25, 403 - 6, 490 - 10, 555 - 14, 635 - 10, NA, NA, NA, NA, NA, NA),
  xmax = c(350 + 25, 403 + 6, 490 + 10, 555 + 14, 635 + 10, NA, NA, NA, NA, NA, NA),
  Wavelength = c(350, 403, 490, 555, 635, 350, 405, 488, 555, 647, 506),
  value = c(1, 0.9, 1, 0.8, 1, 0.5, 0.25, 0.25, 0.25, 0.25, 0.25),
  ymin = rep(0, 11), ymax = rep(1, 11)
)

# palette = c("gray40", "gray40", "gray40", "gray40", "gray40", 
#             "darkorchid", "deepskyblue", "green", "darkorange", "deeppink", "yellow")
# palette_ = c("gray40", "gray40", "gray40", "gray40", "gray40", 
#             "darkorchid", "deepskyblue", "darkorange", "deeppink", "green", "yellow")

palette = c("gray50", "gray50", "gray50", "gray50", "gray50", 
            "darkslateblue", "darkslateblue", "darkslateblue", "darkslateblue", "darkslateblue", "deeppink")
palette_ = palette

p <- ggplot(dat, aes(x = Wavelength, y = value, fill = variable, color = variable)) 
p <- p + ggtitle("Excitation")
p <- p + geom_rect(data = exfilters, color = "transparent",
                   aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf, fill = variable), alpha = 0.5,
                   show.legend = F)
p <- p + geom_area(alpha = 0.3, position = position_identity(), show.legend = F)
p <- p + geom_label_repel(data = exfilters, aes(label = variable), show.legend = F, seed = 1001, fill = "white", hjust = 0.5)
p <- p + scale_color_manual(values = palette_)
p <- p + scale_fill_manual(values = palette)
p <- p + labs(fill = "", x = "Wavelength [nm]", y = "Relative intensity")
p <- p + theme_classic() 
p <- p + theme(strip.background = element_blank(), strip.text = element_blank())
p <- p + facet_null()

pdf("figure11b.pdf", width = 5, height = 3.5)
p
dev.off()







em = data.frame(Wavelength = 200:800)
for (mol in c("alexa350", "alexa405", "egfp", "alexa555", "alexa647", "fm464")) {
  data = read.table(sprintf("%sem.txt", mol), header = F)
  names(data) = c("Wavelength", "Density")
  data = rbind(data, data.frame(Wavelength = c(200, 800), Density = c(0, 0)))
  em[mol] = approx(x = data$Wavelength, y = data$Density, xout = em$Wavelength)$y
}

apply(em, 2, max)

dat = melt(em, id = "Wavelength")

levels = c("457 nm filter", "528 nm filter", "617 nm filter", "685 nm filter", 
           "DyLight 350", "Alexa 405", "eGFP", "Alexa 555", "Alexa 647", "FM4-64FX")

dat$variable = as.character(dat$variable)
dat$variable[dat$variable == "alexa350"] = "DyLight 350"
dat$variable[dat$variable == "alexa405"] = "Alexa 405"
dat$variable[dat$variable == "egfp"] = "eGFP"
dat$variable[dat$variable == "alexa555"] = "Alexa 555"
dat$variable[dat$variable == "alexa647"] = "Alexa 647"
dat$variable[dat$variable == "fm464"] = "FM4-64FX"
dat$variable = factor(dat$variable, 
                      levels = levels)

emfilters = data.frame(
  variable = factor(levels, 
                    levels = levels),
  xmin = c(457 - 5, 528 - 19, 617 - 36.5, 685 - 20, NA, NA, NA, NA, NA, NA),
  xmax = c(457 + 5, 528 + 19, 617 + 36.5, 685 + 20, NA, NA, NA, NA, NA, NA),
  Wavelength = c(457, 528, 617, 685, 420, 445, 508, 567, 670, 720),
  value = c(1, 1, 1, 1, 0.25, 0.25, 0.25, 0.25, 0.25, 0.4),
  ymin = rep(0, 10), ymax = rep(1, 10)
)

# palette = c("gray50", "gray50", "gray50", "gray50", 
#             "darkorchid", "deepskyblue", "green", "darkorange", "deeppink", "yellow")
# palette_ = c("gray50", "gray50", "gray50", "gray50", 
#             "darkorchid", "deepskyblue", "darkorange", "deeppink", "green", "yellow")

palette = c("gray50", "gray50", "gray50", "gray50", 
            "darkslateblue", "darkslateblue", "darkslateblue", "darkslateblue", "darkslateblue", "deeppink")
palette_ = palette

p <- ggplot(dat, aes(x = Wavelength, y = value, fill = variable, colour = variable)) 
p <- p + ggtitle("Emission")
p <- p + geom_rect(data = emfilters, color = "transparent",
                   aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf, fill = variable), alpha = 0.5,
                   show.legend = F)
p <- p + geom_area(alpha = 0.3, position = position_identity(), show.legend = F)
p <- p + geom_label_repel(data = emfilters, aes(label = variable), 
                          show.legend = F, seed = 1001, fill = "white", hjust = 0.5)
p <- p + scale_fill_manual(values = palette, breaks = levels) + scale_colour_manual(values = palette_, breaks = levels) 
p <- p + labs(fill = "", x = "Wavelength [nm]", y = "Relative intensity")
p <- p + theme_classic() 
p <- p + theme(strip.background = element_blank(), strip.text = element_blank())
p <- p + facet_null()

pdf("figure11c.pdf", width = 5, height = 3.5)
p
dev.off()





mix = as.matrix(read.table("A11HB_16_e.txt", header = F))

# mix = mix[2:nrow(mix),]

fluorophore_names = c("Background", "DyLight 350", "Alexa 405", "eGFP", "Alexa 555", "Alexa 647", "FM4-64FX")
filterset_names = c("350 nm/457 nm", "403 nm/457 nm", "490 nm/528 nm", "555 nm/617 nm", "635 nm/685 nm", "555 nm/685 nm")

rownames(mix) = fluorophore_names
colnames(mix) = filterset_names

# mix = t(solve(unmix))

mixNorm = sweep(mix, 1, rowMeans(mix), `/`)

dat = melt(mixNorm)
names(dat) = c("Fluorochrome", "Filter setting [Excitation/Emission]", "Value [AU]")

p <- ggplot(dat, aes(x = Fluorochrome, y = `Filter setting [Excitation/Emission]`, fill = `Value [AU]`)) 
p <- p + ggtitle("Empiricial mixing matrix")
p <- p + geom_tile()
p <- p + theme_minimal() 
p <- p + theme(strip.background = element_blank(), strip.text = element_blank(),
               axis.text.x = element_text(angle = 45, hjust = 1),
               axis.text.y = element_text(angle = 45, hjust = 1))
# p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "white",
#                      midpoint = 0, limit = c(-1e4, 1e4), space = "Lab")
# p <- p + scale_fill_gradientn(colours = c("blue","cyan","white", "yellow","red"), limits = c(-10, 10),
#                               values = rescale(c(-1, -0.1, 0, 0.1, 1)))
p <- p + scale_fill_gradientn(colours = c("white", "yellow","red"), limits = c(0, 6),
                              values = rescale(c(0, 0.05, 1)))
p <- p + facet_null()

pdf("figure11d.pdf", width = 5.5, height = 3.5)
p
dev.off()


fluorophore_names = c("DyLight 350", "Alexa 405", "eGFP", "Alexa 555", "Alexa 647", "FM4-64FX")

emfiltersm = matrix(T, ncol = nrow(em), nrow = 4)
for (i in 1:4) {
  data = read.table(sprintf("%dnm.txt", emfilters[i,"Wavelength"]), header = F)
  names(data) = c("Wavelength", "Density")
  data = rbind(data, data.frame(Wavelength = c(200, 800), Density = c(0, 0)))
  emfiltersm[i, ] = approx(x = data$Wavelength, y = data$Density, xout = em$Wavelength)$y
  # emfiltersm[i, 1:(emfilters[i,"xmin"]-200)] = F
  # emfiltersm[i, round(emfilters[i,"xmax"]-200):ncol(emfiltersm)] = F
}
rownames(emfiltersm) = emfilters$variable[1:4]

exfiltersm = matrix(T, ncol = nrow(ex), nrow = 5)
for (i in 1:5) {
  data = read.table(sprintf("%dnm.txt", exfilters[i,"Wavelength"]), header = F)
  names(data) = c("Wavelength", "Density")
  data = rbind(data, data.frame(Wavelength = c(200, 800), Density = c(0, 0)))
  exfiltersm[i, ] = approx(x = data$Wavelength, y = data$Density, xout = em$Wavelength)$y
  # exfiltersm[i, 1:(exfilters[i,"xmin"]-200)] = F
  # exfiltersm[i, round(exfilters[i,"xmax"]-200):ncol(exfiltersm)] = F
}
rownames(exfiltersm) = exfilters$variable[1:5]

data = read.table("quad.txt", header = F)
names(data) = c("Wavelength", "Density")
data = rbind(data, data.frame(Wavelength = c(200, 800), Density = c(0, 0)))
quad = approx(x = data$Wavelength, y = data$Density, xout = em$Wavelength)$y
plot(quad, type = "l")

exfiltersquad = sweep(exfiltersm, 2, 1 - quad, `*`)
emfiltersquad = sweep(emfiltersm, 2, quad, `*`)

a = (emfiltersquad %*% as.matrix(em[,2:ncol(em)]))
b = (exfiltersquad %*% as.matrix(ex[,2:ncol(ex)]))

filtersets = list(c(1, 1), c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(3, 5))
names(filtersets) = filterset_names
w = matrix(NA, ncol = length(filtersets), nrow = ncol(a))
rownames(w) = filterset_names
colnames(w) = fluorophore_names
for (k in 1:length(filtersets)) {
  cc = matrix(c(0, 0, 0, 0, 
                0, 0, 0, 0,  
                0, 0, 0, 0,  
                0, 0, 0, 0,
                0, 0, 0, 0), nrow = nrow(emfiltersm), ncol = nrow(exfiltersm))
  cc[filtersets[[k]][1], filtersets[[k]][2]] = 1
  rownames(cc) = rownames(emfiltersm)
  colnames(cc) = rownames(exfiltersm)
  w[, k] = diag(t(a) %*% cc %*% b)
}
unmix = solve(w)

mix = t(w)
mixNorm = sweep(mix, 1, rowMeans(mix), `/`)

write.table(x = unmix * 4096, file = "A_theoretical.txt", row.names = F, col.names = F)
write.table(x = mix, file = "Ainv_theoretical.txt", row.names = F, col.names = F)

dat = melt(mixNorm)

names(dat) = c("Fluorochrome", "Filter setting [Excitation/Emission]", "Value [AU]")

p <- ggplot(dat, aes(x = Fluorochrome, y = `Filter setting [Excitation/Emission]`, fill = `Value [AU]`)) 
p <- p + ggtitle("Mixing matrix")
p <- p + geom_tile()
p <- p + theme_minimal() 
p <- p + theme(strip.background = element_blank(), strip.text = element_blank(),
               axis.text.x = element_text(angle = 45, hjust = 1),
               axis.text.y = element_text(angle = 45, hjust = 1))
# p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "white",
#                      midpoint = 0, space = "Lab")
p <- p + scale_fill_gradientn(colours = c("white", "yellow","red"), limits = c(0, 6),
                              values = rescale(c(0, 0.05, 1)))
p <- p + facet_null()

pdf("figure11e.pdf", width = 5, height = 3.5)
p
dev.off()









em_ = cbind(em, 0.1)
colnames(em_)[ncol(em_)] = "Background"
ex_ = cbind(ex, 0.1)
colnames(ex_)[ncol(ex_)] = "Background"

a = (emfiltersquad %*% as.matrix(em_[,2:ncol(em_)]))
b = (exfiltersquad %*% as.matrix(ex_[,2:ncol(ex_)]))

filtersets = list(c(1, 1), c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(3, 5))
names(filtersets) = filterset_names
w = matrix(NA, nrow = length(filtersets), ncol = ncol(a))
rownames(w) = filterset_names
colnames(w) = c(fluorophore_names, "Background")
for (k in 1:length(filtersets)) {
  cc = matrix(c(0, 0, 0, 0, 
                0, 0, 0, 0,  
                0, 0, 0, 0,  
                0, 0, 0, 0,
                0, 0, 0, 0), nrow = nrow(emfiltersm), ncol = nrow(exfiltersm))
  cc[filtersets[[k]][1], filtersets[[k]][2]] = 1
  rownames(cc) = rownames(emfiltersm)
  colnames(cc) = rownames(exfiltersm)
  w[k, ] = diag(t(a) %*% cc %*% b)
}
unmix = ginv(w)
rownames(unmix) = colnames(w)
colnames(unmix) = rownames(w)

unmix * 4096

# unmixNorm = t(sweep(unmix, 1, abs(rowSums(unmix)), `/`))
