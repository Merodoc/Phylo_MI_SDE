library(phytools)

data(cordylid.tree)
data(cordylid.data)

cordylid.armor_score <- setNames(cordylid.data$pPC1, rownames(cordylid.data))

cordylid.mcmc <- anc.Bayes(cordylid.tree, cordylid.armor_score, ngen = 500000)
cordylid.ace <- summary(cordylid.mcmc)

cordylid.contMap <- contMap(cordylid.tree, cordylid.armor_score, anc.states = cordylid.ace, plot = FALSE)
cordylid.contMap <- setMap(cordylid.contMap, viridisLite::viridis(n=10, direction = 1))
plot(cordylid.contMap, ftype = "i", fsize = c(0.6, 0.7), leg.txt = "PC 1 (increasing armor)", lwd =3)
nodelabels(frame = "circle", bg = "white", cex = 0.6)

contMap(cordylid.tree, cordylid.armor_score)

data(tropidurid.tree)
data(tropidurid.data)

print(tropidurid.tree, printlen = 2)

tropidurid.tree <- mergeMappedStates(tropidurid.tree, "n_rock", "non-rock dwelling")
tropidurid.tree <- mergeMappedStates(tropidurid.tree, "rock", "rock-dwelling")

cols <- setNames(c("white", "black"), c("non-rock dwelling", "rock-dwelling"))

sigmoidPhylogram(tropidurid.tree, direction = "upwards", outline = TRUE, colors = cols, direction = "upwards",
                 outline = TRUE, lwd = 2, fsize = 0.4, ftype = "i", offset = 1)

legend("bottomright", c("non-rock dwelling", "rock-dwelling"), pch = 22, pt.bg = cols, cex = 0.8, pt.cex = 1.2)

tropidurid.fits <- evolvcv.lite(tropidurid.tree, tropidurid.data)

anova(tropidurid.fits)

data(primate.tree)
data(primate.data)
primate.lnSkull <- setNames(log(primate.data$Skull_length), rownames(primate.data))

par(mfrow = c(1,1))
primate.widthMap <- edge.widthMap(primate.tree, primate.lnSkull)
plot(primate.widthMap, color = palette()[4], 
     legend = "log(skull length)", border = TRUE, fsize = 0.4, mar = c(4.1, 1.1, 2.1, 0.1))
mtext("a)", adj = 0, line = 0, cex = 1.4)
phenogram(primate.tree, primate.lnSkull, fsize = 0.4, ftype = "i",
          spread.cost = c(1, 0), mar = c(4.1, 4.1, 2.1, 0.1),
          quiet = TRUE, las = 1, cex.axis = 0.8)



mars_tree <- read.tree("newtree.nwk")
Mars <- read.csv("mars.csv", header=TRUE)
Mars_Avg <- read.csv("mars_avg.csv")

mars_avg2 <- Mars_Avg[which(Mars_Avg$X %in% mars_tree$tip.label),]
rownames(mars_avg2) <- mars_avg2$X
mars_avg2 <- mars_avg2[match(mars_tree$tip.label, rownames(mars_avg2)),]

mars.data <- as.matrix(mars_avg2$femur)
mars.data
mars.data <- t(mars.data)
mars.data <- setNames(mars.data, mars_avg2[,1])
mars.data
names(mars.data) <- mars_avg2[,1]
mars.data
fit <- fastAnc(mars_tree, mars.data, vars = TRUE, CI = TRUE)
fit

par(mfrow = c(1,1))
mars.contMap <- contMap(mars_tree, mars.data, anc.states = fit$ace, plot = FALSE)
mars.contMap <- setMap(mars.contMap, viridisLite::viridis(n=10, direction = 1))
plot(mars.contMap, ftype = "i", fsize = c(0.6, 0.7), leg.txt = "Trochanteric Fossa Length", lwd =3)

mars.mcmc <- anc.Bayes(mars_tree, mars.data, ngen = 500000)




fit$ace
mars.data
mars_tree

plotTree(mars_tree)

bbridge <- read.csv("FemurMeans.csv", header = TRUE)

mars.contMap <- contMap(mars_tree, bbridge.leaves, anc.states = bbridge.anc, plot = FALSE)
mars.contMap <- setMap(mars.contMap, viridisLite::viridis(n=10, direction = 1))
plot(mars.contMap, ftype = "i", fsize = c(0.6, 0.7), leg.txt = "Trochanteric Fossa Length", lwd =3)

fit <- fastAnc(mars_tree, bbridge.leaves, vars = TRUE, CI = TRUE)
fit

par(mfrow = c(1,1))
mars.contMap <- contMap(mars_tree, bbridge.leaves, anc.states = fit$ace, plot = FALSE)
mars.contMap <- setMap(mars.contMap, viridisLite::viridis(n=10, direction = 1))
plot(mars.contMap, ftype = "i", fsize = c(0.6, 0.7), leg.txt = "Trochanteric Fossa Length", lwd =3)


bbridge
bbridge.leaves <- bbridge$femur[1:35]
bbridge.leaves
names(bbridge.leaves) <- bbridge$Species[1:35]

bbridge.anc <- bbridge$femur[36:69] 
names(bbridge.anc) <- bbridge$Species[36:69]
