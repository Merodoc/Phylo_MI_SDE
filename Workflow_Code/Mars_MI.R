library("mice")
library("plyr")
Mars <- read.csv("mars.csv", header=TRUE)
Mars <- Mars[,-1]
Mars$Species <- as.factor(Mars$Species)
md.pattern(Mars)

imp <- mice(Mars, m=25)

fit <- lm.mids(dentary ~ dentary.1, imp)

imp.dat <- complete(imp, "long")

imp.dat$dentary

MEANS <- ddply(imp.dat, .(Species), summarize, mean=mean(dentary))

write.csv(MEANS, "Imp_Mars25.csv", row.names=FALSE)

library(ape)

mars_tree <- read.tree("MarS_TimeTree.nwk")

