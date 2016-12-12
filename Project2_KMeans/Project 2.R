mydata <- read.table(file="C:/Users/talem/workspace_Mars/Project2_KMeans/KMeansPCA_new_dataset_1.csv", header=FALSE,row.names=1, sep=",")
fit <- princomp(mydata)
plot(fit,type="lines")
biplot(fit)


mydata <- read.table(file="C:/Users/talem/workspace_Mars/Project2_KMeans/AlgmClusters_new_dataset_2.csv", header=FALSE,row.names=1, sep=",")
fit <- princomp(mydata)
plot(fit,type="lines")
biplot(fit)


mydata <- read.table(file="C:/Users/talem/workspace_Mars/Project2_KMeans/DBScan_new_dataset_1.csv", header=FALSE,row.names=1, sep=",")
fit <- princomp(mydata)
plot(fit,type="lines")
biplot(fit)


library(devtools)
install_github("ggbiplot", "vqv")

library(ggbiplot)
g <- ggbiplot(fit, obs.scale = 1, var.scale = 1, 
              ellipse = TRUE, 
              circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)


install.packages('ggfortify')
library(ggfortify)
autoplot(fit)
