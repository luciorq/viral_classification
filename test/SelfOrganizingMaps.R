########################################
#SELF-ORGANIZING MAPS
# Change the data frame with training data to a matrix
# Also center and scale all variables to give them equal importance during
# the SOM training process. 
db_matrix <- as.matrix(scale(db[,212:291]))
head(db_matrix)

# Create the SOM Grid - you generally have to specify the size of the 
# training grid prior to training the SOM. Hexagonal and Circular 
# topologies are possible
som_grid <- somgrid(xdim = 20, ydim=20, topo="hexagonal")

# Finally, train the SOM, options for the number of iterations,
# the learning rates, and the neighbourhood are available
som_model <- som(db_matrix, 
                 grid=som_grid, 
                 rlen=100, 
                 alpha=c(0.05,0.01), 
                 keep.data = TRUE,
                 n.hood= "circular" )
#Nodes count
plot(som_model, type="count")
#Neighbour Distance
plot(som_model, type="dist.neighbours")
#Heatmaps
coolBlueHotRed <- function(n, alpha = 1) {rainbow(n, end=4/6, alpha=alpha)[n:1]}
plot(som_model, type = "property", property = som_model$codes[,4], main=names(som_model$data)[4], palette.name=coolBlueHotRed)


##################
#Clustering
mydata <- som_model$codes 
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var)) 
for (i in 2:15) {
  wss[i] <- sum(kmeans(mydata, centers=i)$withinss)
}
plot(wss)
## use hierarchical clustering to cluster the codebook vectors
som_cluster <- cutree(hclust(dist(som_model$codes)), 6)

source('coolBlueHotRed.R')
# Colour palette definition
pretty_palette <- c("#1f77b4", '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2')
# plot these results:
plot(som_model, type="mapping", bgcol = pretty_palette[som_cluster], main = "Clusters") 
add.cluster.boundaries(som_model, som_cluster)
