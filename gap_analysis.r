##################################################################################################### 
## Assessing the performance of protected area networks in representing a regional species pool            
##                                                 
## It computes the number of taxa represented in a reserve network for a given conservation        
## target and assesses whether this level of representation is significantly lower or greater than    
## expected by chance. An analysis for multiple thresholds to consider a grid cell as protected is 
## undertaken                                      
##                                                 
## Pedro Abellán                                   
## July 2014                                       
##################################################################################################### 

### Load data 
dat <- read.table("matrix.txt", header=TRUE, row.names=1)  # load presence/absence matrix; cells in rows and species in columns; it must contain row and columns names 
cells <-  read.table("cells.txt", header=TRUE, row.names=1) # load table with information about grid cells; it must contain row and columns names and a field with the percentage of overlap with protected areas ("protected" in this example) 

### Set parameters for analysis 
pro <- "protected" # name of the field in table "cells" with the percentage of overlap with protected areas 
occ <- 1 # set conservation target (number of occurrences within the network to consider a species as covered) 
runs <- 999 # set number of runs for randomizations
thresholds <- c(1:100) # set the set of thresholds (percentage of overlap with reserves) to consider a grid cell as protected. In this example, we use 100 possible thresholds (from 1 to 100% in increments of 1) 

## function to compute species richness in an equivalent number of random cells 
# inputs: data matrix, number of protected cells and conservation target 
sr.rand.occ <- function(dat, n.pro, occ){ 
  sam <- sample(c(1:nrow(dat)), n.pro) # sample an equivalent number of random cells from data matrrix 
  sub.ran <- dat[sam,] # submatrix with just random cells 
  spp.ran <- colSums(sub.ran) # the number of cells in which each species occur 
  occ.ran <- spp.ran # duplicate spp.ram to modificate it 
  occ.ran[] <- 0 
  occ.ran[which(spp.ran>=occ)] <- 1  # identify those species that meet conservation target 
  sr.random <- length(which(occ.ran > 0)) # number of species that meet conservation target 
  return(sr.random) 
} 

### Analyses for different thresholds to consider a grid cell as protected 
# prepare outputs 
out.null.test <-  data.frame() 
out.occurrences <-  data.frame(taxa=names(dat)) 

for (i in thresholds){ # loop for different thresholds to considered a cell as protected
threshold <- i  # threshold ith 

# which and how many cells are protected under that threshold?
protected <- which(cells[,pro]>=threshold)  # identify protected cells
n.pro <- length(protected) # number of protected cells

# Compute total species richness included in protected cells
sub.pro <- dat[protected,] # submatrix with just protected cells 
spp.pro <- colSums(sub.pro)   # the number of cells in which each species occur
occ.pro <- spp.pro 
occ.pro[] <- 0 
occ.pro[which(spp.pro>=occ) ] <- 1  # convert to presence/absence
sr.network <- length(which(occ.pro > 0)) # species richness in protected cells

# Null model for species richness 
sr.null <- replicate(runs, sr.rand.occ(dat, n.pro, occ)) # apply function sr.rand.occ multiple times to obtain a null distribution
ma.sr <- c(sr.network,sr.null) # add sr.network to null distribution
r1 <- rank(ma.sr)[1] # rank of sr.network in null distribution
ltpv <- r1/(runs+1) 
p <- 2 * min(ltpv, 1 - ltpv) # compute  two-tailed p-value
results.sr <- data.frame(threshold, n.pro, sr.network, mean.sr.rand=mean(sr.null), min.rand=min(sr.null), max.rand=max(sr.null), L95=quantile(sr.null,0.025), U95=quantile(sr.null,0.975),rank=r1, p, runs)
out.null.test <- rbind(out.null.test, results.sr) 

# Obtain outputs with representation of each species at that conservation target 
b <- data.frame(occ.pro); names(b) <- i 
out.occurrences <-  cbind(out.occurrences, b) 
} # end loop 

### Export results as tables 
write.table(out.null.test, "out.null.test.txt", row.names=F) # results of randomization tests
write.table(out.occurrences, "out.occurrences.txt") # for each taxon, presence/absence in reserves at each threshold under that conservation target

### Export results as a figure 
pdf("out.null.test.pdf") 
plot(out.null.test$threshold, out.null.test$sr.network*100/dim(dat)[2], pch=16, cex=1, xlab="Threshold", ylab="% of species", , ylim=c(0,100))
polygon(c(thresholds,rev(thresholds)),c(out.null.test$L95*100/dim(dat)[2],rev(out.null.test$U95*100/dim(dat)[2])),col = "grey80", border = FALSE)
lines(thresholds, out.null.test$mean.sr.rand*100/dim(dat)[2], type = 'l',  lty = 2)
points(out.null.test$threshold, out.null.test$sr.network*100/dim(dat)[2], pch=16, cex=1, xlab="Threshold", ylab="% of species")
dev.off() 
