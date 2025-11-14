#################################################################################################
## Assessing species’ representation in protected area networks by calculating mean            ##
## percentage overlap between species distribution and protected areas (grid cells)            ##        
##                                                                                             ##
## Pedro Abellán                                                                               ##
## April 2014                                                                                  ##
## Citation: Sánchez-Fernández D, Abellán P. 2015. Using null models to identify               ##
## under-represented species in protected areas: A case study using European amphibians        ##
## and reptiles. Biological Conservation 184: 290–299.                                         ##
#################################################################################################

### Load data and set parameters
dat <- read.table("matrix.txt", header=TRUE, row.names=1)  # load presence/absence matrix; cells in rows and species in columns;
# it must contain row and colums names
cells <-  read.table("cells.txt", header=TRUE, row.names=1) # load table with information about grid cells;
# it must contain a field with the percentaje of overlap with protected areas ("protected" in this example)

runs <- 999 # set number of runs for randomizations
pro <- "protected" # name of the field in table "cells" with the percentaje of overlap with protected areas

### function to compute mean % protected area in k random cells
propor <- function(k, pro){
sam <- sample(c(1:nrow(cells)), k)
prop <- mean(cells[sam,pro])
return(prop)
}

### compute mean % protected area for all species in dataset
out <- data.frame() # prepare empty data frame to store results
for (j in 1:ncol(dat)){ # start loop
cells.sp <- row.names(dat[which(dat[,j]==1),]) # identify cells in which the jth species occur
n.cells <- length(cells.sp) # number of cells in which the jth species occur
mean.pro <- mean(cells[which(row.names(cells)%in%cells.sp), pro]) # mean percentage overlap for jth species
result <- data.frame(species=names(dat)[j], n.cells, mean.pro) # data frame with results for jth species
out <- rbind(out, result) # combine results among loops
} # end loop

### Null distribution of expected % protected area in k random cells given the dataset
min.cells <- min(out$n.cells) # minimum number of cells in which species occur in the dataset
max.cells <- max(out$n.cells) # maximum number of cells in which species occur in the dataset
m <- matrix(ncol=length(min.cells:max.cells),nrow=runs)  # create a empty matrix to store results
for (k in min.cells:max.cells){ # start loop
m[,k-(min.cells-1)] <- replicate(runs, propor(k,pro)) # apply function propor multiple times for k cells
} # end loop


### Export resuls as a table

# add new colunms to store results
out$mean.random <- NA
out$CI.lower <- NA
out$CI.upper <- NA
out$rank <- NA
out$p.value <- NA

for (i in 1:nrow(out)){ # start loop
null.i <- m[,out[i,"n.cells"]-(min.cells-1)]  # select relevant column in matrix m
ma.i<-c(out$mean.pro[i],null.i) # combine mean percentage overlap and null distribution in k random cells
r1<-rank(ma.i)[1] # rank of observed mean % protected area in null distribution of equivalent number of cells
ltpv <- r1/(runs+1)
p.value <- 2 * min(ltpv, 1 - ltpv) # compute  two-tailed p-value
out[i,4] <- mean(null.i)
out[i,5] <- quantile(null.i,0.025)
out[i,6] <- quantile(null.i,0.975)
out[i,7] <- r1
out[i,8] <- p.value
} # end loop

write.table(out, "out_mpo.txt", row.names=F)

### Export results as a figure

L <- apply(m, 2, quantile, 0.025)
U <- apply(m, 2, quantile, 0.975)
M <-  apply(m, 2, mean)
df <- data.frame(x=min.cells:max.cells, M, L, U)

pdf("out_mpo.pdf")
plot(out$n.cells, out$mean.pro, pch=16, cex=1, xlab="Number of cells", ylab="% of protected area")
polygon(c(df$x,rev(df$x)),c(df$L,rev(df$U)),col = "grey70", border = FALSE)
lines(min.cells:max.cells, M, type = 'l',  lty = 5)
points(out$n.cells, out$mean.pro, pch=16, cex=1)
text(x=(out$n.cells), y=out$mean.pro+((max(out$mean.pro)-(out$mean.pro))*0.02), labels=out$sp, cex=0.7) # (optional) show species labels
dev.off()
