#############################################
# Plotting Recombination Break Points
#############################################

### Read in the break points.txt file produced by the
### SNP binning script
brkpt.intervals <- read.table("P85_genotypes_wBrkPts2.txt", header=TRUE, stringsAsFactors=FALSE)

### Generate a column of transformed positions
### (Only need this if you want to plot chromosomes side by side in the same plot)
brkpt.intervals$TRANS <- brkpt.intervals$POS
chromosomes <- as.numeric(levels(as.factor(brkpt.intervals$CHR)))
to.add <- c(0)
for (c in 1:max(chromosomes)) {
    brkpt.intervals$TRANS[which(brkpt.intervals$CHR==c)] <- brkpt.intervals$TRANS[which(brkpt.intervals$CHR==c)] + to.add[c]
    chr.length <- max(brkpt.intervals$TRANS[which(brkpt.intervals$CHR==c)])
    to.add <- c(to.add, chr.length)
}

### Create a factor with 3 levels (1 for each genotype block to plot: A,B,H)
### Create a list of colors to use for plotting
temp <- factor(c("A","B","H"))
colors <- c("red","blue","yellow")

### If plotting all chromosomes together, create an empty plot NOW
### If plotting chromosome separately, create this plot inside the next loop
plot(c(0, max(brkpt.intervals$TRANS)), c(0, ((ncol(brkpt.intervals) - 3) * 1.5)), type= "n", xlab="", ylab="")

### Start looping through the chromosomes (even if plotting together, need to find end points sep.)
for (c in 1:max(chromosomes)) {
    chr.intervals <- subset(brkpt.intervals, brkpt.intervals$CHR==c)

    ## Next, plot the intervals for each individual separately
    for (ind in 3:(ncol(chr.intervals)-1)) {

        # Plotting coordinates will be stored in a "list of lists":
        # Each list corresponds to "A" "B" or "H",
        # and each list contains the all start and end points of each block
        all.coord <- list(A=list(c(0),c(0)), B=list(c(0),c(0)), H=list(c(0), c(0)))

        #**# If plotting chromosomes together, coordinates need to come from TRANS #**#        
        #**# If plotting chromosome individuall, coordinates need to come from POS #**#

        # Get the first (non-missing) position, determine which genotype it is,
        # and store it in the appropriate "start" list (use the factor levels to figure out the right list)
        first <- head((grep("[ABH]", chr.intervals[,ind], perl=TRUE)), n=1)
        first.pos <- chr.intervals$TRANS[first]
        gen <- substr((chr.intervals[first, ind]),1 ,1)
        t <- unlist(all.coord[[as.numeric(temp[grep(gen,temp)])]][1])
        t <- c(t,first.pos)
        all.coord[[as.numeric(temp[grep(gen,temp)])]][[1]] <- t

        # Get all of the break points (look for "/"), then find the genotypes on either side of "/"
        # Use the left gen. to assign position to "end" list, and right gen. to assign to "start" list
        bp.positions <- chr.intervals$TRANS[grep("/", chr.intervals[,ind])]
        if (length(bp.positions) > 0) {
            for (p in 1:length(bp.positions)) {
                code <- chr.intervals[(which(chr.intervals$TRANS==bp.positions[p])),ind]
                gens <- unlist(strsplit(code, "/"))
                end <- unlist(all.coord[[as.numeric(temp[grep(gens[1],temp)])]][2])
                end <- c(end, bp.positions[p])
                all.coord[[as.numeric(temp[grep(gens[1],temp)])]][[2]] <- end
                start <- unlist(all.coord[[as.numeric(temp[grep(gens[2],temp)])]][1])
                start <- c(start, bp.positions[p])
                all.coord[[as.numeric(temp[grep(gens[2],temp)])]][[1]] <- start
            }
        }
        # Put the last (non NaN) position into the appropriate "end" list
        # This will have to correspond to the last "start" gen used in the previous loop
        last <- chr.intervals$TRANS[which(chr.intervals$TRANS==max(chr.intervals$TRANS))]
        t2 <- unlist(all.coord[[as.numeric(temp[grep(gens[2],temp)])]][2])
        t2 <- c(t2, last)
        if (length(bp.positions) > 0) {
            all.coord[[as.numeric(temp[grep(gens[2],temp)])]][[2]] <- t2
        } else {
            all.coord[[as.numeric(temp[grep(gen,temp)])]][[2]] <- t2
        }

        # Get the y-coordinates (same for the entire length of individual)
        ytop <- (ind - 2) * 1.5
        ybottom <- ytop - 1.5

        # Finally, go through the list of coordinates and plot all rectangle for each
        # genotype block
        for (i in 1:3) {
            xleft <- all.coord[[i]][[1]]
            xright <- all.coord[[i]][[2]]
            for (j in 1:length(xleft)) {
                rect(xleft[j], ybottom, xright[j], ytop, col=colors[i], border=NA)
            }
        }
    }

    # Add a vertical line between chromosomes (if plotting them all together)
    abline(v=last)
}
