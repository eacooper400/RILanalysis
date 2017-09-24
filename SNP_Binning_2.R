setwd("/Users/boylesre/Documents/GeneticMap")

install.packages("zoo")
library(zoo)

#########################################################
#*#*#*#   PREDEFINED FUNCTIONS  #*#*#*#*#*#*#*#*#*#*#*#*#
#########################################################
calcProp <- function(x) {
    A <- length(grep("A", x))
    B <- length(grep("B", x))
    prop <- A/B
    return(prop)
}
defineTrans <- function(x) {
    # Define the possible strings that indicate a break point
    trans.strings <- c("40,100,250","-250,-100,-40","140,-100,-40","40,100,-100,-40","-250,-100,100,250","-250,-140","140,250","40,100,-100,250","40,100,-140","-250,-100,100,-40","140,-100,100,250","40,100,-100,-140","140,-100,100,-140","140,-100,100,-40","-250,-100,100,-140","140,-140")
    bp.start.type <- c("A","B","A/H","A/H","B/H","B","A","A/H","A/H","B/H","A/H","A/H","A/H","A/H","B/H","A/H")
    bp.end.type <- c("A/B","B/A","H/A","H/A","H/B","B/A","A/B","H/B","H/A","H/A","H/B","H/A","H/A","H/A","H/A","H/A")

    # For each string in seqs, look up to see if it corresponds to a break point type
    s.type <- c()
    e.type <- c()

    for (i in 1:length(x)) {
        type <- grep(x[i], trans.strings)
        if (length(type) == 0) {
            s.type <- c(s.type, 0)
            e.type <- c(e.type, 0)
        } else {
            s.type <- c(s.type, bp.start.type[type])
            e.type <- c(e.type, bp.end.type[type])
        }
    }
    return(list(s.type,e.type))
}
fixEnds <- function(x,y) {
    starts <- x
    ends <- y
    new.ends <- c()
    for (s in 1:length(starts)) {
        dist <- c()
        for (e in 1:length(ends)) {
            d <- ends[e] - starts[s]
            dist <- c(dist,d)
        }
        best <- which(dist==min(dist[which(dist>0)]))
        new.ends <- c(new.ends, ends[best])
    }
    return(new.ends)
}
fixStarts <- function(x,y) {
    starts <- x
    ends <- y
    new.starts <- c()
    for (e in 1:length(ends)) {
        dist <- c()
        for (s in 1:length(starts)) {
            d <- ends[e] - starts[s]
            dist <- c(dist,d)
        }
        best <- which(dist==min(dist[which(dist>0)]))
        new.starts <- c(new.starts, starts[best])
    }
    return(new.starts)
}

#####################################################################
#####################################################################
#hmp <- read.table(file="BTX_ABH_for_Bins.hmp.txt", header=TRUE, stringsAsFactors=FALSE)
hmp <- read.table(file="P85_ABH_for_Bins.hmp.txt", header=TRUE, stringsAsFactors=FALSE)

#Adding columns for chrom and pos to match function
#hmp <- hmp[,c(3,4,281,2,5:280)]
hmp$pos <- hmp$id
hmp <- hmp[,c(1,281,2:280)]
hmp$chrom <- hmp$id
#hmp <- hmp[,c(1,282,2:281)]
hmp <- hmp[,c(282,2:281)]
hmp$chrom <- gsub("\\_.*", "", hmp$chrom)
hmp$chrom <- gsub("S", "", hmp$chrom)
hmp$pos <- gsub("^.*?_","",hmp$pos)

genoWins <- data.frame()

#Changed max(chromosomes) to 10 since it was excluding chromosome 10
for (c in 1:10) {
    chr.hmp <- subset(hmp, hmp$chrom==c)
    
    gCalls <- rollapply(chr.hmp[,-c(1:3)], width=15, by=1, FUN=calcProp, partial=TRUE)
    chrWins <- cbind(chr.hmp$chrom, chr.hmp$pos)
    colnames(chrWins) <- c("CHR", "POS")
    chrWins <- cbind (chrWins, gCalls)
    genoWins <- rbind(genoWins, chrWins)
}

genoWins.mat <- as.matrix(genoWins)
genoWins.mat[genoWins.mat=="Inf"] <- "15"
genoWins <- data.frame(genoWins.mat)

for(i in 1:ncol(genoWins)) {
    genoWins[,i] <- as.numeric(as.character(genoWins[,i]))
}

#### Altering the next part slightly, so now there is an H1 and H2, to keep track of when we cross the 8:7/7:8 threshold in the hets
#### Also, I'm going to encode A,B,H1,and H2 numerically so that I can manipulate them more easily to find transition points later

genoWins.recode <- genoWins
genoWins.recode[genoWins>2] <- 10
genoWins.recode[genoWins<0.25] <- 400

genoWins.recode[(genoWins<=2) & (genoWins>=(8/7))] <- 50
genoWins.recode[(genoWins<(8/7)) & (genoWins>=0.25)] <- 150

# After recoding, replace chromosome and position columns w/ original values
genoWins.recode[,1:2] <- genoWins[,1:2]

### Create a new data frame to hold break point values for all chromosomes
brkpt.intervals <- data.frame()

#### Find the break point intervals separately for each chromosome, then add them to the full brkpt.intervals table afterwards
for (c in 1:10) {
    genoWins.chr <- subset(genoWins.recode, genoWins.recode$CHR==c)

    # This first step creates a matrix of pairwise comparisons between rows (positions)
    # If the 2nd position has the same genotype code as the first, the comparison will equal zero
    # If it is non-zero, some kind of transition has occurred (but it is still unknown if it is a break point)
    tmp.matrix <- matrix(0, nrow=nrow(genoWins.chr), ncol=ncol(genoWins.chr))

    for (i in 2:nrow(genoWins.chr)) {
        trans.code <- unname(unlist(genoWins.chr[i,3:ncol(genoWins.chr)])) - unname(unlist(genoWins.chr[i-1,3:ncol(genoWins.chr)]))
        trans.code <- c(genoWins.chr[i,1], genoWins.chr[i,2], trans.code)
        tmp.matrix[i,] <- trans.code
    }
    tmp.matrix[1,] <- unname(unlist(genoWins.chr[1,]))

    #### Next is the most complicated part, where we need to figure out where transitions from homozygous/homozygous are vs. transitions that are hom/het
    #### This will depend on if the "line" crosses over the 8:7 threshold more than once, hence the earlier encoding

    # First, create a new matrix to hold break point codes for just the current chromosome
    test.bp <- matrix(0, nrow=nrow(genoWins.chr), ncol=ncol(genoWins.chr))
    test.bp[,1] <- unname(unlist(genoWins.chr[,1]))
    test.bp[,2] <- unname(unlist(genoWins.chr[,2]))

    # Then, loop through each column (individual) of tmp.matrix to find start and end points of potential recombination intervals
    for (i in 3:ncol(tmp.matrix)) {
        tmp <- tmp.matrix[,i]

        # First, I need to get the "starts" of all possible transitions; I'm defining these are points where an individual starts as 1 parental haplotype,
        # and then changes to any other state (the way the encoding works is too involved to type out here; I will have to write up something separate to explain it)
        starts <- which(tmp == 40 | tmp == 140 | tmp == -350 | tmp == -250)

        # Next, I'll get the "ends" of all transitions: these are points where an individual transitions INTO a parental haplotype from something else
        ends <- which(tmp == -40 | tmp == 350 | tmp == -140 | tmp == 250)

        # Check that there is at least one of each start and/or end; otherwise, no break points
        # Check that there are the same number of start and end points; if not, fix whichever vector has too many values
        if ((length(starts) == 0) | (length(ends) == 0)) {
            starts <- c()
            ends <- c()
        }
        if (length(starts) > length(ends)) {
            starts <- fixStarts(starts,ends)
        }
        else if (length(ends) > length(starts)) {
            ends <- fixEnds(starts,ends)
        }

        # Codes 390 and -390 are special codes where 1 hap transitions directly into the other hap (so they are unambiguous hom/hom break points)
        directBP1 <- which(tmp == 390)
        directBP2 <- which(tmp == -390)

        # This next loop will find all transitions between every start and end point, then determine if the transition is from a
        # hom/hom break point or a hom/het break point (or is not a break point at all)
        intervals <- genoWins.chr[,i]
        if (length(starts) > 0) {
            seqs <- c()
            for (j in 1:length(starts)) {
                all.values <- tmp[starts[j]:ends[j]]
                no.zero <- all.values[all.values!=0]
                if ((length(no.zero) > 3)) {
                    x <- paste(c(no.zero[1:3], tail(no.zero, n=1)), collapse=",")
                    seqs <- c(seqs, x)
                } else {
                    x <- paste(no.zero, collapse=",")
                    seqs <- c(seqs, x)
                }
            }
            new.codes <- defineTrans(seqs)
            intervals[starts] <- unlist(new.codes[1])
            intervals[ends] <- unlist(new.codes[2])
        }
        # Now, replace the values from genoWins.recode with a break point term (if they correspond to one)
        # or the genotype symbol "A", "B", "H1" or "H2"
        intervals[directBP1] <- "A/B"
        intervals[directBP2] <- "B/A"

        intervals[grep("^10$",intervals)] <- "A"
        intervals[grep("^150$",intervals)] <- "H2"
        intervals[grep("^400$",intervals)] <- "B"
        intervals[grep("^50$",intervals)] <- "H1"

        # Replace the column in the current chromosome break point matrix with the new values
        test.bp[,i] <- intervals
    }

    # Once all individuals have been checked for break points, add the chromosome matrix to the main data frame
    brkpt.intervals <- rbind(brkpt.intervals, test.bp)
}
colnames(brkpt.intervals) = colnames(genoWins)
rownames(brkpt.intervals) = rownames(genoWins)

for(i in 1:ncol(brkpt.intervals)) {
    if (i <= 2) {
        brkpt.intervals[,i] <- as.numeric(as.character(brkpt.intervals[,i]))
    } else {
        brkpt.intervals[,i] <- as.character(brkpt.intervals[,i])
    }
}

# Count the break points (without counting the same position twice, so go by row)
num.bp <- 0
for (i in 1:nrow(brkpt.intervals)) {
    check <- grep("/", brkpt.intervals[i,])
    if (length(check) != 0) {
        num.bp <- num.bp + 1
    }
}

num.het <- 0
for (i in 1:nrow(brkpt.intervals)) {
    check <- grep("H/", brkpt.intervals[i,])
    if (length(check) != 0) {
        num.het <- num.het + 1
    }
}
num.het <- num.het*2

# Write this file to a table, since it also takes some time to generate
write.table(brkpt.intervals, "P85_genotypes_wBrkPts.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)


