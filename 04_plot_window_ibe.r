# no scientific notation
options(scipen=999)

# minimum sites per window to include
n_sites <- 50000
# minimum variable sites per window to include
n_v_sites <- 50

# read in geographic and environmental distances
distances <- read.table("sitta_distances_new.txt", header=T, stringsAsFactors=F)

# read in unsorted fst values for all windows
fst_values <- read.table("window_fst.txt", header=T, stringsAsFactors=F)

# scaffold order
scaffolds <- paste0("Scaffold_", c(seq(from=1, to=10), seq(from=18, to=33), c(34, 12, 35, 36, 11, 13, 17, 37, 16, 14, 38, 39, 15)))

# define parts of output object
chrom <- list()
start <- list()
end <- list()
ibe_r2 <- list()
ibe_p <- list()
ibd_r2 <- list()
ibd_p <- list()
window_pattern <- list()
counter <- 1
# loop for each scaffold
for(a in 1:length(scaffolds)) {
	print(a)
	a_rep <- fst_values[fst_values$chr == scaffolds[a], ]
	# loop for each window
	window_starts <- sort(unique(a_rep$start))
	for(b in 1:length(window_starts)) {
		b_rep <- a_rep[a_rep$start == window_starts[b], ]
		# remove any comparisons without the number of sites needed
		b_rep <- b_rep[b_rep$number_sites >= n_sites & b_rep$number_variable_sites >= n_v_sites, ]
		# continue only if all pairwise comparisons have been kept
		if(nrow(b_rep) == nrow(distances)) {
			# add window to output
			chrom[[counter]] <- b_rep$chr[1]
			start[[counter]] <- b_rep$start[1]
			end[[counter]] <- b_rep$end[1]
			
			# calculate IBE and IBD for this window
			dist_rep <- c()
			# reorder Fst values to match distance matrix
			for(d in 1:nrow(distances)) {
				dist_rep <- c(dist_rep, b_rep$calculated_stat[(distances$Pop1[d] == b_rep$pop1 & distances$Pop2[d] == b_rep$pop2) | (distances$Pop2[d] == b_rep$pop1 & distances$Pop1[d] == b_rep$pop2)])
			}
			
			# IBD?
			regression_summary <- summary(lm(dist_rep ~ distances$g_dist))
			ibd_r2[[counter]] <- regression_summary$r.squared
			ibd_p[[counter]] <- regression_summary$coefficients[2,4]
			# IBE?
			regression_summary <- summary(lm(dist_rep ~ distances$new_e_distances))
			ibe_r2[[counter]] <- regression_summary$r.squared
			ibe_p[[counter]] <- regression_summary$coefficients[2,4]
			
			# add window pattern: IBD, IBE, both, or neither
			if(ibd_p[[counter]] <= 0.05 & ibe_p[[counter]] <= 0.05) {
				window_pattern[[counter]] <- "Both"
			} else if(ibd_p[[counter]] <= 0.05) {
				window_pattern[[counter]] <- "IBD"
			} else if(ibe_p[[counter]] <= 0.05) {
				window_pattern[[counter]] <- "IBE"
			} else {
				window_pattern[[counter]] <- "Neither"
			}
			
			# add to counter
			counter <- counter + 1
		}
	}
}

output <- data.frame(chrom=as.character(unlist(chrom)), start=as.numeric(unlist(start)), end=as.numeric(unlist(end)), ibe_r2=as.numeric(unlist(ibe_r2)), ibe_p=as.numeric(unlist(ibe_p)), ibd_r2=as.numeric(unlist(ibd_r2)), ibd_p=as.numeric(unlist(ibd_p)), window_pattern=as.character(unlist(window_pattern)))

ibe_p_corrected <- p.adjust(output$ibe_p, method="BH")
ibd_p_corrected <- p.adjust(output$ibd_p, method="BH")

output <- cbind(output, ibe_p_corrected, ibd_p_corrected)

# find what the window pattern is after p value correction
window_pattern_corrected <- c()
for(a in 1:nrow(output)) {
	if(output$ibd_p_corrected[a] <= 0.05 & output$ibe_p_corrected[a] <= 0.05) {
		window_pattern_corrected <- c(window_pattern_corrected, "Both")
	} else if(output$ibd_p_corrected[a] <= 0.05) {
		window_pattern_corrected <- c(window_pattern_corrected, "IBD")
	} else if(output$ibe_p_corrected[a] <= 0.05) {
		window_pattern_corrected <- c(window_pattern_corrected, "IBE")
	} else {
		window_pattern_corrected <- c(window_pattern_corrected, "Neither")
	}
}
output <- cbind(output, window_pattern_corrected)

# add a window number to the output
window <- seq(nrow(output))
output <- cbind(output, window)

################################################
################################################
################################################
################################################
################################################
# plot window IBE and IBD

# what are the unique chromosomes and their bounding areas for plotting?
total_windows <- nrow(output)
chr <- scaffolds
chr_polygons <- list()
# make the plotting polygons
for(a in 1:length(chr)) {
	a1 <- output$window[output$chrom == chr[a]]
	a2 <- a1[length(a1)]
	a1 <- a1[1]
	chr_polygons[[a]] <- rbind(c(a1, 0), c(a2, 0), c(a2, 4), c(a1, 4), c(a1, 0))
}

# set up plotting dimensions
par(mfrow=c(2,1))
par(mar=c(0.5,5,1,0))

# plot the corrected p-values (those with p <= 0.05 in black)
# IBE
plot(c(-1,-1), ylim=c(0,4), xlim=c(1, total_windows), xaxt="n", col="white", bty="n", cex.axis=1.1, cex.lab=1.3, ylab="-log10(p)", xlab="")
# plot polygons of scaffolds
odd <- 0
for(a in 1:length(chr_polygons)) {
	if(odd == 1) {
		polygon(chr_polygons[[a]], col="snow2", border="white")
		odd <- 0	
	} else {
		odd <- 1
	}
}
# plot points
points(window, -log10(output$ibe_p_corrected), pch=19, cex=0.1, col="darkgray")
points(window[output$ibe_p_corrected <= 0.05], -log10(output$ibe_p_corrected[output$ibe_p_corrected <= 0.05]), pch=19, cex=0.2, col="black")
title("IBE", adj=0.01, line=-1, font.main=1)

# IBD
plot(c(-1,-1), ylim=c(0,4), xlim=c(1, total_windows), xaxt="n", col="white", bty="n", cex.axis=1.1, cex.lab=1.3, ylab="-log10(p)", xlab="")
# plot polygons of scaffolds
odd <- 0
for(a in 1:length(chr_polygons)) {
	if(odd == 1) {
		polygon(chr_polygons[[a]], col="snow2", border="white")
		odd <- 0	
	} else {
		odd <- 1
	}
}
points(window, -log10(output$ibd_p_corrected), pch=19, cex=0.1, col="darkgray")
points(window[output$ibd_p_corrected <= 0.05], -log10(output$ibd_p_corrected[output$ibd_p_corrected <= 0.05]), pch=19, cex=0.2, col="black")
title("IBD", adj=0.01, line=-1, font.main=1)

# end plotting windows
################################################
################################################
################################################
################################################
################################################

# test the full genome relationship
fst <- c()
# reorder Fst values to match distance matrix
for(d in 1:nrow(distances)) {
	d_rep <- fst_values[(distances$Pop1[d] == fst_values$pop1 & distances$Pop2[d] == fst_values$pop2) | (distances$Pop2[d] == fst_values$pop1 & distances$Pop1[d] == fst_values$pop2), ]
	d_rep <- sum(d_rep$number_variable_sites * d_rep$calculated_stat) / sum(d_rep$number_variable_sites)
	fst <- c(fst, d_rep)
}
distances <- cbind(distances, fst)

par(mfrow=c(1, 3))
plot(distances$g_dist, distances$e_dist, xlab="Geographic distance (km)", ylab="Environmental distance", pch=19, cex=0.8, xlim=c(0,500), ylim=c(0,3))
abline(lm(distances$e_dist ~ distances$g_dist))
plot(distances$g_dist, distances$fst, xlab="Geographic distance (km)", ylab="Fst", pch=19, cex=0.8, ylim=c(0,0.2))
abline(lm(distances$fst ~ distances$g_dist))
plot(distances$e_dist, distances$fst, xlab="Environmental distance", ylab="Fst", pch=19, cex=0.8, ylim=c(0,0.2))
abline(lm(distances$fst ~ distances$e_dist))












