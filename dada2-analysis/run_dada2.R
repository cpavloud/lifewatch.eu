#!/usr/bin/env Rscript
library(argparser)
library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(vegan); packageVersion("vegan")
library(tools)
library(jsonlite)

sessionInfo()

	#####################################################################
	# Graphing function for NMDS plot
	#####################################################################

nice.plot.scatter = function(title, subtitle=NULL, xLabel=NULL, yLabel=NULL, data.as.matrix) {
    series = list()
    for (row in 1:nrow(data.as.matrix)) {
        data = list(c(data.as.matrix[row, 1], data.as.matrix[row, 2]))
        name = row.names(data.as.matrix)[row]
        thing = list()
        thing$name = name
        thing$data = data
        series[[row]] = thing
    }
    plot = list()
    plot$title$text = unbox(title)
    if (!is.null(subtitle)) {
        plot$subtitle$text = unbox(subtitle)
    }
    if (!is.null(xLabel)) {
        plot$xAxis$title$text = unbox(xLabel)
    }
    if (!is.null(yLabel)) {
        plot$yAxis$title$text = unbox(yLabel)
    }
    plot$chart$type = unbox("scatter")
    plot$chart$zoomType = unbox("xy")
    plot$series = series
    toJSON(plot)
}

	#####################################################################
	# Options
	#####################################################################

parser <- arg_parser("Arguments")

parser <- add_argument(parser, "--truncate_length_fwd", short="-tf", type="integer", default=240, 
		help="Truncate length of forward reads [default %default]")

parser <- add_argument(parser, "--truncate_length_rev", short="-tr", type="integer", default=160, 
		help="Truncate length of reverse reads [default %default]")

parser <- add_argument(parser, "--truncate_quality", short="-tq", type="integer", default=11, 
		help="Truncates based on quality score [default %default]")

parser <- add_argument(parser, "--fwd_reads", short="-f", type="character", 
		help="dataset file name")

parser <- add_argument(parser, "--rev_reads", short="-r", type="character", 
		help="dataset file name")

parser <- add_argument(parser, "--output", short="-o", type="character", 
		help="output file name")

parser <- add_argument(parser, "--path_to_silva", short="-s", type="character", 
		help="Path to silva database")

parser <- add_argument(parser, "--max_n", short="-n", type="integer", default=0, 
		help="maximum of ambiguous bases (Ns) [default %default]")

parser <- add_argument(parser, "--rank", short="-r", type="character", default="Family", 
		help="Taxonomic rank to display in abundance bar plot (note that 'Genus' may display oddly)")

opt <- parse_args(parser)

    #####################################################################
    # clean and match
    #####################################################################

	## Removes all sequences containing any ambiguous bases (N).
	## Tries to match up all pairs of forward and reverse reads by ID,
	## Removes all sequences that have been orphaned by filtering.
	## Creates new directory for clean and matched files

	## set up for interactive testing of script
	#path_to_input <- "/Users/alper/repos/test"
	#opt$output <- "/Users/alper/repos/test/output"
#
	#opt$fwd_reads <- list.files(path_to_input, pattern = "R1_001.fastq.gz", all.files = TRUE, full.names = TRUE, recursive = FALSE)
	#opt$rev_reads <- list.files(path_to_input, pattern = "R2_001.fastq.gz", all.files = TRUE, full.names = TRUE, recursive = FALSE)
	#opt$path_to_silva <- "/Users/alper/repos/test/silva/silva_nr_v132_train_set.fa"
	#path_to_fwd <- opt$fwd_reads
	#path_to_rev <- opt$rev_reads

	#####################################################################
	# Setting file paths
	#####################################################################

	path_to_fwd <- unlist(strsplit(opt$fwd_reads, ','))   # skip when testing interactively
	path_to_rev <- unlist(strsplit(opt$rev_reads, ','))   # skip when testing interactively
	path_to_output <- opt$output

	fwd_table <- read.table(text = basename(path_to_fwd), sep = ".", as.is = TRUE)
	rev_table <- read.table(text = basename(path_to_rev), sep = ".", as.is = TRUE)

	filtered_path <- file.path(path_to_output, "filtered")
	denoised_path <- file.path(path_to_output, "denoised")

	compressed_input = TRUE
	for(i in path_to_fwd) {
		if (file_ext(i) != "gz") {
			compressed_input <- FALSE
			}
		}

	if (compressed_input == TRUE) {
		filtered_fwd <- file.path(filtered_path,paste(fwd_table[,1],"filtered",fwd_table[,2], fwd_table[,3], sep="."))
		filtered_rev <- file.path(filtered_path,paste(rev_table[,1],"filtered",rev_table[,2], rev_table[,3], sep="."))
	} else {
		filtered_fwd <- file.path(filtered_path, paste(fwd_table[,1],"filtered",fwd_table[,2], "gz", sep="."))
		filtered_rev <- file.path(filtered_path, paste(rev_table[,1],"filtered",rev_table[,2], "gz", sep="."))
	}

	denoised <- file.path(denoised_path, paste(fwd_table[,1],"denoised.fasta", sep="."))

	#####################################################################
	# Checking the raw input FASTQ files
	#####################################################################

	write(("\n Generating raw read quality profiles..."), stdout())

	png(file.path(path_to_output, "RawQualityProfile_fwd.png"))
	plotQualityProfile(file.path(path_to_fwd[1:length(fwd_table)]))
	garbage <- dev.off()
	png(file.path(path_to_output, "RawQualityProfile_rev.png"))
	plotQualityProfile(file.path(path_to_rev[1:length(rev_table)]))
	garbage <- dev.off()
	write(("Done. \n"), stdout())

	if(length(fwd_table) != length(rev_table))
        # DevNote: Returns error to stderr
        write(paste0("The number of input files looks dodgy.\n",
              "number of forward reads: ", length(path_to_fwd),
              "\number of reverse reads: ", length(path_to_rev)),
              stderr())

	write(file.path(path_to_fwd), stdout())
	write(file.path(filtered_fwd), stdout())
	write(file.path(path_to_rev), stdout())
	write(file.path(filtered_rev), stdout())

	write(file.path(denoised), stdout())

	#####################################################################
	# Filtering and trimming of reads
	#####################################################################

	out <- filterAndTrim(fwd = file.path(path_to_fwd),
			filt = file.path(filtered_fwd),
			rev = file.path(path_to_rev),
			filt.rev = file.path(filtered_rev),
			truncLen = c(opt$truncate_length_fwd, opt$truncate_length_rev),
			maxEE = 2,
			truncQ = opt$truncate_quality,
			maxN = opt$max_n,
			rm.phix=TRUE,
			compress=TRUE,
			verbose=TRUE,
			multithread=TRUE)

	write(("\n Generating filtered read quality profiles..."), stdout())
	png(file.path(path_to_output, "FilteredQualityProfile_fwd.png"))
	plotQualityProfile(file.path(filtered_fwd))
	garbage <- dev.off()
	png(file.path(path_to_output, "FilteredQualityProfile_rev.png"))
	plotQualityProfile(file.path(filtered_rev))
	garbage <- dev.off()
	write(("Done. \n"), stdout())

	if (length(path_to_fwd)*2 != length(list.files(filtered_path))) {
		write(paste0("Filtering did not result in any sequences from some files.\n",
			"The number of filtered sequence files is ",
			length(list.files(filtered_path)),
			" \n",
			"while the number of raw files is ",
			length(path_to_fwd) + length(path_to_rev),
			" \n"),
			stdout())
	}

	#####################################################################
	# Prepare variables for further sequence preparation
	#####################################################################

	sample.names <- sapply(strsplit(basename(filtered_fwd), ".fastq"), `[`, 1)	# Assumes filename = samplename_XXX.fastq.gz
	sample.names.rev <- sapply(strsplit(basename(filtered_rev), ".fastq"), `[`, 1)	# Assumes filename = samplename_XXX.fastq.gz
	names(filtered_fwd) <- sample.names
	names(filtered_rev) <- sample.names

	#####################################################################
	# Learn error rates
	#####################################################################

	write(("\n Learning error rates..."), stdout())
	set.seed(100)
	error_fwd <- learnErrors(filtered_fwd, nbases=1e8, multithread=TRUE)
	error_rev <- learnErrors(filtered_rev, nbases=1e8, multithread=TRUE)
	write(("Done. \n"), stdout())

	#####################################################################
	# Dereplication, denoising and merger of paired-end reads
	#####################################################################

	mergers <- vector("list", length(sample.names))
	names(mergers) <- sample.names
	dada_fwd <- vector("list", length(sample.names))
	dada_rev <- vector("list", length(sample.names))

	for(sam in sample.names) {
		cat("Processing:", sam, "\n")
			derepF <- derepFastq(filtered_fwd[[sam]])
			ddF <- dada(derepF, err = error_fwd, multithread=TRUE)
			dada_fwd[[sam]] <- ddF
			derepR <- derepFastq(filtered_rev[[sam]])
			ddR <- dada(derepR, err = error_rev, multithread=TRUE)
			dada_rev[[sam]] <- ddR
			merger <- mergePairs(ddF, derepF, ddR, derepR)
			mergers[[sam]] <- merger
	}

	rm(derepF); rm(derepR); rm(merger); rm(ddF); rm(ddR)

	#####################################################################
	# Sequence table construction and chimera removal
	#####################################################################

	seqtab <- makeSequenceTable(mergers)
	write(paste0("Merging resulted in ", ncol(seqtab), " sequences."), stdout())

	if(ncol(seqtab) == 0){
		# DevNote: Returns error to stderr
		write(paste0("Merging did not result in any sequences.\n", "You may consider a different truncation for forward and reverse reads. For now the program will continue running the concatenate option instead."),
		stdout())
		rm(seqtab)
		for(sam in sample.names) {
			cat("Processing:", sam, "\n")
			derepF <- derepFastq(filtered_fwd[[sam]])
			ddF <- dada(derepF, err = error_fwd, multithread=TRUE)
			dada_fwd[[sam]] <- ddF
			derepR <- derepFastq(filtered_rev[[sam]])
			ddR <- dada(derepR, err = error_rev, multithread=TRUE)
			dada_rev[[sam]] <- ddR
			concatenate <- mergePairs(ddF, derepF, ddR, derepR, justConcatenate=TRUE)
			mergers[[sam]] <- concatenate
		}
		seqtab <- makeSequenceTable(mergers)
		write(paste0("Concatenating the sequences resulted in ", ncol(seqtab), " sequences."))
		rm(derepF); rm(derepR); rm(concatenate); rm(ddF); rm(ddR)
		seqtab.final <- seqtab
	} else {
		seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)  ### should be used for plotting taxonomy
		saveRDS(seqtab.nochim, file.path(path_to_output, "seqtab_nochim.rds"))
		seqtab.final <- seqtab.nochim
	}

	saveRDS(seqtab.final, file.path(path_to_output, "seqtab.rds"))
	write(paste0("Overall abundance of chimeric sequences: ", 1-(sum(seqtab.nochim)/sum(seqtab)), "%"), stdout())
	write(("Done. \n"), stdout())

	#####################################################################
	# Print denoised, de-chimera-ed sequences by sample
	#####################################################################

	# DevNote - see if there is a more elegant way of achieving this
	# DevNote - denoised sample names currently displaying with fwd read filename, even though reads have been merged

	dir.create(denoised_path)

	for (i in rownames(seqtab.final)) {
		fasta_name <- paste(denoised_path, "/" ,gsub("filtered","denoised",i), ".fasta", sep="")
		seq_names <- vector()
		my_seqs <- vector()
		seqno <- 1
		for (j in 1:dim(seqtab.final)[2]) {
			if (seqtab.final[i,j] != 0) {
				seq_names[seqno] <- paste(">", seqno, "_", seqtab.final[i,j], "x", sep="")
				my_seqs[seqno] <- colnames(seqtab.final)[j]
				seqno <- seqno + 1
			}
		}
	my_fasta <- c(rbind(seq_names, my_seqs))
	write(my_fasta, fasta_name)
	}

	#####################################################################
	# Assign taxonomy
	#####################################################################

	write(("\n Starting taxonomy analysis."), stdout())

        tax <- assignTaxonomy(seqtab.final, opt$path_to_silva, multithread=TRUE, tryRC=TRUE)  ### should be used for plotting taxonomy

	# Write to file
        saveRDS(tax, file.path(path_to_output, "tax_final.rds"))
	# tax = readRDS(file.path(path_to_output, "tax_final.rds"))

	write(("Finalized taxonomy analysis. \n"), stdout())

	#####################################################################
	# Track sequences through analysis
	#####################################################################
	getN <- function(x) sum(getUniques(x))
	track <- cbind(out,
			  	   sapply(mergers, getN),
			  	   rowSums(seqtab.final))
	# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
	colnames(track) <- c("input", "filtered", "merged", "final")
	rownames(track) <- sample.names

	# write to file
	write.csv(track, file.path(path_to_output, "track_sequence_stats.csv"))
	write(("\n Write tracking file to output. \n"), stdout())

	#####################################################################
	# Run stats and visualization
	#####################################################################
	theme_set(theme_bw())

	# Parsed matrix for taxonomy plots
	ps <- phyloseq(otu_table(seqtab.final, taxa_are_rows=FALSE), tax_table(tax))

	otu_matrix <- as(otu_table(ps), "matrix")
	tax_matrix <- as(tax_table(ps), "matrix")
	otutax_matrix <- cbind(t(otu_matrix), tax_matrix)

	# Plot sequence abundance by taxon
	# DevNote - Needs refining, particularly all of the 'else if' arguments

	write(("\n Generating per-sample taxon abundance plot. \n"), stdout())

	png(file.path(path_to_output, "TaxonAbundance.png"))
	p <- plot_bar(ps, fill=opt$rank, title="Per-sample abundance of taxa")

	if (opt$rank == "Genus") {
		p + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
	} else if (opt$rank == "Family") {
		p + geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack")
	} else if (opt$rank == "Order") {
		p + geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack")
	} else if (opt$rank == "Class") {
		p + geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack")
	} else if (opt$rank == "Phylum") {
		p + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
	}

	garbage <- dev.off()

	# rarefyed_ps <- rarefy_even_depth(ps)

	# Rarefy multiple times
	seqtab_for_nmds <- seqtab.final[rowSums(seqtab.final)!=0,]
        # DevNote: Returns error to stderr
	if(nrow(seqtab_for_nmds) != nrow(seqtab.final))
        # DevNote: Returns error to stderr
		write(paste0("The analyses did not result in any sequences in ", length(which(rowSums(seqtab.final)==0)), " sample(s). \nThe sample(s) without sequences is/are following: ", toString(names(which(rowSums(seqtab.final)==0)))),
		stdout())

	rarefied <- rrarefy(seqtab_for_nmds, min(rowSums(seqtab_for_nmds)))
	rarefied_set <- rarefied[,colSums(rarefied)!=0] # final 16s set:  samples and  OTUs

	# Plots ordination (only with more than 3 samples)
	if(length(fwd_table) < 3){
	# DevNote: Returns error to stderr
		write(paste0("The number of samples is ", length(fwd_table), ".\n", "You need to provide at least 4 samples!"),
		stdout())
	} else {
		# plot_ordination(ps, ordinate(ps, "MDS")) + geom_point(size = 5)


		## provide NMDS plot and stress value - these should be plotted!!!
	   		# nmds (only with more than 3 samples)
		nmds <- metaMDS(rarefied_set, distance='bray', trymax=200)
#		save.image("fubar.rda") # Load this image with
		scores_nmds <- scores(nmds) ### coordinates for xy plot - by clicking on points one should see sample name
		# Plot 'scores_nmds' (data as a two column table). Include the Int 'stress_nmds' somewhere in teh plot.

		json = nice.plot.scatter(title="NMDS plot", subtitle="Something something", data.as.matrix=scores_nmds)
		write(json, file.path(path_to_output, "nmds.graph"))
		stress_nmds <- nmds$stress  ### stress value should be included somewhere in the plot
		write.csv(scores_nmds, file.path(path_to_output, "scores_nmds.csv"))
		write.csv(stress_nmds, file.path(path_to_output, "stress_nmds.csv"))
	}

  	## Diversity estimates
	R5 = as.data.frame(t(estimateR(rarefied_set)))	 # chao1 and ACE

	H5  =  diversity(rarefied_set)                   # shannon

	S  =  specnumber(rarefied_set)
	J5  =  H5/log(S)                         		 #Pielou's
	diversity_estimates = cbind(R5, J5, H5)

	colnames(diversity_estimates) = c("observed", "chao1", "se.chao1", "ACE", "se.ACE", "pielou_evenness", "shannon")

	write.csv(diversity_estimates, file.path(path_to_output, "diversity_estimates_stats.csv"))

	write(paste0("Analysis run to the end."), stdout())

	### ACE print with Std.dev. per sample - pielou and shannon per sample each in one graph - clicking on points should give sample name - names should also be given at bottom.

	### Other option might be printing the table diversity estimates

#	save.image(file.path(path_to_output, "rrun_final.RData"))
