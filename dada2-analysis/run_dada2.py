#!/usr/bin/env python

# Licence:

import sys
import argparse
import logging as log
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
dada2 = importr("dada2")
ggplot2 = importr("ggplot2")
phyloseq = importr("phyloseq")
plyr = importr("plyr")
reshape = importr("reshape")

parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        description="Run a dada2 analysis on fastq files."
)
parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Be more verbose"
)
parser.add_argument(
        "--truncate_length_fwd",
        type=int,
        default=240,
        help="Truncate forward reads"
)
parser.add_argument(
        "--truncate_length_rev",
        type=int,
        default=160,
        help="Truncate reverse reads"
)
parser.add_argument(
        "--truncate_quality",
        type=int,
        default=11,
        help="Truncates based on quality"
)
parser.add_argument(
    "--input-files",
    help=(
        "Comma separated list of input files. If this argument is omitted, "
        "all files with file ending '.fastq' will be used instead"
    )
)
parser.add_argument(
    "--path_to_silva",
    help="Alex, please add"
)
parser.add_argument(
        "--max_n",
        default=0,
        type=int,
        help="Alex, please add"
)
args = parser.parse_args()


def main():

    #####################################################################
    # Specific functions
    #####################################################################

    robjects.r("getN = function(x) sum(getUniques(x))")

    #####################################################################
    # default settings
    #####################################################################

    rcode = "path_to_input = '/media/input'"
    robjects.r(rcode)

    rcode = "path_to_output = '/media/output'"
    robjects.r(rcode)

    # Initialise log file
    logfile = '/media/output/edna_out.log'
    log.basicConfig(filename=logfile,
                    format='%(asctime)s %(message)s',
                    level=log.DEBUG)
    log.info("Analysis started.")

    # Store the filtering parameters as R variables
    rcode = "truncate_length_fwd = %s" % args.truncate_length_fwd
    robjects.r(rcode)
    rcode = "truncate_length_rev = %s" % args.truncate_length_rev
    robjects.r(rcode)
    rcode = "truncate_quality = %s" % args.truncate_quality
    robjects.r(rcode)
    rcode = "max_n = %s" % args.max_n
    robjects.r(rcode)

    # Store the path to the Silva data base as an R variable
    rcode = "path_to_silva = \'%s\'" % args.path_to_silva
    robjects.r(rcode)

    #####################################################################
    # loading and manipulation
    #####################################################################

    # Sort the input files to ensures the forward and reverse
    # reads are in same order.

    log.info("Sorting input files")

    if args.input_files:
        infiles = [f.strip() for f in args.input_files.split(",")]
        # TODO: There's probably a better way of doing this:
        file_list = ','.join(['"{}"'.format(f) for f in infiles])
        robjects.r('infiles = c({})'.format(file_list))
    else:
        robjects.r("infiles = list.files(path_to_input, pattern='.fastq')")

    # DevNote: Hard-coded file extention
    robjects.r('fastqFs = sort(infiles[grep("*_R1_001.fastq", infiles)])')
    robjects.r('fastqRs = sort(infiles[grep("*_R2_001.fastq", infiles)])')

    ###########################################################################
    # Sanity check. Is there an equal number of forward and reverse read-files?
    ###########################################################################

    if len(robjects.r("fastqFs")) != len(robjects.r("fastqRs")):
        # DevNote: Return something usefull to the frontend
        sys.exit("Bummer, the number of input files looks dodgy")

    # Extract sample names, assuming the filenames have
    # the format: SAMPLENAME_XXX.fastq
    # DevNote: Make sure file name is in the right format.
    robjects.r("sample.names <- sapply(strsplit(fastqFs, '_'), `[`, 1)")

    #####################################################################
    # quality profiles of forward and reverse reads
    #####################################################################

    log.info("Running quality check on input reads.")

    robjects.r(
            '''
            qual_profile_path <- file.path(path_to_output, "Quality_profiles")
            dir.create(file.path(qual_profile_path))

            for (fastq in c(fastqFs, fastqRs)){

            filename_safe = gsub("/", "_", fastq)
            file_name = paste(qual_profile_path, "/", filename_safe, ".svg", sep="")
            ## set format according to your needs
            svg(file = file_name)
            print(plotQualityProfile(file.path(path_to_input, fastq)))
            dev.off()
            }
            '''
    )

    #####################################################################
    #
    #####################################################################

    log.info("Filtering sequences.")

    robjects.r(
            '''
            # Filtered forward files go into the path/filtered_fwd/subdirectory
            filtpathF = file.path(path_to_input, "filtered_fwd")
            # Filtered forward files go into the path/filtered_rev/subdirectory
            filtpathR = file.path(path_to_input, "filtered_rev")
            fwd = file.path(path_to_input, fastqFs)
            filtFs = file.path(filtpathF, fastqFs)
            rev = file.path(path_to_input, fastqRs)
            filtRs = file.path(filtpathR, fastqRs)
            truncLen = c(truncate_length_fwd, truncate_length_rev)
            out = filterAndTrim(fwd = file.path(path_to_input, fastqFs),
                                filt = filtFs,
                                rev = file.path(path_to_input, fastqRs),
                                filt.rev = filtRs,
                                truncLen = c(truncate_length_fwd, truncate_length_rev),
                                rm.phix = TRUE,
                                compress = TRUE,
                                verbose = FALSE,
                                multithread = FALSE)
            '''
    )

    #####################################################################
    # error rates and remove erroneous sequences
    #####################################################################

    log.info("Removing erronious sequences.")

    robjects.r(
            '''
            # Assumes filename = samplename_XXX.fastq.gz
            sample.namesF = sapply(strsplit(basename(filtFs), "_"), `[`, 1)
            sample.namesR = sapply(strsplit(basename(filtRs), "_"), `[`, 1)
            '''
    )

    ###########################################################################
    # Sanity check. Alex, what is beeing tested here?

    #    NOT WORKING RIGHT NOW!

    ###########################################################################
    #    print robjects.r("sample.namesF")
    #    print robjects.r("sample.namesR")
    #    print robjects.r("(!identical(sample.namesF, sample.namesR))")
    # DevNote: Good enough for now.
    #    if robjects.r("(!identical(sample.names, sample.namesR))"):
    #         sys.exit("Forward and reverse files do not match.")

    robjects.r('''
                names(filtFs) = sample.names
                names(filtRs) = sample.names
                set.seed(100)

                # Learn forward error rates
                errF = learnErrors(filtFs, nread=2e6, multithread=TRUE)

                # Learn reverse error rates
                errR = learnErrors(filtRs, nread=2e6, multithread=TRUE)

                # Sample inference and merger of paired-end reads
                mergers = vector("list", length(sample.names))
                names(mergers) = sample.names
    ''')

    robjects.r('''
                for(sam in sample.names) {
                    cat("Processing:", sam, "\n")
                    derepF = derepFastq(filtFs[[sam]])
                    ddF = dada(derepF, err=errF, multithread=TRUE)
                    derepR = derepFastq(filtRs[[sam]])
                    ddR = dada(derepR, err=errR, multithread=TRUE)
                    merger = mergePairs(ddF, derepF, ddR, derepR)
                    mergers[[sam]] = merger
                    }
                rm(derepF)
                rm(derepR)
                ''')

    # DevNote: Check what is beeing removed on last line.

    # Construct sequence table and remove chimeras
    robjects.r("seqtab = makeSequenceTable(mergers)")
    robjects.r("saveRDS(seqtab, file.path(path_to_output, 'seqtab.rds'))")

    #####################################################################
    # analysis and plotting
    #####################################################################

    log.info("Removing chimeras and assigning taxonomic classification to reads.")

    robjects.r('''
                # Remove chimeras
                seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)

                # Assign taxonomy
                tax = assignTaxonomy(seqtab, file.path(path_to_silva), multithread=TRUE)

                # Write to disk
                saveRDS(seqtab, file.path( path_to_output, 'seqtab_final.rds'))                    # DevNote: Duplication?
                saveRDS(tax, file.path(path_to_output, "tax_final.rds")) # CHANGE ME ...
                ''')

    #####################################################################
    ############## Track reads through pipeline #########################
    #####################################################################

    log.info("Tracking reads through pipeline.")

    robjects.r('''
                getN = function(x) sum(getUniques(x))
                track = cbind(out, sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
                colnames(track) <- c("input", "filtered", "denoised.merged", "tabled", "nonchim")
                head(track)
                write(track, file = file.path(path_to_output, "track.csv"), sep = ",")

                # plot tracking
                track_melt = melt(track, id.vars = rownames)

                for(i in levels(track_melt$Var1)) {
                      p <- ggplot(subset(track_melt, Var1==i), aes(Var2, value,  fill = Var2)) +
                     facet_wrap(~ Var1) +
                     geom_bar(stat="identity", show_guide=FALSE)
                    ggsave(paste0("figure_",i,".svg"), p)
                    }
                ''')

    ####################################################################
    ############# plot taxonomy ########################################
    ####################################################################

    log.info("Generating taxonomic plots.")

    robjects.r('''
                samples.out = rownames(seqtab.nochim)

                # Construct phyloseq object (straightforward from dada2 outputs)
                ps = phyloseq(otu_table(seqtab.nochim,
                        taxa_are_rows=FALSE),
                        tax_table(tax))

                # richness
                richness_path <- file.path(path_to_output, "Richness")
                dir.create(file.path(richness_path))

                svg(file = file.path(richness_path, "richness.svg"))
                print(plot_richness(ps, measures=c("Shannon", "Simpson"), color=NULL) + theme_bw())
                dev.off()

#                # top 20 taxa
#                top20 = names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
#                ps.top20 = transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
#                ps.top20 = prune_taxa(top20, ps.top20)
#
#                svg(file = "top20.svg")
#                print(plot_bar(ps.top20, fill="Family") + facet_wrap(~When, scales="free_x"))
#                dev.off()
                ''')

    log.info("Analysis finished.")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:        # DevNote: Not working right now
        print('Interrupted')
        sys.exit(0)
