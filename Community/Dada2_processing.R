setwd('/Users/periwal/ShikiFactory/WP3/SW_SF100/Community/')

##loading packages

library(dada2)
library(DECIPHER)
library(ShortRead)
library(Biostrings)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(phyloseq)

set.seed(100)

### Reading demultiplexed files and primers/adapters

#The 16S sequencing data comprise of paired-end sequencing reads from V4 region of the 
#bacterial genomes.Plot qualtiy of first few reads and identify trimming parameters. 
#The adapters were checked with cutadapt tool but there were no adaptors found.
  
#reading reads
  
path<-"./sweet_community/16Sreads"
list.files(path)

forward_reads<-sort(list.files(path, pattern = "_1_sequence.txt.gz", full.names = TRUE))
reverse_reads<-sort(list.files(path, pattern = "_2_sequence.txt.gz", full.names = TRUE))

plotQualityProfile(forward_reads[1:2])
plotQualityProfile(reverse_reads[1:2])

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_1_sequence")[[1]][1]
sample.names <- unname(sapply(forward_reads, get.sample.name))
head(sample.names)

### Filter and trim reads (Denoising)

#Quality plot showed read quality drops at ends (F-230,R-240), thus we discard reads with 
#lower quality score. Standard filtering parameters: 
#maxN=0 (DADA2 requires no Ns) are used.
  
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#takes very long time (30min+)
out <- data.frame(filterAndTrim(forward_reads, filtFs, reverse_reads, filtRs,
                                maxN=0,truncLen=c(230,240),maxEE=2,multithread=TRUE, matchIDs=TRUE))
head(out)

#percentage of reads filtered
out_percent<-out %>% mutate(percent_filtered = reads.out/reads.in)
#View(out_percent)

### Learn error rates

errF <- learnErrors(filtFs, multithread=TRUE)
#107548920 total bases in 467604 reads from 7 samples will be used for learning the error rates.
errR <- learnErrors(filtRs, multithread=TRUE)
#112224960 total bases in 467604 from 7 samples will be used for learning the error rates.

#plot errors
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

### Sample inference

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#dadaFs[[20]]
#dadaRs[[20]]

### Merge paired reads

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

### Construct sequence table

seqtab <- makeSequenceTable(mergers)
#head(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

### Remove chimeras

#The core dada method corrects substitution and indel errors, but chimeras remain. 
#Chimeric sequences are identified if they can be exactly reconstructed by combining
#a left-segment and a right-segment from two more abundant “parent” sequences.
  
#used 'per-sample' method to remove chimeras, other methods discards majority of the sequences as chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="per-sample", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
head(seqtab.nochim)
sum(seqtab)
sum(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

### Track reads through the pipeline

#As a final check of our progress, we’ll look at the number of reads that made it 
#through each step in the pipeline:
  
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim), rowSums(seqtab.nochim)/sapply(mergers, getN))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim","%")
rownames(track) <- sample.names
View(track)

### Assign taxonomy

ref<-"~/ShikiFactory/WP3/SW_SF100/Community/sweet_community/16S_25bugs.fa"

set.seed(999)
taxa <- assignTaxonomy(seqtab.nochim, ref, multithread = TRUE, tryRC = TRUE, minBoot = 50)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
View(taxa.print)

#sample metadata
samples.out <- rownames(seqtab.nochim)
samples.out

samdf <- data.frame(read.table(file = '~/ShikiFactory/WP3/SW_SF100/Community/sweet_community/metadata_sampbiol.tsv', header = TRUE, sep = '\t'))
rownames(samdf) <- samples.out
head(samdf)

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

###save files
asv_table<-data.frame(t(otu_table(ps)))
#View(asv_table)
df0 <- tibble::rownames_to_column(asv_table, "ASV")
#View(df0)
write.table(df0, file = "~/ShikiFactory/WP3/SW_SF100/Community/sweet_community/ASV_25bugs.tsv", sep = '\t', row.names = FALSE)

#change sample names of ASV_25bugs in text editor

taxa_table<-data.frame(tax_table(ps))
df <- tibble::rownames_to_column(taxa_table, "ASV")
head(df)
write.table(df, file = "~/ShikiFactory/WP3/SW_SF100/Community/sweet_community/TAXA_25bugs.tsv", sep = '\t', row.names = FALSE)
