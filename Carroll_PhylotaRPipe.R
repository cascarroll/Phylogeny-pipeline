rm(list = ls())

setwd("~/Desktop/skipperPhylo")
#---------------------------------------------------------------
library(phylotaR)
library(seqinr)

# https://docs.ropensci.org/phylotaR/articles/phylotaR.html

#----------------------
# Hesperiinae : NCBI #: 40096
wd2 <- './Hesperiinae' 

#Only need to run below once, it sets up the directories in the Hesperiinae folder that we made
#setup(wd = wd2,
#      txid = 40096,
#        ncbi_dr = '/usr/local/ncbi/blast/bin/')
#run(wd = wd2)

Hesperiinae <- read_phylota(wd2)
print(Hesperiinae)
#-----------------------------------------------------------------
#     --- finding clusters---
Hesperiinae_all_clusters <- Hesperiinae
print(Hesperiinae_all_clusters)

cids <- Hesperiinae_all_clusters@cids 
#       cids are our 'cluster IDs'
#       like in python, they begin at zero
#       but R starts at 1, so call zero at 1
cids #to confirm that R starts cid '0' at 1.

n_taxa <- get_ntaxa(phylota = Hesperiinae_all_clusters, cid = cids)
#       n_taxa tells us the number of taxa per cluster,
# we range from 1 to 784
sort(unique(n_taxa))
# if we want to see the spread (since we have 570 clusters):
hist(n_taxa) # most are 0-50.
#since we have so many clusters and a max of 784 taxa in a cluster...
# let's drop any clusters with fewer than 100 taxa

keep <- cids[n_taxa > 99] # this drops any cluster with less than 4 taxa
selected <- drop_clstrs(phylota = Hesperiinae_all_clusters, cid = keep)
selected_summary <- summary(selected)
selected_summary
#       you'll notice the output is for the seed, number of sequences, and so on
# for Hesperiinae, 11 clusters are kept (see summary output below)
#-----------------------------------------------------------------------------------------------
#     selected_summary:
# ID    Type       Seed  Parent N_taxa N_seqs Med_sql       MAD
# 1   0 subtree KT125949.1   40096    784   7690     658 0.4518462
# 2   1 subtree KT125949.1 2839394    430   5193     658 0.4448292
# 3   9 subtree MF555259.1   40096    276    434    1066 0.7741341
# 4  12 subtree KT582671.1   40096    222    361     400 0.9342155
# 5  16 subtree KY028554.1   40096    159    231     610 0.9684270
# 6  20 subtree KY028371.1   40096    142    197     411 0.9817333
# 7  22 subtree KY027495.1   40096    136    195     667 0.9055327
# 8  23 subtree MF555345.1   40096    163    192     688 0.7912138
# 9  25 subtree KY045515.1   40096    145    173     849 0.9156953
# 10 30 subtree KY027732.1   40096    147    154     600 0.9372189
# 11 38 subtree KX947181.1   40096    119    121     555 0.8762414
# Definition     Feature
# 1  cytochrome (0.09), oxidase (0.09) barcode (1)
# 2     cytochrome (0.09), gene (0.09) barcode (1)
# 3      cds (0.09), elongation (0.09)           -
#   4           cds (0.1), partial (0.1)           -
#   5            cds (0.09), gene (0.09)           -
#   6            cds (0.09), gene (0.09)           -
#   7              cds (0.1), gene (0.1)           -
#   8     cds (0.1), dehydrogenase (0.1)           -
#   9               cad (0.1), cds (0.1)           -
#   10    cds (0.1), dehydrogenase (0.1)           -
#   11      arginine (0.1), argkin (0.1)           -

#     { some important things provided here }
# > median sequence length
# > mean alignment density (MAD)
# > common words in sequence description
# > feature names
#-----------------------------------------------------------------------------------------------

# Now that clusters are parsed, we can look at each cluster by calling it by its number
cid0 <- selected_summary[1, 'ID'] #I'd just like to point out: 
#                                 we are using python (and/or every language but R)
#                                 notation ; so the first cluster has ID == 0
#                                 BUT! since this is an R package, we call cid 0 as 1.

#we are going to extract CID 0 (ID=1) and CID 1 (ID=2) from the Hesperiinae set
#both are cytochrome sequences
cluster_record <- selected@clstrs[[cid0]]
seq_records <- selected@sqs[cluster_record@sids]
seq_record <- seq_records[[seq_records@ids[[1]]]] # this gives us a single record ID
summary(seq_record)
seq <- rawToChar(seq_record@sq)
reduced <- drop_by_rank(phylota = selected, rnk = 'species', n = 10) #10 is the default
txids <- get_txids(phylota = reduced, cid = cid0, rnk = 'species')
taxonomy <- get_tx_slot(phylota = reduced, txid = txids, slt_nm = 'scnm')
# cleaning taxonomy so it's cleaner
taxonomy <- gsub('\\.','', taxonomy)
taxonomy <- gsub('\\s+', '_', taxonomy) # better for linux (rm whitespace)
print(taxonomy)
# sequence IDs (for a given cluster)
sids <- reduced@clstrs[[cid0]]@sids
#write_sqs(phylota = reduced,
#            sid = sids,
#            sq_nm = taxonomy,
#          outfile = file.path('./Hesperiinae', 'Hesperiinae_cid0.fasta'))


#Doing the exact same thing as above for our cid 1 (ID=2)
cid1 <- selected_summary[2, 'ID'] #I'd just like to point out: 
#                                 we are using python (and/or every language but R)
#                                 notation ; so the first cluster has ID == 0
#                                 BUT! since this is an R package, we call cid 1 as 2.

#cluster_record1 <- selected@clstrs[[cid1]]
#seq_records1 <- selected@sqs[cluster_record1@sids]
#seq_record1 <- seq_records1[[seq_records@ids[[2]]]] # this gives us a single record ID
#summary(seq_record1)
seq1 <- rawToChar(seq_record@sq)
reduced1 <- drop_by_rank(phylota = selected, rnk = 'species', n = 10) #10 is the default
txids1 <- get_txids(phylota = reduced1, cid = cid1, rnk = 'species')
taxonomy1 <- get_tx_slot(phylota = reduced1, txid = txids, slt_nm = 'scnm')
# cleaning taxonomy so it's cleaner
taxonomy1 <- gsub('\\.','', taxonomy)
taxonomy1 <- gsub('\\s+', '_', taxonomy) # better for linux (rm whitespace)
print(taxonomy1)
# sequence IDs (for a given cluster)
sids1 <- reduced1@clstrs[[cid1]]@sids
write_sqs(phylota = reduced,
          sid = sids1,
          sq_nm = taxonomy,
          outfile = file.path('./Hesperiinae', 'Hesperiinae_cid1.fasta'))



#---------------------------------------------------------------------------------

library(ape)
library('phangorn')

Hesperiinae0 <- read.dna('./Hesperiinae/aliHesperiinae_cid0.fst', format='fasta')
print(Hesperiinae0)
Hesperiinae0_phyDat <- phyDat(Hesperiinae0, type = "DNA", levels = NULL)
phyDat

#devtools::install_github('Bioconductor/Biostrings') -------- [ Use Biostrings to read in fasta generated in wd/blast]
library(Biostrings)
Hesperiinae_fasta <- readDNAStringSet('./Hesperiinae/blast/taxon-40096-typ-subtree-db.fa')
Hesperiinae_fasta
#--------------------------------------------------------------
# Pipeline > open in seaview
# I opened up the files we write at lines 104 and 130 in seaview for alignment