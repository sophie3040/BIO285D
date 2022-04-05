
#S¿THIS ARE THE LIBRARIES FOR ALL THE WORKFLOW
library(ggplot2)
library(phyloseq); packageVersion("phyloseq")
library(ShortRead)
library(dada2)
library(ape) #library for creating  tree
library(dplyr)
library(vegan)
library(DECIPHER)
library(rlang)
library(ggplot2); packageVersion("ggplot2")
library(msa)
library(phangorn)
library(microbiome) # data analysis and visualisation
library(microbiomeutilities) # some utility tools 
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(DT) # interactive tables in html and markdown
library(data.table) # alternative to data.frame
library(picante)
library("scales")
library("grid")
library(tidyverse)
library(BiocGenerics)


####DADA2#########
## Importe de secuenciasd para analisis de calidad

```{r, echo=TRUE, message=FALSE, warning=FALSE}
path <- "C:/Users/sophi/Desktop/Tesis/Dataset/petB/"
fns <- list.files(path)
fastqs <- fns[grepl(".fastq", fns)] 
fastqs <- sort(fastqs)
fnFs <- fastqs[grepl("_R1.fastq", fastqs)] # Fws
fnRs <- fastqs[grepl("_R2.fastq", fastqs)] # Rvs

# Get sample names from the first part of the forward read filenames
#sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names <- c("LOM_70", "LOM_71", "LOM_72", "LOM_73", "LOM_74", "LOM_75", "LOM_76", "LOM_77", "LOM_78", "LOM_79", "LOM_80", "LOM_81", "LOM_82", "LOM_83", "LOM_84")
sample.names <- sapply(strsplit(basename(fnFs), "_"), '[',3)
str(sample.names)
# Para especificar de manera concreta donde se encuentras las secuencias Fw y Rv
fnFs <- paste0(path, fnFs)
fnRs <- paste0(path, fnRs)
```

#Finally, we visualize the quality profiles of the reads:
# In **gray-scale** is a heat map of the frequency of each quality score at each base position. The mean quality score at each position is shown by the **green line**, and the quartiles of the quality score distribution by the **orange lines**. The **red line** shows the scaled proportion of reads that extend to at least that position (this is more useful for other sequencing technologies, as Illumina reads are typically all the same lenghth, hence the flat red line).

#### Quality profiles visualization
```{r, echo=TRUE, message=FALSE, warning=FALSE}
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
```

## Sequence filter and trimming
#Si más adelante aparecen "problemas" o se aprecian ciertas irregularidades no es mala opción evaluar la exclusión de las muestras "LOM 78" y "LOM 79" ya que en el proceso podrían mejorar el rendimiento de las funciones de dada2.

#Lo siguiente varía según tu criterio, de acuerdo a lo que yo aprecio un largo de 230 nucleotidos para el marco de lectura "Forward" significa una calidad buena de secuencias. Para las secuencias de marco de lectura complementario un tamaño con un puntaje de calidad adecuado parece estar entre los 200 y 210 nucleotidos. Finalmente realizaremos un corte al comienzo y termino de las secuencias buscando estabilizar la curva de calidad de las secuencias (10 - 20 nucleotidos). 

#Based on quality profiles we can trim **NTs** from left (*trimLeft*) and anything below **Q30** (*truncLen*) on the right. Also we can set filtering parameters according to our needs (computing power, strictness and so on).

```{r, echo=TRUE, message=TRUE, warning=TRUE}
ptm <- proc.time()
filtpth <- file.path(path, "Filtered_")
filtFs <- paste0(filtpth, sample.names, "_F_filt.fastq.gz")
filtRs <- paste0(filtpth, sample.names, "_R_filt.fastq.gz")

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(230,200), trimLeft=c(15,15), 
                     maxN=0, maxEE=c(2,2), truncQ=2, 
                     compress=TRUE, multithread=FALSE, verbose=TRUE)

rinout <- as.data.frame(out)
write.csv2(rinout, file = "reads_out.csv")

proc.time() - ptm
```

#The *filterAndTrim(...)* function filters the forward and reverse reads jointly, outputting only those pairs of reads that both pass the filter. In this function call we did four things: We removed the first *trimLeft=10* nucleotides of each read. We truncated the forward and reverse reads at *truncLen=c(240, 200)* nucleotides respectively. We filtered out all reads with more than *maxN=0* ambiguous nucleotides. And we filtered out all reads with more than two expected errors. The filtered output files were stored as gzipped fastq files *(compress=TRUE)*.

#This represents a fairly standard set of filtering/trimming parameters. However, it is always worth evaluating whether the filtering and trimming parameters you are using are appropriate for your data. One size does not fit all! (And are you sure you have removed your primers?)

#An important consideration: If using paired-end sequencing data, you must maintain a suitable overlap (>20nts) between the forward and reverse reads after trimming! This is especially important to keep in mind for mult-V-region amplicions (such as V3-V4) in which there may be relatively little overlap to begin with, and thus little read-truncation is possible if reads are to be merged later on.


## Dereplicate

#The next thing we want to do is “dereplicate” the filtered fastq files. Finding the set of unique sequences, equivalently, the process of finding duplicated (replicate) sequences. Dereplication combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance” equal to the number of reads with that unique sequence. During dereplication, we condense the data by collapsing together all reads that encode the same sequence, which significantly reduces later computation times.

#Dereplication in the DADA2 pipeline has one crucial addition from other pipelines: DADA2 retains a summary of the quality information associated with each unique sequence. The consensus quality profile of a unique sequence is the average of the positional qualities from the dereplicated reads. These quality profiles inform the error model of the subsequent sample inference step, significantly increasing DADA2’s accuracy.

```{r,echo=TRUE, message=FALSE, warning=FALSE}
ptm <- proc.time()
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names 
names(derepRs) <- sample.names
proc.time() - ptm
```


## Error rate estimation
Ojo que puedes de-replicar previo a este punto (como lo hice yo), o posteriormente. Lo he visto de las formas y no se generan diferencias perceptibles en el analisis posterior de las secuencias.

Generalmente el calculo de la taza de error se realiza leyendo las muestras en el orden dado hasta lograr un numero suficiente de lecturas. Yo creo que es mejor que el calculo lecturas obtenidas de forma aleatoria desde el "pool" de secuencias.

The DADA2 algorithm makes use of a parametric error model *(err)* and every amplicon dataset has a different set of error rates. The $learnErrors$ method learns the error model from the data, by alternating estimation of the error rates and inference of sample composition until they converge on a jointly consistent solution. As in many machine-learning problems, the algorithm must begin with an initial guess, for which the maximum possible error rates in this data are used (the error rates if only the most abundant sequence is correct and all the rest are errors).

```{r, echo=TRUE, message=FALSE, warning=FALSE}
ptm <- proc.time()
errF <- learnErrors(filtFs, multithread=TRUE, randomize = TRUE)
errR <- learnErrors(filtRs, multithread=TRUE, randomize = TRUE)
proc.time() - ptm
```

It is always worthwhile to visualize the estimated error rates:
  
  #### Error rate graphs.
  
  ```{r, echo=TRUE, message=FALSE, warning=FALSE}
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
```
The error rates for each possible transition (A→C, A→G, …) are shown. **Points** are the observed error rates for each consensus quality score. The **black line** shows the estimated error rates after convergence of the machine-learning algorithm. The **red line** shows the error rates expected under the nominal definition of the Q-score.


## Joint sample inference and error rate estimation.

```{r, echo=TRUE, message=FALSE, warning=FALSE}
ptm <- proc.time()
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
proc.time() - ptm
```

By default, the dada function processes each sample independently. However, pooling information across samples can increase sensitivity to sequence variants that may be present at very low frequencies in multiple samples. The dada2 package offers two types of pooling. *dada(..., pool=TRUE)* performs standard pooled processing, in which all samples are pooled together for sample inference. *dada(..., pool="pseudo")* performs pseudo-pooling, in which samples are processed independently after sharing information between samples, approximating pooled sample inference in linear time.


## Merge paired reads.

We now merge the forward and reverse reads together to obtain the full denoised sequences. Merging is performed by aligning the denoised forward reads with the reverse-complement of the corresponding denoised reverse reads, and then constructing the merged “contig” sequences. By default, merged sequences are only output if the forward and reverse reads overlap by at least 12 bases, and are identical to each other in the overlap region (these conditions can be changed via function arguments). Non-overlapping reads are supported, but not recommended

```{r, echo=TRUE, message=FALSE, warning=FALSE}
mergers <- mergePairs(dadaFs, derepFs, verbose=TRUE)
class(mergers) #list
length(mergers) #15 elements in the list, one for each sample
names(mergers)
head(mergers[[1]])
```

The mergers object is a list of data.frames from each sample. Each data.frame contains the merged *$sequence*, its *$abundance*, and the indices of the *$forward* and *$reverse* sequence variants that were merged. Paired reads that did not exactly overlap were removed by *mergePairs*, further reducing spurious output.


## Construct ASV table

The sequence table is a *matrix* with rows corresponding to (and named by) the samples, and columns corresponding to (and named by) the sequence variants.

```{r, echo=TRUE, message=FALSE, warning=FALSE}
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab))) # Inspect distribution of sequence lengths
```


## Remove chimeras

#The core dada method corrects substitution and indel errors, but chimeras remain. Fortunately, the accuracy of sequence variants after denoising makes identifying chimeric ASVs simpler than when dealing with fuzzy OTUs. Chimeric sequences are identified if they can be exactly reconstructed by combining a left-segment and a right-segment from two more abundant “parent” sequences.
#The frequency of chimeric sequences varies substantially from dataset to dataset, and depends on factors including experimental procedures and sample complexity.
#Most of your **reads** should remain after chimera removal (it is not uncommon for a majority of **sequence variants** to be removed though). If most of your reads were removed as chimeric, upstream processing may need to be revisited. In almost all cases this is caused by primer sequences with ambiguous nucleotides that were not removed prior to beginning the DADA2 pipeline.
#By default the *method* is set in **"consensus"** where the samples in a sequence table are independently checked for bimeras, and a **consensus decision on each sequence variant is made**. If it sets as **"per-sample"**, samples in a sequence table are independently checked for bimeras, **and sequence variants are removed (zeroed-out) from samples independently**. Another alternative is **"pooled"** the samples fot bimera identification.

```{r, echo=TRUE, message=FALSE, warning=FALSE}
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) #how much have been lost in the process
```


## Track reads through the pipeline

#As a final check of our progress, we’ll look at the number of reads that made it through each step in the pipeline:
  ```{r, echo=TRUE, message=FALSE, warning=FALSE}
getN <- function(x) sum(getUniques(x))
summary_tab <- data.frame(row.names = )
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nochim")
rownames(track) <- sample.names
print(track)
trackdf <- as.data.frame(track)
write.csv2(trackdf, file = "process_track.csv")
#Alt
Samples <- sample.names
summary_tab <- data.frame(row.names = Samples, 
                          dada2_input = out[, 1], 
                          filtered = out[, 2], 
                          dada_Fw = sapply(dadaFs, getN), 
                          dada_Rv = sapply(dadaRs, getN),
                          merged = sapply(mergers, getN),
                          nonchim = rowSums(seqtab.nochim),
                          final_perc_reads_retained = round(rowSums(seqtab.nochim)/out[, 1]*100, 1)
)
write.csv2(summary_tab, file = "summary_tab_dada2.csv")
#If error in some point, go back to those and roll again
```

Outside of filtering, there should no step in which a majority of reads are lost. If a majority of reads failed to merge, you may need to revisit the *truncLen* parameter used in the filtering step and make sure that the truncated reads span your amplicon. If a majority of reads were removed as chimeric, you may need to revisit the removal of primers, as the ambiguous nucleotides in unremoved primers interfere with chimera identification.

#### Remove sequence variants seen less than a given number of times 
```{r, echo=TRUE, message=FALSE, warning=FALSE}
seqtab.nochim = seqtab.nochim[,colSums(seqtab.nochim) > 10]
# Opcional, de acuerdo a tus criterios y teniendo en cuenta el impacto que tendrán en tus resultados (el porqué)
```


## Write sequence table to file
```{r, echo=TRUE, message=FALSE, warning=FALSE}
write.csv2(seqtab.nochim, file = "seqtab_nochim.csv")
```


## Taxonomy assignation (Silva, RDP, Greengenes, etc.)

It is common at this point, especially in 16S/18S/ITS amplicon sequencing, to assign taxonomy to the sequence variants. The DADA2 package provides a native implementation of the naive Bayesian classifier method for this purpose. The *assignTaxonomy* function takes as input a set of sequences to be classified and a training set of reference sequences with known taxonomy, and outputs taxonomic assignments with at least *minBoot* bootstrap confidence.

```{r Taxonomy assignation, echo=TRUE, message=FALSE, warning=FALSE, include=FALSE}
taxa <- assignTaxonomy(seqtab.nochim, "Euk_database.fa", taxLevels = c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species"), multithread = TRUE)
```


```{r Inspection, echo=TRUE,message=FALSE, warning=FALSE, include=TRUE}
# Removing sequence rownames for display only
strict <- taxa
rownames(strict) <- NULL
head(strict)
```


## Write taxonomy assignments to file (PR2)

### PR2 4.14.0
```{r Taxonomy table, message=FALSE, warning=FALSE, include=TRUE}
write.csv2(taxa, file = "taxa_Bact.csv")
```

###Extracting the goods(to differentiate from usual method goods)
```{r Goods, echo=TRUE, message=TRUE, warning=TRUE, include=TRUE}
#giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochin)
asv_headers <- vector(dim(seqtab.nochin)[2], mode = "character")

for (i in 1:dim(seqtab.nochin)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep = "_")
}

#making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "seqs_Bact_lom.fasta")

#count table
asv_tab <- t(seqtab.nochin) %>% as.data.frame()
rownames(asv_tab) <- sub(">", "", asv_headers)
colnames(asv_tab) <- sample.names
asv_tab <- mutate(asv_tab, Seqs = all_of(colnames(seqtab.nochin)), .before = "F1")
write.csv2(asv_tab, file = "ASVs_counts.csv")
```
asv_tab


### END OF DADA2###


####PHYLOSEQ####


