# Python scripts associated with yam genome assembly

Siadjeu C.+, **Pucker B.+**, Viehoever P., Albach D. and Weisshaar B. (2020). High contiguity de novo genome sequence assembly of Trifoliate yam (Dioscorea dumetorum) using long read sequencing. Genes. doi:[10.3390/genes11030274](https://www.mdpi.com/2073-4425/11/3/274/htm).


Generation of a stacked barplot based on the summary string provided by [BUSCO](https://busco.ezlab.org/).

python BUSCO_plot.py \
--in <INPUT_FILE_WITH_BUSCO_SUMMARY_STRINGS> \
--fig <FIGURE_FILE (PDF)>


Analyses a coverage file and produces a histogram to illustrate the frequency of different coverage values.

python analyze_cov_per_contig.py \
--cov <COVERAGE_FILE> \
--out <OUTPUT_FOLDER>


Repeat maskings with different reference databases can be combined. A base is masked if it is masked in at least one input file.

python combine_repeat_masking.py \
--in1 <INPUT_FILE2> \
--in2 <INPUT_FILE1> \
--out <OUTPUT_FILE>


A histogram is produced based on a coverage file.

python cov2hist.py \
--cov <COVERAGE_FILE> \
--fig <FIGURE_FILENAME>


A coverage histogram is produced, but using the average value per contig instead of per base values.

python cov_hist_per_contig.py \
--cov <COVERAGE_FILE> \
--out <OUTPUT_FOLDER>


Selection of genes (coding sequences and peptides) after homology_based_gene_selection.py was applied. The results of different cutoff values can be evaluated to find the most appropriate settings.

python gene_selector.py \
--list <HOMOLOGY_BASED_ANALYSIS_RESULTS (TXT)> \
--gff <INPUT_FILE (FASTA)> \
--out <OUTPUT_FOLDER> \
optional: \
--fasta <MATCHING_ASSEMBLY_FILE> \
--name <NAME_OF_OUTPUT_FILES>["x"] \
--sr <FLOAT, score ratio cutoff>[0.1] \
--cov <FLOAT, percentage of query aligned>[0.1]


This scripts performs BLAST-based comparison of genes on the peptide level to identify similarity of predicted gene products to existing annotations in databases. The score of hits against the database is normalized by the score against the sequence itself.

python homology_based_gene_selection.py \
--query <INPUT_FILE (FASTA)> \
--subject <INPUT_FILE (FASTA)> \
--out <OUTPUT_FOLDER>


Functional annotations can be transferred from Arabidopsis thaliana ([Araport11](https://www.araport.org/)) after running [OrthoFinder2](https://github.com/davidemms/OrthoFinder) for the identification of orthologs.

python map_annotation_via_orthologs.py \
--in <INPUT_FILE> \
--out <OUTPUT_FILE>
      

Average coverage per gene is calculated and visualized in a histogram. Based on the known average sequecing depth, genes can be classified into 'phased' and 'merged'.

python phasing_status_per_gene.py \
--cov <COVERAGE_FILE> \
--gff <ANNOTATION_FILE> \
--out <OUTPUT_FOLDER> \
--vcf <VCF_FILE>


Contig names can be reduced to tigXXX where XXX represents the number. This allows to remove all additional entries from the header which might have accumulated during the assmbly and polishing process.
WARNING: All contig numbers may only occur ones in the file. Otherwise, contigs would be lost.

python reduce_to_tig.py \
--in <INPUT_FILE> \
--out <OUTPUT_FILE>



Additional scripts are located here:

https://github.com/bpucker/script_collection

https://github.com/bpucker/variant_calling

https://github.com/bpucker/At7

https://github.com/bpucker/banana

https://github.com/bpucker/GenomeAssemblies2018

https://github.com/bpucker/Nd1_PacBio



# References

Siadjeu C.+, **Pucker B.+**, Viehoever P., Albach D. and Weisshaar B. (2020). High contiguity de novo genome sequence assembly of Trifoliate yam (Dioscorea dumetorum) using long read sequencing. Genes. doi:[10.3390/genes11030274](https://www.mdpi.com/2073-4425/11/3/274/htm).

Pucker, B., Holtgräwe, D., Rosleff Sörensen, T., Stracke, R., Viehöver, P., and Weisshaar, B. (2016). A de novo Genome Sequence Assembly of the Arabidopsis thaliana Accession Niederzenz-1 Displays Presence/Absence Variation and Strong Synteny. PloS-ONE 11:e0164321. doi:[10.1371/journal.pone.0164321](https://doi.org/10.1371/journal.pone.0164321).

Pucker, B. and Brockington, S.F. (2018). Genome-wide analyses supported by RNA-Seq reveal non-canonical splice sites in plant genomes. BMC Genomics. 2018;19(1). doi:[10.1186/s12864-018-5360-z](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-5360-z).

Baasner, J.-S., Howard, D., Pucker, B. (2019). Influence of neighboring small sequence variants on functional impact prediction. bioRxiv. doi:[10.1101/596718](https://www.biorxiv.org/content/10.1101/596718v2.full).

Pucker, B., Holtgraewe, D., Stadermann, K. B., Frey, K., Huettel, B., Reinhardt, R., & Weisshaar, B. (2019). A Chromosome-level Sequence Assembly Reveals the Structure of the Arabidopsis thaliana Nd-1 Genome and its Gene Set. PLOS ONE: e0216233. doi: [10.1371/journal.pone.0216233](https://doi.org/10.1371/journal.pone.0216233).

Pucker, B., Feng, T., Brockington, S. (2019). Next generation sequencing to investigate genomic diversity in Caryophyllales. bioRxiv 646133; doi:[10.1101/646133](https://www.biorxiv.org/content/10.1101/646133v2.full).

Pucker B, Rückert C, Stracke R, Viehöver P, Kalinowski J, Weisshaar B. Twenty-Five Years of Propagation in Suspension Cell Culture Results in Substantial Alterations of the Arabidopsis Thaliana Genome. Genes. 2019. doi:[10.3390/genes10090671](https://www.mdpi.com/2073-4425/10/9/671/htm).
