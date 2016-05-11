# QBB_Puerto_Rico_2016
An introduction to bioinformatics and a tutorial on RNA-seq analysis in halophilic Archaea

Knowledge of and proficiency in bioinformatics are more important than ever in this day and age in biology. With the advent of high throughput and cheap sequencing technologies huge datasets are becoming more available. In order to work and parse through these datasets, a biologist needs to have a good familiarity with UNIX environments, scripting, and stastistics. While most of the ideas are the same between Eukaryotic and microbial bioinformatic analysis, there are indeed some considerations to take and different software/pipelines to use to specifically study Bacteria and Archaea. This workshop is designed to introduce you to the field of bioinformatics and expose you to the tools available for microbial bioinformatic analysis (focusing on organisms from halophilic environments). The workshop is divided into three modules: 1) Quality control (QC) of raw sequence data; 2) RNA-seq alignment and quantitation; and 3) RNA-seq differential expression and visualization.

###The learning objectives are:
> - Become familiar with working with big data sets (“-omics”), emphasis on microbial genomics

> - Exposure to the most recent bioinformatic programs

> - Ability to work in a UNIX environment, command-line, executing programs

> - Literacy in file structure of big data sets

> - RNA-seq analysis in *Haloarchaea*

#Setting up a UNIX environment for non-Mac/Linux/Ubuntu people

insert Virtual box .img

This virtual box .img has Ubuntu v. installed, along with the whole skew of programs we need for RNA-seq analysis. If you can succesfully install this you can skip the installation steps in the Tools section.

#Tools

After setting up a Unix-based (Mac/Linux) environment, I need the following software to teach RNA-seq analysis:

> - Homebrew: [Mac Version](http://brew.sh/) [Linux Version](http://linuxbrew.sh/)
> - [Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
> - [Python 2.7.1](https://www.python.org/download/releases/2.7.1/)
> - [Pip](https://pip.pypa.io/en/stable/installing/)

*Installing the first 4 bolded software will make it very easy to install everything else. With homebrew, you can install most of the other software using `brew install <name_of_program>` by tapping into the [science homebrew github page](https://github.com/Homebrew/homebrew-science)*

**RNA-seq analysis software:**

> - [Samtools](http://www.htslib.org/)
> - [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml)
> - [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
> - [EDGEpro](http://ccb.jhu.edu/software/EDGE-pro/)
> - [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/getting_started/)
> - [Trim_galore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
> - [Stringtie](https://ccb.jhu.edu/software/stringtie/)
> - [BEDtools](http://bedtools.readthedocs.org/en/latest/content/installation.html)
> - [BLAST](http://www.ncbi.nlm.nih.gov/books/NBK279690/)
> - [Gffcompare](https://ccb.jhu.edu/software/stringtie/gff.shtml)
> - [HTSeq](http://www-huber.embl.de/HTSeq/doc/install.html#install)
> - [Salmon](http://salmon.readthedocs.org/en/latest/building.html#installation)
> - [R](https://www.r-project.org/)
> - [DEseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
> - [IGV](https://www.broadinstitute.org/software/igv/log-in) *you must register for free*
> - [Sleuth](https://github.com/pachterlab/sleuth)
> - [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

Another great resource of scientific software is found on one of my professor's, Dr. James Taylor, [github page](https://github.com/jxtx/mac-dev-playbook) ---> this is Mac specific though.

#Data

> - *Haloferax volcanii* NCBI reference genome: 
> `wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000025685.1_ASM2568v1/GCF_000025685.1_ASM2568v1_genomic.fna.gz`
> - *Haloferax volcanii* NCBI reference gene annotations: 
> `wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000025685.1_ASM2568v1/GCF_000025685.1_ASM2568v1_genomic.gff.gz`
> - *Haloferax volcanii* RNA-seq raw reads: 
> `wget `

#Literature
[Count-based differential expression analysis of RNA sequencing data using R and Bioconductor](http://www.nature.com/nprot/journal/v8/n9/pdf/nprot.2013.099.pdf)
Simon Anders, Davis J McCarthy, Yunshun Chen, Michal Okoniewski, Gordon K Smyth,
Wolfgang Huber & Mark D Robinson

#*Module 1: Introduction and Quality Control of RNA-seq .fastq Reads*
[[file:///Users/DRG/Dropbox/Screenshots/Screenshot%202016-05-06%2016.43.28.png]]

#*Module 2: RNA-seq Alignment, Transcript Assembly, and Processing*
We will be doing a reference-based RNA-seq alignment analysis today. Inherent in this, we require a reference genome to do this sort of analysis. Fortunately, the *Haloferax volcanii* genome has been sequenced. Not everyone is so lucky, and because so few genomes are actually sequenced and annotated there are reference-indepenedent approaches as well (ie *de novo* assembly of transcripts) but this is beyond the scope of this workshop.

###*Step 1: Filter out Ribosomal RNA reads*
So we've got our reads trimmed of adapters and quality filtered but there is another filtering step required: removal of rRNA reads. You might ask, 
"Why worry about rRNA?" In all of life as we know it, rRNA is the most transcribed RNA species in the cell. It can attribute for as much as 90% of total RNA in a cell! Moreover, rRNA is typically constituitively expressed at high-levels regardless of growth state (translation is occuring all the time). Because of this skew in rRNA abundance, a typical RNA-seq run will result in ~90% of total reads mapping to rRNA. Now-a-days bioinformaticians do not need worry about rRNA due to awesome techniques in removing rRNA during the library preparation stage (see Illumina Ribo-zero kit), effectively depleting rRNA levels to <1% and saving all your sequence reads to target mRNA and other RNA species! Unfortunately, these techniques are well established for Bacteria and Eukarya but not for Archaea... :( 

A couple of other reasons to remove rRNA, more important to bioinformatic analysis, is that something as highly expressed as rRNA can mask expression of more lowly expressed RNA species. So our first task in our RNA-seq adventure is to remove rRNA from our reads. We will do this by mapping our reads to 5S, 16S, and 23S rRNA. Let's pull some reference files from NCBI:

**Get rRNA fasta files**
> `Open Chrome.app and go to 

**Build rRNA hisat2 index**
> `$ hisat2-build /path/to/HFX_NCBI_rRNA.fa /path/to/index/*prefix*`

We first need to build an index of our rRNA reference sequences. The *prefix* can be anything you want, but usually something logical is best like *HFX_NCBI_rRNA*


**Align reads to rRNA & extract reads that do not align (eg mRNA)**

Next we will align our reads against the HFX rRNAs, effectively "sticking" our rRNA reads onto an alignment file (.sam), and allowing all the other reads that don't map to rRNA to "flow" pass and allow us to capture them. Since this is a paired-end RNA-seq data set we will have to map each mate pair separately (ie Read1, Read2).

Read1:
> `$ hisat2 --verbose *--un /path/to/read1_rRNA_removed.fq*`

First we are calling hisat2 on the command line and beginning to use options the program has. To see all the options run hisat2 --help. We want to set the option *--un* to return unaligned reads (ie our mRNA and other RNAs) and the pathway to the directory where we want our rRNA filtered reads to be written

> `hisat2 --verbose --un /path/to/read1_rRNA_removed.fq *--no-spliced-alignment*` 

A great option we want to specify for Prokaryotic (sorry for the *archaic* term ;) but microbial is not appropriate to say) RNA-seq analysis is to block the aligner from looking for splice sites. There is usually no splicing in Bacteria and Archaea.

> `hisat2 --verbose --un /path/to/read1.fq --no-spliced-alignment *--rna-strandness RF*`

We want to set the *--rna-strandness RF* option because this is strand-specific RNA-seq data. RNA-seq on its own does not preserve strand specificity of transcripts, but wizards have developed trickery to preserve the orientation of transcripts (ie dUTP method). How do we know whether to specificy RF or FR? Follow the rabbit hole of manual citations (hisat2 manual --> tophat2 manual).
 
> `hisat2 --verbose --un /path/to/read1_rRNA_removed.fq --no-spliced-alignment --rna-strandness RF *--dta -I 0 -X 500*`

Some more important options: *--dta* is an option to put to make sure the output data is in the correct format to be piped into a transcript assembler stringtie (which we will be using). *-I and -X* are options specific to the insert sizes of the reads. This is information you get from library prepartion (ie the average fragment size of your RNA). For this data, the average insert size is ~350 bp. Setting a bigger max insert size usually yields in better alignments but causes slower aligning. 

> `hisat2 --verbose --un /path/to/read1_rRNA_removed.fq --no-spliced-alignment --rna-strandness RF --dta -I 0 -X 500 *-x /path/to/hisat2/index/prefix -U /path/to/read1.fq*`

Next we set the path to the hisat2 index we built earlier and the path to your first mate reads (Read1).

> `hisat2 --verbose --un /path/to/read1_rRNA_removed.fq --no-spliced-alignment --rna-strandness RF --dta -I 0 -X 500 -x /path/to/hisat2/index/prefix -U /path/to/read1.fq *-S /path/to/rRNA_alignment_R1.sam*`

Lastly we'll set the output pathway to where we want our alignment file to be put (.sam file)

Do the same for read2:
> `$ hisat2 --verbose --un /path/to/read2_rRNA_removed.fq --no-spliced-alignment --rna-strandness RF --dta -I 0 -X 500 -x /path/to/hisat2/index/prefix -U /path/to/read2.fq -S /path/to/rRNA_alignment_R2.sam`


**Synchronize paired-end files (thus removing singlets)**
> `$ python Downloads/Scripts-master/fastqCombinePairedEnd.py /Users/DRG/Desktop/HFX_rRNA/rRNA_removed_trimmed_reads/HFX_C1_IR_rRNA_removed_reads_1.fq /Users/DRG/Desktop/HFX_rRNA/rRNA_removed_trimmed_reads/HFX_C1_IR_rRNA_removed_reads_2.fq`


###*Step 2: Align rRNA-filtered reads against NCBI reference genome*
Now that we have rRNA-clean reads, let us do the real alignment.

**Build hisat2 NCBI refseq index**
> `$ hisat2-build /path/to/HFX_NCBI_reference_genome.fa`


**Align filtered reads**
> `$ hisat2 --verbose  --no-spliced-alignment --rna-strandness RF --dta -I 0 -X 500 -x /path/to/hisat2/HFX_index -1 /path/to/synced/rRNA-removed_read1.fq -2 /path/to/synced/rRNA-removed_read2.fq -S /path/to/rRNA-removed_hisat2_alignment.sam`

You now have your first alignment file, which is in the .sam format. You can learn more about the format here: https://samtools.github.io/hts-specs/SAMv1.pdf

> #Task #2

###*Step 3: Assemble transcripts from aligned reads and quantitate*

**Convert sam to bam**
> `$ samtools view -bS /Users/DRG/Desktop/HFX_rRNA/rRNA_removed_alignments/HFX_O3_mRNA_rRNA_mapped_hisat2_rRNA_removed_alignment.sam > /Users/DRG/Desktop/HFX_rRNA/rRNA_removed_alignments/HFX_O3_mRNA_rRNA_mapped_hisat2_rRNA_removed_alignment.bam`

**Sort by coordinate & by name**

coordinate:

> `$ samtools sort /Users/DRG/Desktop/HFX_rRNA/rRNA_removed_alignments/bam_files/HFX_O1_mRNA_rRNA_mapped_hisat2_rRNA_removed_alignment.bam -o /Users/DRG/Desktop/HFX_rRNA/rRNA_removed_alignments/bam_files/HFX_O1_mRNA_rRNA_mapped_hisat2_rRNA_removed_alignment.cordsorted.bam`

name:
> `$ samtools sort -n /Users/DRG/Desktop/HFX_rRNA/rRNA_removed_alignments/bam_files/HFX_O1_mRNA_rRNA_mapped_hisat2_rRNA_removed_alignment.bam -o /Users/DRG/Desktop/HFX_rRNA/rRNA_removed_alignments/bam_files/HFX_O1_mRNA_rRNA_mapped_hisat2_rRNA_removed_alignment.namesorted`


**Assemble reads into transcripts and quantify abundance**
> `$ stringtie /Users/DRG/Desktop/HFX_rRNA/rRNA_removed_alignments/bam_files/cordsorted_bam/HFX_C1_mRNA_rRNA_mapped_hisat2_rRNA_removed_alignment.cordsorted -G /Users/DRG/Desktop/sRNA_in_Archaea/Data/RNA-seq1-H2O2/EDGEpro_out_HFX/HFX_genome.gff -o /Users/DRG/Desktop/HFX_rRNA/rRNA_removed_alignments/stringtie_out/HFX_C1_rRNA_removed_stringtie_out/HFX_C1_rRNA_removed_transcriptome_stringout.gtf -A /Users/DRG/Desktop/HFX_rRNA/rRNA_removed_alignments/stringtie_out/HFX_C1_rRNA_removed_stringtie_out/HFX_C1_rRNA_removed_gene_abund.tab -C /Users/DRG/Desktop/HFX_rRNA/rRNA_removed_alignments/stringtie_out/HFX_C1_rRNA_removed_stringtie_out/HFX_C1_rRNA_removed_cov_refs.gtf -m 30`

**Merge transcriptomes**
> `$ cuffcompare -o /Users/DRG/Desktop/HFX_rRNA/rRNA_removed_alignments/stringtie_out/cuffcompare/HFX_IR_RNAseq3_rRNA_removed_cuffcompare -i /Users/DRG/Desktop/HFX_rRNA/rRNA_removed_aligments/all_gtf_out_list.txt  -r /Users/DRG/Desktop/sRNA_in_Archaea/Data/RNA-seq1-H2O2/EDGEpro_out_HFX/HFX_genome.gff`



#*Module 3: RNA-seq Differential Expression Analysis and Visualization*



