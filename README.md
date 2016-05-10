# QBB_Puerto_Rico_2016
An introduction to bioinformatics and a tutorial on RNA-seq analysis in halophilic Archaea

Knowledge and proficiency in bioinformatics is more important than ever in this day and age in biology. With the advent of high throughput and cheap sequencing technologies huge datasets are becoming more available. In order to work and parse through these datasets, a biologist needs to have a good familiarity with UNIX environments, scripting, and stastistics. While most of the ideas are the same between Eukaryotic and microbial bioinformatic analysis there are indeed some considerations to take and different software/pipelines to use to specifically study Bacteria and Archaea. This workshop is designed to introduce you to the field of bioinformatics and expose you to the tools available for microbial bioinformatic analysis (focusing on organisms from halophilic environments). The workshop is divided into three modules: 1) Quality control (QC) of raw sequence data, 2) RNA-seq alignment and quantitation 3) RNA-seq differential expression and visualization.

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
> `##Build hisat2 index##

hisat2-build ../Desktop/sRNA_in_Archaea/Data/RNA-seq1-H2O2/EDGEpro_out_HFX/HFX_genome_edited_STAR.fa hisat2_HFX_genome_index_out/HFX_NCBI`

> `##Align##
hisat2 --verbose --no-spliced-alignment --rna-strandness RF --dta -I 0 -X 500 -x hisat2_HFX_genome_index_out/HFX_NCBI -1 HTSeq_subsample.py_reads/trimmed/C1_mRNA_10percent_subsample_trimmed_1.fq  -2 HTSeq_subsample.py_reads/trimmed/C1_mRNA_10percent_subsample_trimmed_2.fq -S HFX_hisat2_practice_alignments/HFX_C1_mRNA_HTSeq_subsample_hisat2_alignment.sam`

> `##Convert sam to bam##
samtools view -bS /Users/DRG/QBB_Puerto_Rico_2016_testdata/HFX_hisat2_practice_alignments/HFX_C1_mRNA_HTSeq_subsample_hisat2_alignment.sam > /Users/DRG/QBB_Puerto_Rico_2016_testdata/HFX_hisat2_practice_alignments/HFX_C1_mRNA_HTSeq_subsample_hisat2_alignment.bam`


#*Module 3: RNA-seq Differential Expression Analysis and Visualization*

