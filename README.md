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

Tools/commands for this section:

`wget`

`hisat2`

`less`

`wc -l`

`cat`

`echo`

`python`

`fastqCombinePairedEnd.py`

`samtools`

`stringtie`

`cuffcompare`


We will be doing a reference-based RNA-seq alignment analysis today. Inherent in this, we require a reference genome to do this sort of analysis. Fortunately, the *Haloferax volcanii* genome has been sequenced. Not everyone is so lucky, and because so few genomes are actually sequenced and annotated there are reference-indepenedent approaches as well (ie *de novo* assembly of transcripts) but this is beyond the scope of this workshop.

###*Step 1: Filter out Ribosomal RNA reads*
So we've got our reads trimmed of adapters and quality filtered but there is another filtering step required: removal of rRNA reads. You might ask, 
"Why worry about rRNA?" In all of life as we know it, rRNA is the most transcribed RNA species in the cell. It can attribute for as much as 90% of total RNA in a cell! Moreover, rRNA is typically constituitively expressed at high-levels regardless of growth state (translation is occuring all the time). Because of this skew in rRNA abundance, a typical RNA-seq run will result in ~90% of total reads mapping to rRNA. Now-a-days bioinformaticians do not need worry about rRNA due to awesome techniques in removing rRNA during the library preparation stage (see Illumina Ribo-zero kit), effectively depleting rRNA levels to <1% and saving all your sequence reads to target mRNA and other RNA species! Unfortunately, these techniques are well established for Bacteria and Eukarya but not for Archaea... :( 

A couple of other reasons to remove rRNA, more important to bioinformatic analysis, is that something as highly expressed as rRNA can mask expression of more lowly expressed RNA species. So our first task in our RNA-seq adventure is to remove rRNA from our reads. We will do this by mapping our reads to 5S, 16S, and 23S rRNA. Let's pull some reference files from NCBI:

**Get rRNA fasta files**
> `Open Chrome.app and go to 

We will be doing all our alignments with a very fast aligner for both DNA & RNA called hisat2, written by Daehwan Kim in Dr. Steven Salzberg's lab. This is the succesor to tophat2 and should be incorporated into pipelines that use full alignments. Take a look at the manual for explanation of options: https://ccb.jhu.edu/software/hisat2/manual.shtml

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

Let's take a look at what the reads look like. A terminal command to do this is `less /path/to/your/reads`. The format looks awful! Take a look at the man page for it `man less`. `less -S` allows us to view files in proper format. The reads look fine so far. Let's count how many reads there are in each file. `wc -l` allows us to count each line in a file. Remember, fastq files have four lines for each read so we have to divide the `wc -l` output by 4. Or you can do it all in one shot like so: `cat /path/to/you/reads | echo $((`wc -l`/4))

The paired-end read files have different numbers of reads. This will cause problems when we try to align against the reference genome because the aligner is written to take paired-end reads that have equal numbers of mated pairs. Our next task is to synchronize these reads.

**Synchronize paired-end files (thus removing singlets)**
> `$ python /Scripts-master/fastqCombinePairedEnd.py /path/to/reads_1.fq /path/to/reads_2.fq`

Take a quick look at the output files:

`cat /path/to/reads_R1.fq | echo $((`wc -l`/4))`
`cat /path/to/reads_R2.fq | echo $((`wc -l`/4))`

Perfect, they're synced. 

###*Step 2: Align rRNA-filtered reads against NCBI reference genome*
Now that we have rRNA-clean reads, let us do the real alignment.

**Build hisat2 NCBI refseq index**
> `$ hisat2-build /path/to/HFX_NCBI_reference_genome.fa`


**Align filtered reads**
> `$ hisat2 --verbose  --no-spliced-alignment --rna-strandness RF --dta -I 0 -X 500 -x /path/to/hisat2/HFX_index -1 /path/to/synced/rRNA-removed_read1.fq -2 /path/to/synced/rRNA-removed_read2.fq -S /path/to/rRNA-removed_hisat2_alignment.sam`

You now have your first alignment file, which is in the .sam format. You can learn more about the format here: https://samtools.github.io/hts-specs/SAMv1.pdf

> #Task #2
> You may have noticed that the reads started with a letter "C". This stands for control, meaning no treatment to *Haloferax volcanii*. The letter following "C" (eg C1, C2, C3) corresponds to a biological replicate. It is always important to do biological replicates to insure the change you are seeing is indeed real and not a one-time fluke. The reads starting with "O" are *Haloferax volcanii* that has been treated with ionizing radiation, thus causing oxidative stress (hence the "O"). Your task now is to:

> 1. Remove rRNA from the remaining reads

> 2. Synchronize the reads so that each paired-end pair has equal mates.

> 3. Align the rRNA-removed reads against the reference HFX genome


###*Step 3: Assemble transcripts from aligned reads and quantitate*

Because Illumina reads are so short (<250 bp), the alignment of those reads to a reference genome does not directly tell us about the abundance of whole transcript. This is one of the hardest problems in RNA-seq software development: the assembly of mapped reads into full transcripts. In theory, it is very easy to imagine any overlap between where reads overlap can count as an assembled region and then you grow outwards from there. Paired-end reads make this job even easier by filling in gaps where there are no overlaps, but because of the orientation of the way the reads mapped you can confidently fill in the gap between mated pairs. 

	(R1) ---->    [[[[GAP]]]]     <----(R2)
	
			o----------------------o
	
			  	  Fill in gap

Unfortunately, it is not so straightforward when you take into account 5' & 3' UTRs, operons, and truncated or elongated transcripts. Thus, assembling full transcripts is a guesstimate of the full length transcript and must be verified experimentally (ie northern blot). Our next task is to assemble full length transcripts from the aligned reads and then quantitate how many reads fall within the assembled transcript boundaries.

We will be using one of the best transcript assemblers available called `stringtie` which was developed by in Dr. Steven Salzberg's lab. This is the succesor to cufflinks2. In order use `stringtie` with our alignments we have to convert the files into the correct file types.

**Convert sam to bam**
> `$ samtools view -bS /Users/DRG/Desktop/HFX_rRNA/rRNA_removed_alignments/HFX_O3_mRNA_rRNA_mapped_hisat2_rRNA_removed_alignment.sam > /Users/DRG/Desktop/HFX_rRNA/rRNA_removed_alignments/HFX_O3_mRNA_rRNA_mapped_hisat2_rRNA_removed_alignment.bam`

Samtools is a great tool available for parsing and filtering aligned reads in .sam format. It has a whole host of options available and their manual warrants reading extensively: http://www.htslib.org/doc/samtools.html

Samtools view allows use to convert .sam files into a .bam file, which is a binary format (and thus compressed) of a .sam file. Many tools require .bam format input.

**Sort by coordinate & by name**
Next, stringtie requires the aligned reads to be sorted by the coordinates in which they aligned onto the genome. For example, transcript A mapped to coordinates 1..400 on the Chromomsome while transcript B mapped to coordinates 550..637. We want to sort these reads in order: Transcript A, then Transcript B.

We will use another samtools command: `samtools sort`

To sort by coordinate:

> `$ samtools sort /path/to/alignment.bam -o /path/to/output/alignment.cordsorted.bam`

For another tool to count the reads that fall under assembled transcripts (HTSeq-count) that we will use downstream, we need alignments that are sorted by name instead of coordinates. Let's take care of this now.

To sort by name:
> `$ samtools sort -n /path/to/alignment.bam -o /path/to/output/alignment.namesorted`

Now we can assemble the transcripts based on the alignments and a reference annotation. `stringtie` calculates the coverage (ie the number of reads in a given area) and normalizes these counts into an expression value per assembled transcript termed RPKM/FPKM. This value stands for: *Fragments Per Kilobase of exon per Million reads*. This basically means it is expression of a transcript normalized by total length of the assembled transcript and normalized by the total size of the RNA-seq library, thus making it possible to compare Gene A in Sample 1 to Sample 2 even if Sample 1′s RNA-seq library has 60 million pairs of reads and Sample 2′s library has only 30 million pairs of reads. For further information: http://www.cureffi.org/2013/09/12/counts-vs-fpkms-in-rna-seq/

'stringtie' essentially builds a transcriptome based on where reads have aligned in a tab-delimited format called .gtf. You will find this format or varients (.gff, .gff3) used practically all the time for reference gene annotations. We can also provide a reference transcriptome to help with the assembly calling, which we will do today to continue with our reference-based RNA-seq pipeline. 

We will call `stringtie` with a few options to get out tab-delimited files of coverage `-C` and gene abundances `A`. It is worth noting that stringtie builds whole transcriptomes, meaning it will build transcripts for genes and anything else (ie non-coding RNA). We will ignore anything that is not a gene for today. We will also specify the minimum distance (30 nt in this case) allowed between two reads to be called a single transcript with `-m 30`. For more information on stringtie read the manual: https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual

**Assemble reads into transcripts and quantify abundance**
> `$ stringtie /path/to/cordsorted.bam -G /path/to/reference_genome_annotation.gff -o /path/to/assemble_transcriptome_stringtie_out.gtf -A /path/to/gene_abund.tab -C /path/to/cov_refs.gtf -m 30`

> #Task #3
> All of these tasks must be done for each individual alignment in order to compare the expression between samples. Your task now is to:

> 1. Convert alignments into file types usable in downstream transcript assembly and quantitation applications.

> 2. Assemble transcriptome for each sample.

> 3. Find the highest expressed gene in each sample. Hint: cut, sort

> Answer thought question: what does an RPKM/FPKM value of 0 mean? How do you asses this?


Now that we have the individual transcriptomes for each sample, we have to compare each transcriptome against one another to get a consensus transcriptome. We want to make sure we count reads that fall under transcripts that are present in all of the transcriptomes to make sure we don't miss anything. To do this we will use a tool called `cuffcompare`, which is also packaged as `gffcompare`.

First we will need to make an accession file that contains the paths to each stringtie-built transcriptome. You can open a text editor like so: `nano`

In the text editor you want to have a line for each transcriptome. Your file should look similar to this when done:
> `./Desktop/HFX_rRNA/rRNA_removed_alignments/stringtie_out/HFX_C1_rRNA_removed_stringtie_out/HFX_C1_rRNA_removed_transcriptome_stringout.gtf
./Desktop/HFX_rRNA/rRNA_removed_alignments/stringtie_out/HFX_C2_rRNA_removed_stringtie_out/HFX_C2_rRNA_removed_transcriptome_stringout.gtf
./Desktop/HFX_rRNA/rRNA_removed_alignments/stringtie_out/HFX_C3_rRNA_removed_stringtie_out/HFX_C3_rRNA_removed_transcriptome_stringout.gtf
./Desktop/HFX_rRNA/rRNA_removed_alignments/stringtie_out/HFX_O1_rRNA_removed_stringtie_out/HFX_O1_rRNA_removed_transcriptome_stringout.gtf
./Desktop/HFX_rRNA/rRNA_removed_alignments/stringtie_out/HFX_O2_rRNA_removed_stringtie_out/HFX_O2_rRNA_removed_transcriptome_stringout.gtf
./Desktop/HFX_rRNA/rRNA_removed_alignments/stringtie_out/HFX_O3_rRNA_removed_stringtie_out/HFX_O3_rRNA_removed_transcriptome_stringout.gtf`

To close the editor and save the file press `ctrl + x` --> press `y` to save --> name `gtf_acession_list.txt`

Next to get a consensus trancriptome:

**Merge transcriptomes**
> `$ cuffcompare -o /path/to/conensus_transcriptome_cuffcompare_output.gtf -i /path/to/gtf_acession_list.txt -r /path/to/reference_genome_gene_annotation.gtf`

The output of this is a consensus transcriptome called `.combined.gtf` With the `r` option we compared out individual transcriptomes against the reference gene annotation, which gave a category code for each assembled transcript based on the relationship it has with the reference gene annotation. For example, if an assembled transcript STRG.1 matches a reference annotation it will be given a category code of "=". 

A cheatsheat for category markers:
> Class codes:

> "u" is novel intergenic trancscript, (you may check that whether some of them could be expressing retro-viral elements).

> "i" is intronic transcript (there are some intronic lncRNA and repetitive elements get expressed in cell-type specific manner)

> "j"  is potential novel isoform (alternately spliced, since has no splice site match, you might check the coding potential)

> "x" is cis-antisense transcript (exonic overlap but opposite strand)

> "o" is "other overlap" - that is an exonic overlap with the reference transcript that doesn't fall in any other, "more interesting" overlap categories - e.g. no splice sites match ('j' class), no containment ('c' code) etc. These 'o' codes could be assigned, for example, to assembled single-exon fragments that happen to overlap one of the terminal exons of a reference transcript (but not enough to make it "contained")

> "p" is "polymerase run" - it's supposed to signal that the relative positioning of the transcripts to the reference transcript suggests a potential polymerase read-through downstream of the 3' UTRs of the reference.  ( I would ignore it)

> "=" is complete match.

> "e" is single exon transfag with some basepairs of intron retention ( could be premRNA contaminant)

Since we are focused on reference gene expression, we will only work with assembled transcripts that have the category code "=". To pull out these transcripts from everything else we will use the command `awk`. `awk` is extremely useful, but a full explanation and demonstration of how it works is far beyond the scope of this workshop so please just have faith and execute the command:

awk '$22 ~ /=/ { print }' /path/to/consensus_transcriptome_combined.gtf > /path/to/consensus/transcriptome/genes_only.gtf

#*Module 3: RNA-seq Differential Expression Analysis and Visualization*





