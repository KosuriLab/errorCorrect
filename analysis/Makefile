# PATHS
SHELL := /bin/bash
SCRIPTS := ./scripts
DATA := /data/projects/errorCorrect/rawData

# PROGRAMS
BBDUK := bbduk.sh
BBMERGE := bbmerge.sh
BBMAP := bbmap.sh
BOWTIE := bowtie2
ZIP := pigz

# Some string processing to make dependencies easier
READS := $(sort $(wildcard $(DATA)/*.fastq))
BASE := $(subst _R1,, $(filter-out %_R2, $(notdir $(READS:.fastq=))))

#================================================================================
# Recipies

all: csv secondary stats alt finalFigs.pdf paper/lubock-2017.pdf

csv: $(addprefix pipeline/, $(addsuffix .csv, $(BASE)))

secondary: output/errs-all-samples.csv.gz\
    output/char-counts.txt.gz\
    output/read-counts.txt

stats: $(addprefix pipeline/, $(addsuffix .baseq.txt, $(BASE)))

alt: pipeline/1_nonDoped.bowtie.csv pipeline/1_nonDoped.bbmap.csv

finalFigs.pdf: finalFigs.Rmd \
    $(addprefix pipeline/, 1_nonDoped.bowtie.csv 1_nonDoped.bbmap.csv) \
    $(addprefix output/, errs-all-samples.csv.gz char-counts.txt.gz read-counts.txt)
	R -e "rmarkdown::render('$<')"

paper/lubock-2017.pdf: paper/lubock-2017.tex finalFigs.pdf
	cd paper; latexmk -pdf

clean:
	rm -f output/* pipeline/* finalFigs.pdf ref-seq.fasta.fai
	cd paper; latexmk -C

.PHONY: clean

#================================================================================
# Primary Processing Pipeline

pipeline/%.trim.fastq: adapters.fasta $(sort $(DATA)/%_R1.fastq) $(sort $(DATA)/%_R2.fastq)
	@# Actions:
	@#    use bbduk to trim adapters from the raw reads
	@#    will interleave reads into single fastq for easier tracking
	@# Dependencies:
	@#    fasta file with adapter sequence(s) to trim
	@#    all of the raw sequencing fastq's (NOTE will sort lexicographically)
	@echo Trimming - $(word 2, $^) and $(word 3, $^)
	@$(BBDUK) ref=$< in1=$(word 2, $^) in2=$(word 3, $^) out=$@ \
	    k=21 mink=8 hdist=2 hdist2=1 ktrim=r tpe=t tbo=f -Xmx1g \
	    overwrite=t threads=5 2> $(@:.fastq=.err)

pipeline/%.filter.fastq: pipeline/eColi-phiX.fasta pipeline/%.trim.fastq
	@# Actions:
	@#    filter out any reads mapping to these genomes with bbduk
	@#    get rid of any reads with N's in them
	@# Dependencies:
	@#    combined eColi/PhiX genome (dynamically pulled from ncbi)
	@#    trimmed (and interleaved) reads from previous step
	@echo Filtering - $(word 2, $^)
	@$(BBDUK) ref=$< in=$(word 2, $^) out=$@ stats=$(@:.fastq=.stats.txt)\
	    interleaved=t k=27 hdist=1 ktrim=f maxns=0 -Xmx12g \
	    overwrite=t threads=1 2> $(@:.fastq=.err)

pipeline/%.merge.fastq: pipeline/%.filter.fastq
	@# Actions:
	@#    only merge forward and reverse reads that perfectly overlap
	@# Dependencies:
	@#    filtered reads from last step
	@echo Merging - $<
	@$(BBMERGE) in=$< outm=$@ -Xmx4g\
	    interleaved=t t=5 pfilter=1 overwrite=t 2> $(@:.fastq=.err)

pipeline/%.csv: ref-seq.fasta pipeline/%.merge.fastq
	@# Actions:
	@#    Align merge reads and parse the results for errors
	@# Dependencies:
	@#    reference sequence and merged reads from previous step
	@# Note:
	@#    can add python aligner to dependencies to track any changes
	@#    in that file (useful for testing alignment algos)
	@echo Parsing - $(word 2, $^)
	@python $(SCRIPTS)/parse.py -p10 $^ > $@

#===============================================================================
# Secondary Processing Steps

pipeline/eColi-phiX.fasta:
	@# Actions:
	@#    Grab eColi and PhiX genomes from the ncbi
	curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&amp;id=U00096.3&;rettype=fasta&;retmode=text" > $@
	curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&amp;id=NC_001422.1&amp;rettype=fasta&amp;retmode=text" >> $@


pipeline/%.baseq.txt: pipeline/%.merge.fastq
	@# Actions:
	@#    Generate per-base quality stats for paper
	@# Dependencies:
	@#    Merged reads
	@echo Calculating base qualities for - $<
	@$(BBDUK) in=$< qchist=$@ -Xmx2g t=1 2> /dev/null


output/errs-all-samples.csv.gz: $(addprefix pipeline/, $(addsuffix .csv, $(BASE)))
	@# Actions:
	@#    Combine all of the output csv's into a single flat file
	@#    Useful for downstream R analysis
	@# Dependencies:
	@#    The resulting csv's from make csv. We use this instead a wildcard
	@#    to ensure that things are executed in the correct order
	@echo "Merging all csv's into one"
	@echo "Sample,Name,Pos,Type,Diff" > $(basename $@)
	@$(foreach var, $^,\
	    awk -F, -v file=$(notdir $(basename $(var))) 'NR > 1 {print file FS $$0}' $(var) >> $(basename $@);)
	@$(ZIP) $(basename $@)


output/char-counts.txt.gz: $(addprefix pipeline/, $(addsuffix .merge.fastq, $(BASE)))
	@# Actions:
	@#    Get the length of every read for every sample
	@#    Used in per-base error rate calculations
	@# Dependencies:
	@#    All merged intermediates specified by make csv. Again, we use
	@#    this instead of a wildcard to enure proper execution order
	@echo "Getting length of every read for error rate calcs"
	@echo "Name, Len, Sample" > $(basename $@)
	@$(foreach var, $^,\
	    awk -v samp=$(notdir $(var:.merge.fastq=)) '{if (NR % 4 == 1) print substr($$1, 2)","; else if (NR % 4 == 2) print length($$0)"," FS samp}' $(var) | paste - - >> $(basename $@);)
	@$(ZIP) $(basename $@)

output/read-counts.txt:\
    $(addprefix pipeline/, $(addsuffix .merge.fastq, $(BASE))) \
    $(addprefix pipeline/, $(addsuffix .csv, $(BASE)))
	@# Actions:
	@#    get the number of reads and errs for every treatment
	@#    Note that for our csv's we only care about the existance of an
	@#    error, so we unique the names
	@# Dependencies:
	@#    all merged and parsed reads
	@# Note:
	@#    Requires a relatively new version of (G)AWK in order to
	@#    run length(array)

	@# count number of reads
	@echo "Getting perfect vs all reads"
	@wc -l $(filter %.merge.fastq, $^) |\
	    awk '{gsub(/.*\/|\..*/, "", $$2); print $$1/4, $$2}' > $(@D)/reads.tmp

	@# count number of errors
	@$(foreach var, $(filter %.csv, $^),\
	    tail -n+2 $(var) |\
	    awk -F, -v file=$(notdir $(var:.csv=)) '{!seen[$$1]++} END {print length(seen), file}' >> $(@D)/errs.tmp;)

	@# join together
	@echo "Sample Reads Errs" > $@
	@join -j2 <(sort -k2 $(@D)/reads.tmp) <(sort -k2 $(@D)/errs.tmp) >> $@
	@rm -f $(@D)/reads.tmp $(@D)/errs.tmp

#===============================================================================
# ALIGNER COMPARE

pipeline/1_nonDoped.bbmap.alt.sam: ref-seq.fasta pipeline/1_nonDoped.merge.fastq
	@# Actions:
	@#     Run BBMap at its most sensitive levels
	@#     Generate diff strings and calculate MD tag (for parser)
	@# Options:
	@#     t controls number of threads
	@#     -XmxNg controls how much RAM to allocate to each thread
	@#     k=8 is the smallest kmer (e.g. more sensitive)
	@#     vslow=t is a macro that bumps internal setting to most sensitive
	@#     sam=1.3 ensures compatibility with parsing script
	$(BBMAP) in=$(word 2, $^) outm=stdout.sam ref=$< nodisk sam=1.3\
	    vslow=t k=8 -Xmx8g mdtag=f t=15 ordered=t 2> $(@:.sam=.out) | \
	    samtools calmd -e - $< > $@

pipeline/ref-seq.build.out: ref-seq.fasta
	@# Must build a bowtie reference before aligning
	$(BOWTIE)-build $< $(@:.build.out=) > $@

pipeline/1_nonDoped.bowtie.alt.sam: pipeline/ref-seq.build.out\
    pipeline/1_nonDoped.merge.fastq ref-seq.fasta
	@# Actions:
	@#     Run Bowtie2 at its most sensitive levels
	@#     Generate diff strings for parser
	@# Dependencies:
	@#     ref-seq.build.out ensures it's built before alignment
	@# Options:
	@#     -N 1 allows for a mismatch in the seeding stage (more sensitive)
	@#     -L 4 sets the length of the seeds to the lowest possible setting
	@#     -i C,1,0 ensures that there is a seed at every position
	@#     --gbar 1 allows there to be gaps in the alignment upto the first pos
	@#     --very-sensitive is a macro to max other sensitivity settings
	@#     -p tunes the number of procs
	$(BOWTIE) -x $(<:.build.out=) -U $(word 2, $^) -p 15 -N 1 -L 4 \
	    -i C,1,0 --gbar 1 --very-sensitive --no-unal \
	    --reorder 2> $(@:.sam=.out) | \
	    samtools calmd -e - $(word 3, $^) > $@ 2> $(@:.sam=.err)

pipeline/%.csv: pipeline/%.alt.sam ref-seq.fasta
	@# Actions:
	@#     Run samparse on alternative aligners to compare sensitivities
	@# Options:
	@#     -u distributes errors over equivalent regions
	$(SCRIPTS)/samparse/samparse.py -u $< $(word 2, $^) \
	    $(patsubst .%, %, $(suffix $(<F:.alt.sam=))) > $@

