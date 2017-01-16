# PATHS
SHELL := /bin/bash
DATA := ../analysis/

BBMAP := bbmap.sh
BOWTIE := bowtie2

#================================================================================
# RECEPIES

.PHONY: clean

all: alignments keys merge

# little hackish here, could clean up the hard-coded 'indel.'
alignments: $(addprefix ./test-aligners/indel., bowtie.sam bbmap.sam pyalign.csv)

keys: $(addprefix ./test-aligners/indel.,$(addsuffix .key.txt, bowtie bbmap pyalign))

merge: $(addprefix ./output/indel., $(addsuffix .merge.txt, bowtie bbmap))

clean:
	rm -f indel.txt indel.fasta test-aligners/* output/*

#================================================================================
# PRE-REQ'S

indel.txt: indel-gen.py $(DATA)/ref-seq.fasta
	./$< -f $(word 2, $^) > $@

indel.fasta: indel.txt
	awk -F, '{if (NR > 1) print ">"$$1"\n"$$5}' $< > $@

#================================================================================
# ALIGNMENTS

./test-aligners/%.bbmap.sam: $(DATA)/ref-seq.fasta %.fasta
	@# Run BBMap at its most sensitive levels
	@# t controls number of threads
	@# -XmxNg controls how much RAM to allocate to each thread
	@# k=8 is the smallest kmer (e.g. more sensitive)
	@# vslow=t is a macro that bumps internal setting to most sensitive
	$(BBMAP) in=$(word 2, $^) outm=$@ ref=$< nodisk \
	    vslow=t k=8 -Xmx8g mdtag=t t=4 ordered=t 2> $(@:.sam=.err)

./test-aligners/ref-seq.build.out: $(DATA)/ref-seq.fasta
	@# Must build a bowtie reference before aligning
	$(BOWTIE)-build $< $(@:.build.out=) > $@

./test-aligners/%.bowtie.sam: ./test-aligners/ref-seq.build.out %.fasta
	@# Make ref-seq.build.out a req to ensure it's built before alignment
	@# -N 1 allows for a mismatch in the seeding stage (more sensitive)
	@# -L 4 sets the length of the seeds to the lowest possible setting
	@# -i C,1,0 ensures that there is a seed at every position
	@# --gbar 1 allows there to be gaps in the alignment upto the first pos
	@# --very-sensitive is a macro to max other sensitivity settings
	@# -p 4 tunes the number of procs
	$(BOWTIE) -x $(<:.build.out=) -f -U $(word 2, $^) -S $@ \
	    -p 4 -N 1 -L 4 -i C,1,0 --gbar 1 --very-sensitive \
	    --reorder 2> $(@:.sam=.out)

./test-aligners/%.pyalign.csv: align.py $(DATA)/ref-seq.fasta %.fasta
	./$< -p1 $(word 2, $^) $(word 3, $^) > $@


#================================================================================
# KEYS

./test-aligners/%.pyalign.key.txt: ./test-aligners/%.pyalign.csv
	tail -n+2 $< | tr ',' '\n' > $@

./test-aligners/%.bbmap.key.txt: ./test-aligners/%.bbmap.sam
	awk '{if (NR > 3) print $$1, $$4, $$6, substr($$14, 6)}' $< > $@

./test-aligners/%.bowtie.key.txt: ./test-aligners/%.bowtie.sam
	awk '{if (NR > 3) print $$1, $$4, $$6, substr($$18, 6)}' $< > $@

# merge results with indel key to get final answer
./output/%.bowtie.merge.txt: %.txt ./test-aligners/%.bowtie.key.txt
	tail -n+2 $< | cut -d, -f1-3 | tr ',' ' ' | join -j1 - $(word 2, $^) | \
	    awk '$$2 != $$5' > $@

./output/%.bbmap.merge.txt: %.txt ./test-aligners/%.bbmap.key.txt
	tail -n+2 $< | cut -d, -f1-3 | tr ',' ' ' | join -j1 - $(word 2, $^) | \
	    awk '$$2 != $$5' > $@
