getData:
	# don't want to risk clobbering things if the find fails
	# use && \ to ensure everything is kosher
	find /data/archive/ectemp/error\ correction/sequencing/miseq2/02* -type f -name "*.fastq.gz" | xargs -i cp {} . && \
	./cleanNames.sh && \
	./handClean.sh && \
	pigz -d *.fastq.gz

clean:
	rm -f *fastq*
