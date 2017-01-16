#!/bin/bash
# Clean the names of your sequencing data
# Nathan Lubock
# 11/4/2014

# Remove the pesky header
rename 's/020613_//' *.fastq*

# Trim the _001 from the ends
rename 's/_001\./\./' *

# Remove superfluous lane info
rename 's/_L001//' *

# Delete the primers
rename 's/_[A-Z]*_/_/' *

# Distinguish treatments 1 and 2
rename 's/([a-z,A-Z,-])([1-2])_/$1_$2_/' *

# remove abberant hyphens and change name to nonDoped
rename s/-_/_/ *
rename s/Non-doped[A-Z,a-z]*_/nonDoped_/ *
