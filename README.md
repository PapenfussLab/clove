CLOVE: Classification of genomic fusions into structural variation events
===================

Clove is a program designed to post-process the output of sequencing data structural variation detection algorithms.
Clove scans the output of one or more sets of fusion calls for patterns of complex variations (such as balanced translocations) and summarizes fusion calls as such.
Clove also annotates the variants with read depth information from the alignment file.

----------------
INSTALLING CLOVE
----------------
Copy the latest release from the github page. 
The ...-jar-with-depdendcies.jar should run without any additional requirements

-------------
RUNNING CLOVE
-------------
Invoke CLOVE with java -jar <release.jar>
This prompts you with the follwoing parameters:
	Options (all mandatory -- input can be specified more than once):
	-i <list of breakpoints> <algorithm (Socrates/Delly/Crest/Gustaf/BEDPE)>
	-b <BAM file> 
	-c <mean coverage> <coverage>
	-o <output filename> [default: CLOVE.vcf]
An example run of CLOVE could look like this: 
	java -jar clove-0.11-jar-with-dependencies.jar -i my_results.txt socrates -b my_bam.bam -c 30 7 -o my_calls.vcf
This will take the input in my_results.txt and my_bam.bam to produce calls in my_calls.vcf.

-----------
A FEW NOTES
-----------
1. The input bam file has to be sorted and indexed, as CLOVE is random accessing it.
2. The output VCF file distinguishes calls that have been classified correctly and/or pass the read depth check as "PASS" (or whatever has been provided by the original SV caller) in the filter field. All calls that failed these criteria are indicated with the "FAIL" filter. If you are interested in the "CLOVE approved" calls only, filter for anything that has not "FAIL" in the VCF entry.
3. When constructing the graph of coordinates and fusions, CLOVE discards redundant events (fusions that connect the same two nodes with identical SV type). Therefore, the output vcf is not necessarily complete with respect to the set of inputs. The algorithm will report "Events merged: X" on the command line to indicate if this has happened (X>0).
