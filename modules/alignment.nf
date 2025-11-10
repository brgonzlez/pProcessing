/*
 * ALIGNMENT{} process will aligned the data against human reference genome and output mapped and unmapped reads in separate files.
 */


process ALIGNMENT {

	conda "${projectDir}/envs/alignment.yaml"

	input:
	path collapsedReads
	path humanReferenceGenome
	tuple val(missing), val(gap), val(seed)
	val parallel

	output:
	tuple path('*SortedMappedreads.fastq'), path('*SortedUnmappedreads.fastq'), emit: sortedReads
	path 'alignment.log', emit: alignmentLog

	script:
	"""
	#!/bin/bash

	bwa index $humanReferenceGenome/*

	alignment() {
	file=\$1

    		sample=\$(basename "\${file%_adapterRemovalOutput.fastq.gz}.fastq.gz")
       		# Making read groups
    		rg_id="\${sample%.fastq*}"  # sample name as id
    		rg_sm="\${sample%.fastq*}" # sample name again
    		rg_pl="illumina"        # I dont think this is very important for this pipeline so its going to be just illumina because why not
    		rg_lb="lib1"            # group id
    		rg_pu="unit1"           # not sure what Ill put here

		echo -e "\n[\$(date)] Sample: \$sample , Running alignment against reference . . ."
    		bwa aln -l $seed -n $missing -o $gap -t $task.cpus $humanReferenceGenome/*fa "\${file}" > "\${sample%.fastq*}.sai"
		echo -e "\n[\$(date)] Sample: \$sample , Done!"
	
		echo -e "\n[\$(date)] Sample: \$sample , Converting SAI to SAM . . ."
    		bwa samse -r "@RG\\tID:\${rg_id}\\tSM:\${rg_sm}\\tPL:\${rg_pl}\\tLB:\${rg_lb}\\tPU:\${rg_pu}" $humanReferenceGenome/*fa "\${sample%.fastq*}.sai" "\${file}" > "\${sample%.fastq*}.sam"
		echo -e "\n[\$(date)] Sample: \$sample , Done!"


		echo -e "\n[\$(date)] Sample: \$sample , Converting SAM to BAM and sorting . . ."
		samtools view -bS "\${sample%.fastq*}.sam" > "\${sample%.fastq*}.bam"
		echo -e "\n[\$(date)] Sample: \$sample , Done!"


		echo -e "\n[\$(date)] Sample: \$sample , Checking sanity of BAM file . . ."
		samtools quickcheck "\${sample%.fastq*}.bam"
		echo -e "\n[\$(date)] Sample: \$sample , Done!"

		echo -e "\n[\$(date)] Sample: \$sample , Sorting BAM . . ."
		samtools sort -o "\${sample%.fastq*}Sorted.bam" -O bam -@ $task.cpus "\${sample%.fastq*}.bam"
		echo -e "\n[\$(date)] Sample: \$sample , Done!"

		echo -e "\n[\$(date)] Sample: \$sample , Generating BAM index . . ."
		samtools index "\${sample%.fastq*}Sorted.bam"
		echo -e "\n[\$(date)] Sample: \$sample , Done!"

		echo -e "\n[\$(date)] Sample: \$sample , Getting only mapped reads . . . "
		samtools view -b -@ $task.cpus -F 4 "\${sample%.fastq*}Sorted.bam" > "\${sample%.fastq*}SortedMappedreads.bam"
		samtools index "\${sample%.fastq*}SortedMappedreads.bam"
		echo -e "\n[\$(date)] Sample: \$sample , Done!"

		echo -e "\n[\$(date)] Sample: \$sample , Getting only unmapped reads . . . "
		samtools view -b -@ $task.cpus -f 4 "\${sample%.fastq*}Sorted.bam" > "\${sample%.fastq*}SortedUnmappedreads.bam"
		samtools index "\${sample%.fastq*}SortedUnmappedreads.bam"
		echo -e "\n[\$(date)] Sample: \$sample , Done!"

		echo -e "\n[\$(date)] Sample: \$sample , Extracting reads from BAM . . . "
    	samtools fastq -@ $task.cpus "\${sample%.fastq*}SortedMappedreads.bam" > "\${sample%.fastq*}SortedMappedreads.fastq"
    	samtools fastq -@ $task.cpus "\${sample%.fastq*}SortedUnmappedreads.bam" > "\${sample%.fastq*}SortedUnmappedreads.fastq"
		echo -e "\n[\$(date)] Sample: \$sample , Done!"


		echo -e "\n[\$(date)] Sample: \$sample , Removing temporary files . . . "
		echo -e "\n[\$(date)] Sample: \$sample , All done!"

	}

	export -f alignment
	find ./ -name "*_adapterRemovalOutput.fastq.gz" | parallel -j $parallel alignment 

	cat .command.log >> alignment.log
	"""
}
