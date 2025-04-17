/*
 * ALIGNMENT{} process will aligned the data against human reference genome and output mapped and aligned reads in separate files.
 */


process ALIGNMENT {
	conda "${projectDir}/envs/alignment.yaml"

	input:
	path collapsedReads
	path humanReferenceGenome
	tuple val(THREADS), val(MISSING_PROB), val(GAP_FRACTION), val(SEED), val(MIN_READ_LENGTH), val(MAX_READ_LENGTH)

	output:
	stdout
	

	script:
	"""
	for sample in $reads/*; do
    		bwa index $panRef
    		name=\$(basename "\$sample")
    		softClip=\$(grep "\$name" $configFile | awk '{print \$2}')
    
    		# Making read groups
    		rg_id="\${name%.fastq*}"  # sample name as id
    		rg_sm="\${name%.fastq*}" # sample name again
    		rg_pl="illumina"        # I dont think this is very important for this pipeline so its going to be just illumina because why not
    		rg_lb="lib1"            # group id
    		rg_pu="unit1"           # not sure what Ill put here

    		bwa aln -l $seedAlignment -n $missingProb -o $gapFraction -t $threadsGlobal $panRef "\$sample" > "\${name%.fastq*}.sai"
    		bwa samse -r "@RG\\tID:\$rg_id\\tSM:\$rg_sm\\tPL:\$rg_pl\\tLB:\$rg_lb\\tPU:\$rg_pu" \
		$panRef "\${name%.fastq*}.sai" "\$sample" > "\${name%.fastq*}.sam"
    
    		samtools view -bS "\${name%.fastq*}.sam" > "\${name%.fastq*}.bam"
    		samtools quickcheck "\${name%.fastq*}.bam"
    		samtools sort -o "\${name%.fastq*}_sorted.bam" -O bam -@ $threadsGlobal "\${name%.fastq*}.bam"
    		samtools index "\${name%.fastq*}_sorted.bam"
    		rm "\${name%.fastq*}.bam"
    		samtools view -b -@ $threadsGlobal -F 4 "\${name%.fastq*}_sorted.bam" > "\${name%.fastq*}_sorted_mappedreads.bam"
    		samtools index "\${name%.fastq*}_sorted_mappedreads.bam"
    		bam trimBam "\${name%.fastq*}_sorted_mappedreads.bam" "\${name%.fastq*}_softclipped.bam" -L "\$softClip" -R "\$softClip" --clip
    		samtools view -q 30 -o "\${name%.fastq*}_qc.bam" "\${name%.fastq*}_softclipped.bam"
    		samtools view -e 'length(seq)>$minReadLength && length(seq)<$maxReadLength' -O BAM -o "\${name%.fastq*}_lg.bam" "\${name%.fastq*}_qc.bam"
    		samtools sort -o "\${name%.fastq*}_DMC_P.bam" -O bam -@ $threadsGlobal "\${name%.fastq*}_lg.bam"
    		samtools coverage "\${name%.fastq*}_DMC_P.bam" > "\${name}_genomicsMetrics.txt"
    		samtools fastq -@ $threadsGlobal "\${name%.fastq*}_DMC_P.bam" > "\${name%.fastq*}_final.fastq"
	done

	rm *sam *sai *_lg.bam *_qc.bam *_sorted_mappedreads.bam*
	cat .command.out >> alignment.log
	"""
}
