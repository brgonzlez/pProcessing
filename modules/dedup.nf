/* 
 * DEDUPLICATION{} process will take aligned data and perform read deduplication. It will publish processed fastq to --output
 */

process DEDUPLICATION {

	conda "${projectDir}/envs/deduplication.yaml"

	input:
	tuple path(mapped), path(unmapped)
	val parallel

	output:
	tuple path('*unMappedReadsdeDuplicated.fastq.gz'), path('*MappedReadsdeDuplicated.fastq.gz'), emit: deduplicated
	path 'dedup.log', emit: dedupLog

	script:
	"""
	#!/bin/bash

	dedupMapped() {
	file=\$1

	name=\$(basename "\${file%.fastq}")
	prinseq++ -derep -out_name "\${name}"_MappedReadsdeDuplicated -fastq "\${file}" >> MappedReadsDeDuplicated.log
	}

	dedupUnmapped() {
	file=\$1

	name=\$(basename "\${file%.fastq}")
	prinseq++ -derep -out_name "\${name}"_unMappedReadsdeDuplicated -fastq "\${file}" >> unMappedReadsDeDuplicated.log

	}

	export -f dedupMapped
	export -f dedupUnmapped
	find ./ -name '*SortedUnmappedreads.fastq' | parallel -j $parallel dedupUnmapped
	find ./ -name '*SortedMappedreads.fastq' | parallel -j $parallel dedupMapped

	rm *bad_out.fastq

	rename_and_compress()
	sample=\$1

		name=\$(basename "\${sample%_good_out.fastq}.fastq")
		mv "\$sample" "\${name}"
		gzip "\${name}"
	}
	export -f rename_and_compress
	find ./ -name "*_good_out.fastq" | parallel -j $task.cpus rename_and_compress

	cat *.log > dedup.log
		
	cp *fastq.gz ${params.output}/FASTQ/"
	cp dedup.log "${params.output}/LOG/"
	"""
}
