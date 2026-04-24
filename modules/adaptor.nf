/*
 * ADAPTOR_REMOVAL{} process will take either single or paired-end reads and do two tasks: Find adapter sequences and remove them from the fasta.gz files. Format accepted for paired end reads 
 * is *_1.fastq.gz and *_2.fastq.gz.
 */


process ADAPTOR_REMOVAL {

	conda "${projectDir}/envs/adaptor.yaml"

	input:
	tuple path(data), val(type), val(parallel)
	tuple val(min_len), val(min_qual)

	output:
	path 'collapsed/*.gz', emit: collapsedReads
	path 'adapters/*_Adapters.txt', emit: adapters
	path 'adapters.log', emit: adaptersLog

	script:
	"""
  	#!/bin/bash

	mkdir -p "${params.output}/LOG
	mkdir -p "${params.output}/FASTQ
	# dont crash with empty wildcards
	shopt -s nullglob

  	# first step: make directory structure
  	mkdir -p adapters
  	mkdir -p collapsed



  	# second step: Get adapter sequences and use them.

	######################
	# DEFINING FUNCTIONS #
	######################
	
  	findAndRemovePaired() {
  	file=\$1
	
        	name=\$(basename "\${file%_1.fastq.gz}")
        	AdapterRemoval --identify-adapters --file1 "\${file}"  --file2 $data/"\${name}_2.fastq.gz" > adapters/"\$name"_Adapters.txt
          	awk -F':' '/adapter1/ {print \$2}' adapters/"\${name}"_Adapters.txt > adapters/"\${name}"_adapter1.txt
          	awk -F':' '/adapter2/ {print \$2}' adapters/"\${name}"_Adapters.txt > adapters/"\${name}"_adapter2.txt
          	adapter1cat=\$(cat adapters/"\${name}"_adapter1.txt | tr -d '[:space:]')
          	adapter2cat=\$(cat adapters/"\${name}"_adapter2.txt | tr -d '[:space:]')

		echo -e "For sample \${name} -> \nAdapter 1: \${adapter1cat} \nAdapter 2: \${adapter2cat}\n"

	  	AdapterRemoval \
          	--threads $task.cpus \
          	--adapter1 "\$adapter1cat" \
          	--adapter2 "\$adapter2cat" \
          	--collapse \
          	--minadapteroverlap 1 \
          	--minlength $min_len \
          	--minquality $min_qual \
          	--gzip \
          	--trimns \
          	--trimqualities \
          	--file1 "\${file}" \
          	--file2 $data/"\${name}_2.fastq.gz" \
          	--basename  "\${name}_adapterRemovalOutput"

		
		mv *"\${name}_adapterRemovalOutput"* adapters/
		mv adapters/*.collapsed.gz collapsed/"\${name}_adapterRemovalOutput.fastq.gz"
  	}

	findAndRemoveSingle() {
	file=\$1
		name=\$(basename "\${file%.fastq.gz}")
		AdapterRemoval --identify-adapters --file1 "\${file}" > adapters/"\${name}"_Adapters.txt
		awk -F':' '/adapter1/ {print \$2}' adapters/"\${name}"_Adapters.txt > adapters/"\${name}"_adapter1.txt
          	awk -F':' '/adapter2/ {print \$2}' adapters/"\${name}"_Adapters.txt > adapters/"\${name}"_adapter2.txt

		adapter1cat=\$(cat adapters/"\${name}"_adapter1.txt | tr -d '[:space:]')
          	adapter2cat=\$(cat adapters/"\${name}"_adapter2.txt | tr -d '[:space:]')

		echo -e "For sample \${name} -> \nAdapter 1: \${adapter1cat} \nAdapter 2: \${adapter2cat}\n"

		AdapterRemoval \
          	--adapter1 "\$adapter1cat" \
          	--adapter2 "\$adapter2cat" \
        	--minadapteroverlap 1 \
        	--minlength $min_len \
        	--minquality $min_qual \
        	--gzip \
		--threads $task.cpus \
        	--trimns \
        	--trimqualities \
        	--file1 "\${file}" \
        	--basename "\${name}_adapterRemovalOutput"

		mv *"\${name}_adapterRemovalOutput"* adapters/
		mv adapters/*.truncated.gz collapsed/"\${name}_adapterRemovalOutput.fastq.gz"
	}



  	if [[ $type == PAIRED ]]; then
	  	export -f findAndRemovePaired
	  	find $data/* -name "*_1.fastq.gz" | parallel -j $parallel findAndRemovePaired  
  	else
    		export -f findAndRemoveSingle
    		find $data/* -name "*.fastq.gz" | parallel -j $parallel findAndRemoveSingle
  	fi

	cat .command.log > adapters.log
	
	cp adapters.log "${params.output}/LOG/"
	"""
}
