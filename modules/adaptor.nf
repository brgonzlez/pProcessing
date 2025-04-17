/*
 * ADAPTOR_REMOVAL{} process will take either single or paired-end reads and do two tasks: Find adapter sequences and remove them from the fasta.gz files. Format accepted for paired end reads is R{1-2}_001.fastq.gz. May change.
 */


process ADAPTOR_REMOVAL {

  conda "${projectDir}/envs/adaptor.yaml"

	input:
	tuple path(data), path(output)
	path report

	output:
	stdout

	script:
	"""
  	#!/bin/bash

  	######################################
  	# Remove adaptors and collapse reads #
  	######################################
	
  	# first step: make directory structure
  	mkdir -p $output/adapters
  	mkdir -p $output/collapsed
	
  	# second step: Get adapter sequences and use them.
  	#for file in "$data"/*R1_001.fastq.gz; do
	
  	findAndRemovePaired() {
  	file=\$1
	
        	name=$(basename "${file%1.fastq.gz}")
        	AdapterRemoval --identify-adapters --file1 "${file}"  --file2 "$data"/"${name}002.fastq.gz" > ../results/adapters/"$name"_Adapters.txt
          	awk -F':' '/adapter1/ {print $2}' ../results/adapters/"$name"_Adapters.txt > ../results/adapters/"$name"_adapter1.txt
          	awk -F':' '/adapter2/ {print $2}' ../results/adapters/"$name"_Adapters.txt > ../results/adapters/"$name"_adapter2.txt
          	adapter1cat=$(cat "$output"/adapters/"$name"_adapter1.txt | tr -d '[:space:]')
          	adapter2cat=$(cat "$output"/adapters/"$name"_adapter2.txt | tr -d '[:space:]')
	  	AdapterRemoval \
          	--threads 20 \
          	--adapter1 "$adapter1cat" \
          	--adapter2 "$adapter2cat" \
          	--collapse \
          	--minadapteroverlap 1 \
          	--minlength 25 \
          	--minquality 25 \
          	--gzip \
          	--trimns \
          	--trimqualities \
          	--file1 "${file}" \
          	--file2 ../data/"${name}_2.fastq.gz" \
          	--basename  "${name}_adapterRemovalOutput"
  	}
	
  	if [[ $type == PAIRED ]]; then
	  	export -f findAndRemovePaired
	  	find $data -name "*.gb" | parallel -j $task.cpus findAndRemovePaired  
  	else
    	export -f findAndRemoveSingle
    	find $data -name "*.gb" | parallel -j $task.cpus findAndRemoveSingle
  	fi
	
	
  	# third step: cleaning up
	
  	mv *adapterRemovalOutput* ../results/adapters/
  	mv ../results/adapters/*.collapsed.gz ../results/collapsed/
}
