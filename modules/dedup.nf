/* 
 * DEDUPLICATION{} process will take aligned data and perform read deduplication.
 */

process DEDUPLICATION {

   publishDir "${params.output}/FASTQ",
             mode: 'copy',
             pattern: '*ReadsdeDuplicated.fastq'

   publishDir "${params.output}/LOG",
             mode: 'copy',
             pattern: 'dedup.log'

  conda "${projectDir}/envs/deduplication.yaml"

  input:
  tuple path(mapped), path(unmapped)
  val parallel

  output:
  tuple path('*unMappedReadsdeDuplicated.fastq'), path('*MappedReadsdeDuplicated.fastq'), emit: deduplicated
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

  for sample in *_good_out.fastq; do
    name=\$(basename "\${sample%_good_out.fastq}.fastq")
    mv "\$sample" "\${name}"
  done
  

  cat *.log > dedup.log
  """
}
