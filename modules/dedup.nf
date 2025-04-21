/* 
 * DEDUPLICATION{} process will take aligned data and perform read deduplication.
 */

process DEDUPLICATION {

  conda "${projectDir}/envs/deduplication.yaml"

  input:
  tuple path(mapped), path(unmapped)
  val parallel

  output:
  tuple path('*unmappedDeduplicated.fastq'), path('*mappedDeDuplicated.fastq'), emit: deduplicated
  path 'dedup.log', emit: dedupLog

  script:
  """
  #!/bin/bash

  dedupMapped() {
  file=\$1

  name=\$(basename "\${file}")

  prinseq++ -derep -out_name "\${name}"_unMappedReadsdeDuplicated -fastq "\${file}" >> MappedReadsDeDuplicated.log
  rm *bad_out.fastq
  mv *_good_out.fastq ./"\${name}"_mappedDeDuplicated.fastq
  }

  dedupUnmapped() {
  file=\$1

  name=\$(basename "\${file}")
  prinseq++ -derep -out_name "\${name}"_MappedReadsdeDuplicated -fastq "\${file}" >> unMappedReadsDeDuplicated.log

  rm *bad_out.fastq
  mv *_good_out.fastq ./"\${name}"_unmappedDeduplicated.fastq
  }

  export -f dedupMapped
  export -f dedupUnmapped
  find ./ -name '*SortedUnmappedreads.fastq' | parallel -j $parallel dedupUnmapped
  find ./ -name '*SortedMappedreads.fastq' | parallel -j $parallel dedupMapped

  cat .command.log > dedup.log
  """
}
