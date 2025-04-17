/* 
 * DEDUPLICATION{} process will take aligned data and perform read deduplication.
 */

process DEDUPLICATION {

  conda "${projectDir}/envs/deduplication.yaml"

  input:

  output:

  script:
  """
  #!/bin/bash

  prinseq++ -derep -out_name unMappedReadsdeDuplicated \
          -fastq /home/shared/GenMed/LabPractical2025/exampleRun/results/alignment/nonHuman/unMappedReadsTest.fastq \
          > /home/shared/GenMed/LabPractical2025/exampleRun/results/alignment/nonHuman/deDuplicated/unMappedReadsDeDuplicated.log

  # We need to manually send the output file to the results folder that you can specify
  mv *_good_out.fastq /home/shared/GenMed/LabPractical2025/exampleRun/results/alignment/nonHuman/deDuplicated/unmappedReadsdeDuplicated.fastq
  rm *bad_out.fastq
  """
}
