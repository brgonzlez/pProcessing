/*
 * PUBLISH{} process will simply send important results to the specified output PATH in compressed format.
 */

process PUBLISH {
 
 input:
 path output
 tuple path(adaptersLog), path(alignmentLog), path(dedupLog)
 tuple path(unmappedFastq), path(mappedFastq)
 


 output:
 stdout
 
 script:
 """
 mkdir -p ${output}/LOG
 mkdir -p ${output}/FASTQ

 mv *fastq.gz ${output}/FASTQ
 mv *log ${output}/LOG
 """
}
