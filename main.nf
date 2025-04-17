#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is the main document. To modify CPU usage please go to nextflow.config
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// Enable DSL2
nextflow.enable.dsl=2


// Calling modules

include { ADAPTOR_REMOVAL } from './modules/adaptor.nf'
include { ALIGNMENT } from './modules/alignment.nf'
include { DEDUPLICATION } from './modules/dedup.nf'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Print pipeline metadata: Version and Help
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


def version() {
        println "preProcessing version 1"
        exit 0
}

if (params.version) {
        version()
}

def help() {
      	println "\n\033[1;31mSYNOPSIS\033[0m"

	println "\n\033[1;33mUSAGE\033[0m"
	println "\nnextflow run main.nf --data <PATH> --output <PATH> --type <STRING> [..OPTIONS..]"
	
	println "\n\033[1;33mMANDATORY\033[0m"
	println "  --data <PATH>		Set data file PATH"
	println "  --output <PATH>		Set output directory PATH"
	println "  --type <STRING>		Define reads type: SINGLE or PAIRED"


	println "\n\033[1;33mOPTIONS\033[0m"
	println "  --help			Print help page and exit"
	println "  --version			Print version and exit"

	
	println "\n\033[1;31mDESCRIPTION\033[0m"
	println "\n\033[1;33m--data <PATH>\033[0m"
	println "Please specify the full PATH of your data. Example: /home/user/mydata/data/"
	
	println "\n\033[1;33m--output <PATH>\033[0m"
	println "Please specify the full PATH of your output folder. You need to make the folder first before running the program."

        exit 0
}

if (params.help) {
        help()
}


// Main workflow

workflow {

    // Ensure mandatory parameters are provided
    if !(params.output) {
        throw new Exception("Output directory must be specified using --output")
    }
    if !(params.data) {
        throw new Exception("Data directory must be specified using --data")
    }
    if !(params.type) {
        throw new Exception("Read type must be specified using --type")
    }

    // Running check
    
    println "\n\033[1;33mCHECKING PARAMETERS\033[0m"
    println "\n\033[1;Data\033[0m: "
    println "\n\033[1;33mOutput\033[0m: "
    println "\n\033[1;33mType\033[0m: "
    println "\n\033[1;31mADAPTOR REMOVAL\033[0m"
    println "\n\033[1;33mMin. Length\033[0m: (params.MIN_LENGTH)"
    println "\n\033[1;33mMin. Quality\033[0m: (params.MIN_QUALITY)"
    println "\n\033[1;31mALIGNMENT\033[0m"
    println "\n\033[1;31mDEDUPLICATION\033[0m"


    // Running the workflow
    GET_DATA(params.index, params.type)

    PARSE_GENBANK(GET_DATA.out.downloadedFiles.map { fasta, gb, acc, meta -> tuple(fasta, gb) })

    REMOVE_REDUNDANCY(GET_DATA.out.downloadedFiles.map { fasta, gb, acc, meta -> fasta}, PARSE_GENBANK.out.GenBankQCReport)

    FILTER_REDUNDANT(REMOVE_REDUNDANCY.out.clusteredFasta, GET_DATA.out.downloadedFiles.map { fasta, gb, acc, meta -> tuple(fasta, gb, meta)})

    CLUSTERED_INDEX(FILTER_REDUNDANT.out.clusteredSeqsNonRedundant)

    BUILD_GRAPH(FILTER_REDUNDANT.out.clusteredSeqsNonRedundant, CLUSTERED_INDEX.out.indexedClustered, tuple(params.identity , 
                params.segmentLength , params.nMappings , params.threads , params.poaLength , params.poaParams , params.minMatchLen))

    PATH_DISTANCE(FILTER_REDUNDANT.out.cleanFastaGenBank.map { fasta, gb, metadata -> metadata } , BUILD_GRAPH.out.pathDistance)

    GENE_CLUSTERING(FILTER_REDUNDANT.out.fastaDatabase, params.geneClusterThreshold)

    MAKE_GFF(GENE_CLUSTERING.out.clusteredDatabase, FILTER_REDUNDANT.out.cleanFastaGenBank.map { fasta, gb, metadata -> tuple(fasta, gb)})
}
