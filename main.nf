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
        println "aPG version 0.1"
        exit 0
}

if (params.version) {
        version()
}

def help() {
      	println "\n\033[1;31mSYNOPSIS\033[0m"

	println "\n\033[1;33mUSAGE\033[0m"
	println "\nnextflow run main.nf --data <PATH> --output <PATH> --index <PATH> [..OPTIONS..]"
	
	println "\n\033[1;33mMANDATORY\033[0m"
	println "  --data <PATH>		Set data file PATH"
	println "  --output <PATH>		Set output directory PATH"
	println "  --index <INT>		Set index file PATH"
	println "  --type <STRING>		Define sequence type: plasmid OR chromosome"


	println "\n\033[1;33mOPTIONS\033[0m"
	println "  --threads <INT>		Set number of threads for PGGB (default: 48)"
	println "  --help			Print help page and exit"
	println "  --version			Print version and exit"

	
	println "\n\033[1;31mDESCRIPTION\033[0m"
	println "\n\033[1;33m--data <PATH>\033[0m"
	println "Please specify the full PATH of your data. Example: /home/user/mydata/data"
	
	println "\n\033[1;33m--output <PATH>\033[0m"
	println "Please specify the full PATH of your output folder. You need to make the folder first before running the program."
	
	println "\n\033[1;33m--index <PATH>\033[0m"
	println "Please specify the full PATH of your index file. Example: /home/user/index/index.tsv"
	println "index.tsv file should contain 2 fields separated by tab. First field should have accession number, second field number of sequences to download (or all)"
	println "\nExample: \n	562	50\n	631	all"
        exit 0
}

if (params.help) {
        help()
}


// Main workflow

workflow {
    // Ensure mandatory parameters are provided
    if (params.output == "") {
        throw new Exception("Output directory must be specified using --output")
    }
    if (!params.index) {
        throw new Exception("Index file PATH must be specified --index")
    }

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
