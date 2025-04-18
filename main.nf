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
include { PUBLISH } from './modules/publish.nf'


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
    if (!params.output) {
        throw new Exception("Output directory must be specified using --output")
    }
    if (!params.data) {
        throw new Exception("Data directory must be specified using --data")
    }
    if (!params.type) {
        throw new Exception("Read type must be specified using --type")
    }
    if (!params.ref) {
        throw new Exception("Indexed human reference genome must be specified using --ref")
    }
    // Running check

    println "\n\033[1;33mCHECKING PARAMETERS\033[0m"
    println "\n======================="
    println "\n\033[1;31mMANDATORY PARAMETERS\033[0m"
    println "\n\033[1;37mData\033[0m: ${params.data}"
    println "\n\033[1;37mOutput\033[0m: ${params.output}"
    println "\n\033[1;37mType\033[0m: ${params.type}"
    println "\n======================="
    println "\n\033[1;31mADAPTOR REMOVAL\033[0m"
    println "\n\033[1;37mMin. Length\033[0m: ${params.MIN_LENGTH}"
    println "\n\033[1;37mMin. Quality\033[0m: ${params.MIN_QUALITY}"
    println "\n======================="
    println "\n\033[1;31mALIGNMENT\033[0m"
    println "\n\033[1;37mMissing Prob.\033[0m: ${params. MISSING_PROB}"
    println "\n\033[1;37mGap Fraction\033[0m: ${params.GAP_FRACTION}"
    println "\n\033[1;37mSeed\033[0m: ${params.SEED}"
    println "\n=======================\n\n"



println '''
                          .
                          |~~
                          |~~
                         /L\\
                  ,.---./LLL\\. _.--.
                .\'(|~~`/LLLLL\\` )  )`.,
              .(\\`  |~~/LLLLLLL\\   )  )-. .
             (  ( /L\\/LLLLLLLLL\\ )  )   `|~~
            ((\\`_./LLL\\LLLLLLLLLL\\`.)_),.)|~~
                /LLLLL\\.=.=.=.=|        /L\\
                 |.=.| .-._.-. |       /LLL\\   ~\'~
                 |  [| | | | | |      /LLLLL\\
                 |   | | | | | | _   _|] _=.|
        ~\'~      |  [| |_|_|_| || |_| |_| | |
                 |  |~~        |=.=.=.=.=.| |       .
                 |  |~~        |    |~~   | |       |~~
                 | /L\\ .-._.-. |    |~~   | |       |~~
                 |/LLL\\| | | | |   /L\\    |/       /L\\
                 |].=.|_ | _ | _  /LLL\\   |       /LLL\\
           ,- _--|]] [| |_| |_| |/LLLLL\\  |      /LLLLL\\
          (|_| |_|]---|.=.=.=.=./LLLLLLL\\ _   _ /LLLLLLL\\
           \\\\.=.=.=|\\\\_/           |.=.=.|_| |_| |_|.=.=.|
           /|[]   |              | []  |.=.=.=.=.|  [] |
           ||     |    .-._.-.   |     | .-----. |     |
           \\\\|     |    | | | |   |     |/|||||||\\\\|     |
            |  [] |    | | | |   |     ]|||||||||[     |
            |  __ |    |_|_|_|   |  [] ]|||| ||||[ []  |
            | /<\\\\_\\    ____      |     ]|||| ||||[     |
            |/ |  "\\\\__/  ) \\\\.-.  |     ]|||||||||[     |_
           /"  )\\\\_ >  ) >\\\\__ ")\\\`\\\\_     ]|||||||||[ ,_./\\`.\\\\
        __/ _/ _ ,| \\\\  __  "|_  ) |_   ]|||||||||[/(\\"_ -">\\\\_
       /> )"__/ \\\\___  "  \\\\__  _) \\\\_ -\\\\_.==___===/.<  \\\\__(\\\\_ \\\\
      /  __/ )___   > \\\\_ ) \\\\  \\\\_ "  ).==_____==( <."/ (_<  \\\\)|
     lc_/>.=__.._\\\\"__\\\\_  >_)___\\\\-_/.=________=/___/.__>__"(__/
'''

    // Running the workflow
    ADAPTOR_REMOVAL(tuple(params.data, params.type, params.parallel), tuple(params.MIN_LENGTH, params.MIN_QUALITY))
    ALIGNMENT(ADAPTOR_REMOVAL.out.collapsedReads , params.ref, tuple(params.MISSING_PROB , params.GAP_FRACTION, params.SEED) , params.parallel)
    DEDUPLICATION(ALIGNMENT.out.sortedReads.map { sortedmapped, sortedunmapped -> tuple(mapped , unmapped) }, params.parallel)
    PUBLISH(params.output, tuple(ADAPTOR_REMOVAL.out.adaptersLog, ALIGNMENT.out.alignmentLog, DEDUPLICATION.out.dedupLog), 
	    DEDUPLICATION.out.deduplicated.map { unmapped, mapped -> tuple(unmapped, mapped)})
}
