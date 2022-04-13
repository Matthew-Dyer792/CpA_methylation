#!/usr/bin/env nextflow
/*
========================================================================================
    methylKit output to CpA methylation.bed
========================================================================================
    take in methylKit CHH context text files and produce the trinucleotide context while
    keeping only CpA's using pysam and awk
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE
========================================================================================
*/

// Function to get list of [ meta, [ CHH.txt ] ]
def create_chh_channel(filePairs) {
    def array = []

    def chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
                'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16',
                'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

    chromosomes.each { val ->
        if (!file(filePairs[1][0]).exists()) {
            exit 1, "ERROR: Please check input -> The CHH file does not exist!\n${row.fastq_1}"
        }

        def meta = [id: filePairs[0], chr: "$val"]

        temp_array = [ meta, [ file(filePairs[1][0]) ] ]

        array.add(temp_array)
    }

    return array
}

/*
========================================================================================
    BUILD WORKFLOW CHANNELS
========================================================================================
*/


// possibly start at the sorted bam stage to produce the methylKit output


// Create a channel for input read files
Channel
    .fromFilePairs(params.input, size: 1, checkIfExists: true)
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --single_end on the command line." }
    .map { create_chh_channel(it) }
    .flatten()
    .collate( 2 )
    .set { chh_files }

// chh_files.view()

// channel for index folder
Channel
    .fromPath(params.index, type: 'dir', checkIfExists: true)
    .set { index }

// index.view()

/*
========================================================================================
    PIPELINE STEPS
========================================================================================
*/

// process bamToChh {
//     tag "$meta.id"
//     label 'process_medium'

//     publishDir "${params.outdir}/chh_methylation", mode: 'copy'

//     conda ('bioconda::bioconductor-methylkit=1.20')

//     cpus 4
//     memory '48 GB'
//     executor 'local'

//     script:
//     """
//     bash ${projectDir}/bin/commands_CHH.sh /research/project/shared/benoukraf_lab/large_memory_share/proj/CpA_work/dedupe_sorted/DNMT23-T0-G1-1/DNMT23-T0-G1-1_dup.sort.bam DNMT23-T0-G1-1 ${workDir}
//     """
// }

process splitFiles { 
    tag "$meta.id"
    label 'process_low'

    cpus 2
    memory '4 GB'
    executor 'local'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*${meta.chr}.txt"), emit: individual_chr

    script:
    // Add soft-links to original FastQs for consistent naming in pipeline
    def prefix = "${meta.id}"
    """
    [ ! -f  ${prefix}_CHH.txt ] && ln -s $reads ${prefix}_CHH.txt
    awk 'NR==1 || /^${meta.chr}\$/ {print > "${prefix}_CHH.${meta.chr}.txt"}' ${prefix}_CHH.txt
    """
}

process runPysam { 
    tag "$meta.id"
    label 'process_medium'

    conda ('conda-forge::python=3.9 bioconda::pysam=0.19')

    cpus 4
    memory '8 GB'
    executor 'local'

    input:
    tuple val(meta.id), path(reads)

    output:
    val(meta), path("*.bed"), emit: bed_file

    script:
    // Add soft-links to original FastQs for consistent naming in pipeline
    def prefix = "${meta.id}"
    """
    python ${projectDir}/bin/methylKit_out_to_CpA_bed.py --file ${reads} --outDir ./
    """
}

process catCpAFiles { 
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.outdir}/cpa_context", mode: 'copy'

    cpus 2
    memory '4 GB'
    executor 'local'

    input:
    tuple val(meta), path(bed_files)

    output:
    path("${meta.id}.CpA.bed"), emit: cpa_file

    script:
    // Add soft-links to original FastQs for consistent naming in pipeline
    def prefix = "${meta.id}"
    // def files_to_cat = []
    // files_to_cat.add(header)
    // files_to_cat.add(bed_files)


    // This doesnt work as the files that we then cat still have a header
    def header_file = $bed_files[0]
    """
    head -1 $header_file > header.txt
    find . -type f -name '*.bed' -exec sed '1d' {} +
    cat header.txt $bed_files > ${meta.id}.CpA.bed
    sort -k1,1V -k2,2V ${meta.id}.CpA.bed
    """
    
    // cat $bed_files > headless.bed
    // sort -k1,1V -k2,2V headless.bed
    // cat header.txt headless.bed > ${meta.id}.CpA.bed

    // cat $files_to_cat > ${meta.id}.CpA.bed
}

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

//
// WORKFLOW: Run main bowtie2_align analysis pipeline
//
workflow CpA_ANALYSIS {
    splitFiles (chh_files)
    splitFiles.out.view()
    // runPysam (splitFiles.out)
    // runPysam.out.groupTuple().view()
    // catCpAFiles (runPysam.out.groupTuple())
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
//
workflow {
    CpA_ANALYSIS ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
