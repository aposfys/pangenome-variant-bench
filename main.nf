#!/usr/bin/env nextflow

/*
 * pangenome-variant-bench
 *
 * Linear-reference versus pangenome-aware variant calling on one GIAB sample,
 * scoped to a single chromosome so a full run fits on a laptop.
 *
 * The comparison step is deliberately the last one and the only one that produces a
 * number: everything before it must be identical between the two arms, or the comparison
 * measures pipeline differences rather than reference differences.
 */

nextflow.enable.dsl = 2

params.sample        = 'HG002'
params.region        = 'chr20'
params.outdir        = 'results'
params.reference     = null   // linear reference FASTA
params.graph         = null   // pangenome graph (GBZ); must exclude params.sample
params.truth_vcf     = null   // GIAB benchmark VCF
params.confident_bed = null   // GIAB high-confidence regions
params.strata_dir    = null   // stratification BEDs

process ALIGN_LINEAR {
    tag "${sample}:${region}"

    input:
    tuple val(sample), val(region), path(reads), path(reference)

    output:
    tuple val(sample), val(region), path("${sample}.linear.bam"), path("${sample}.linear.bam.bai")

    script:
    """
    echo "milestone 1: bwa-mem2 / minimap2 alignment to the linear reference" >&2
    exit 1
    """
}

process ALIGN_GRAPH {
    tag "${sample}:${region}"

    input:
    tuple val(sample), val(region), path(reads), path(graph)

    output:
    tuple val(sample), val(region), path("${sample}.graph.bam"), path("${sample}.graph.bam.bai")

    script:
    """
    echo "milestone 1: vg giraffe alignment against the pangenome graph" >&2
    exit 1
    """
}

process CALL_VARIANTS {
    tag "${sample}:${arm}"

    input:
    tuple val(sample), val(arm), path(bam), path(bai), path(reference)

    output:
    tuple val(sample), val(arm), path("${sample}.${arm}.vcf.gz")

    script:
    """
    echo "milestone 2: DeepVariant (${arm} arm)" >&2
    exit 1
    """
}

process NORMALISE {
    tag "${sample}:${arm}"

    input:
    tuple val(sample), val(arm), path(vcf), path(reference)

    output:
    tuple val(sample), val(arm), path("${sample}.${arm}.norm.vcf.gz")

    script:
    // Representation differences are not disagreements. Never skipped.
    """
    echo "milestone 2: bcftools norm -m -any -f reference" >&2
    exit 1
    """
}

process COMPARE {
    tag "${sample}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample), path(linear_vcf), path(graph_vcf), path(truth), path(confident)

    output:
    path "${sample}.per_stratum.tsv"
    path "${sample}.findings.json"

    script:
    """
    echo "milestone 3: hap.py stratified comparison, then panbench.report" >&2
    exit 1
    """
}

workflow {
    if (!params.reference || !params.graph || !params.truth_vcf || !params.confident_bed) {
        error "reference, graph, truth_vcf and confident_bed are all required"
    }

    log.info """
    ------------------------------------------------------
      sample : ${params.sample}
      region : ${params.region}
      graph  : ${params.graph}
      NOTE   : the graph must not contain ${params.sample};
               evaluating on an included sample measures memorisation.
    ------------------------------------------------------
    """.stripIndent()
}
