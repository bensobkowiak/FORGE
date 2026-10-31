// Per-sample process: collects fixed metrics into a JSON.
// AMR genes are handled by COMBINE_AMR (already run) and passed separately
// to WRITE_RESULTS at the end of the pipeline.
process COMPILE_SAMPLE {
    tag "$sample"
    label 'process_low'
    conda "${projectDir}/envs/qc.yml"
    publishDir { "${params.outdir}/${sample}" }, mode: 'copy'

    input:
    tuple val(sample),
          path(read_stats),
          path(kraken_reads_parsed),
          path(kraken_asm_parsed),
          path(genome_eval),
          path(plasmid_summary),
          path(kleborate_file),
          path(parex_file),
          path(depth_file),
          val(maldi_contigs_val),
          val(removed_lengths_val),
          val(query_count_val),
          val(chr_species_val),
          val(maldi_species)

    output:
    path "${sample}.metrics.json", emit: metrics

    script:
    def krasm_arg     = (!kraken_asm_parsed.name.startsWith('NO_FILE')) ? "--kraken_assembly ${kraken_asm_parsed}" : ""
    def geval_arg     = (!genome_eval.name.startsWith('NO_FILE'))        ? "--genome_eval ${genome_eval}"           : ""
    def klebs_arg     = (!kleborate_file.name.startsWith('NO_FILE'))     ? "--kleborate ${kleborate_file}"          : ""
    def parex_arg     = (!parex_file.name.startsWith('NO_FILE'))         ? "--parex ${parex_file}"                  : ""
    def depth_arg     = (!depth_file.name.startsWith('NO_FILE'))         ? "--depth ${depth_file}"                  : ""
    def maldi_arg     = (maldi_species && maldi_species != 'NA')         ? "--maldi_species '${maldi_species}'"     : ""
    def maldi_ctg_arg = (maldi_contigs_val && maldi_contigs_val != 'NA') ? "--maldi_contigs_used '${maldi_contigs_val}'" : ""
    def removed_arg   = (removed_lengths_val && removed_lengths_val != 'NA') ? "--removed_contig_lengths '${removed_lengths_val}'" : ""
    def query_arg     = (query_count_val && query_count_val != 'NA') ? "--num_query_contigs '${query_count_val}'" : ""
    def chrsp_arg     = (chr_species_val && chr_species_val != 'NA') ? "--chr_contig_species '${chr_species_val}'" : ""
    def cap_arg       = params.filtlong_target_bases ? "--coverage_cap_bases ${params.filtlong_target_bases}" : ""
    """
    combine_results.py \\
        --sample          ${sample} \\
        --nanoplot_stats  ${read_stats} \\
        --kraken_reads    ${kraken_reads_parsed} \\
        ${krasm_arg} \\
        --plasmid_summary ${plasmid_summary} \\
        ${geval_arg} \\
        ${klebs_arg} \\
        ${parex_arg} \\
        ${depth_arg} \\
        ${maldi_arg} \\
        ${maldi_ctg_arg} \\
        ${removed_arg} \\
        ${query_arg} \\
        ${chrsp_arg} \\
        ${cap_arg} \\
        > ${sample}.metrics.json
    """
}

// Final process: runs ONCE after all samples complete.
// Builds (or updates) the wide TSV with one column per resistance gene.
process WRITE_RESULTS {
    label 'process_low'
    conda "${projectDir}/envs/qc.yml"

    input:
    path(all_metrics)   // collected list of *.metrics.json
    path(all_amr)       // collected list of *.amr.tsv
    val(results_file)

    output:
    path "results_written.flag"

    script:
    """
    write_results.py \\
        --results_file ${results_file} \\
        --metrics      ${all_metrics} \\
        --amr          ${all_amr}

    touch results_written.flag
    """
}
