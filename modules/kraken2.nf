process KRAKEN2 {
    tag "${sample} (${mode})"
    label 'process_high'
    conda "${projectDir}/envs/qc.yml"
    publishDir { "${params.outdir}/${sample}/${pub_dir}/kraken2" }, mode: 'copy'

    input:
    tuple val(sample), path(input_file)
    val(mode)        // 'reads' or 'assembly'
    val(pub_dir)     // 'reads' or 'chromosome'
    path(kraken_db)

    output:
    tuple val(sample), path("${sample}_${mode}.kraken2.report"), emit: report
    tuple val(sample), path("${sample}_${mode}.kraken2.parsed"), emit: parsed

    script:
    def fasta_flag  = mode == 'assembly' ? '--fasta-input' : ''
    def gz_flag     = (mode == 'reads' && input_file.name.endsWith('.gz')) ? '--gzip-compressed' : ''
    def mmap_flag   = params.kraken_memory_mapping ? '--memory-mapping' : ''
    // Assemblies are summarised by contig LENGTH, not contig count, which needs
    // the per-sequence output. For reads the per-read file would be enormous and
    // Kraken2's own read-weighted percentages are already the right unit.
    def kout_file   = mode == 'assembly' ? "${sample}_${mode}.kraken2.out" : '/dev/null'
    def weight_args = mode == 'assembly'
                      ? "--kraken_out ${sample}_${mode}.kraken2.out --assembly ${input_file}"
                      : ''
    """
    kraken2 --db ${kraken_db} \\
            --threads ${task.cpus} \\
            --report ${sample}_${mode}.kraken2.report \\
            --output ${kout_file} \\
            ${fasta_flag} ${gz_flag} ${mmap_flag} \\
            ${input_file}

    parse_kraken.py ${sample}_${mode}.kraken2.report ${weight_args} \\
        > ${sample}_${mode}.kraken2.parsed
    """
}

// Re-run Kraken2 for contaminated samples with per-read output (needed for read extraction).
// Only invoked when contamination is detected; avoids re-running Kraken2 for all samples.
