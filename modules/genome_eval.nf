process CHECKM2 {
    tag "$sample"
    label 'process_high'
    conda params.checkm2_env ?: "${projectDir}/envs/checkm2.yml"
    publishDir { "${params.outdir}/${sample}/${pub_dir}/checkm2" }, mode: 'copy'

    input:
    tuple val(sample), path(assembly), val(pub_dir)

    output:
    tuple val(sample), path("checkm2_output/quality_report.tsv"), emit: report
    path "checkm2_output/*"

    script:
    def db_flag = params.checkm2_db ? "--database_path ${params.checkm2_db}" : ""
    def checkm2_cmd = System.getProperty('os.name').toLowerCase().contains('mac') ? 'checkm2_fork' : 'checkm2'
    """
    ${checkm2_cmd} predict --input ${assembly} \\
                           --output-directory checkm2_output \\
                           --threads ${task.cpus} \\
                           --lowmem \\
                           --force \\
                           ${db_flag}
    """
}

process BUSCO {
    tag "$sample"
    label 'process_high'
    conda "${projectDir}/envs/busco.yml"
    publishDir { "${params.outdir}/${sample}/chromosome/busco" }, mode: 'copy'

    input:
    tuple val(sample), path(chromosome)

    output:
    tuple val(sample), path("busco_out/short_summary*.txt"), emit: summary
    path "busco_out/*"

    script:
    """
    busco -i ${chromosome} \\
          -o busco_out \\
          -m genome \\
          -l ${params.busco_lineage} \\
          --cpu ${task.cpus} \\
          -f
    """
}

process MLST {
    tag "$sample"
    label 'process_low'
    conda "${projectDir}/envs/mlst.yml"
    publishDir { "${params.outdir}/${sample}/chromosome/mlst" }, mode: 'copy'

    input:
    tuple val(sample), path(chromosome)

    output:
    tuple val(sample), path("${sample}.mlst.tsv"), emit: results

    script:
    def scheme_flag = params.mlst_scheme ? "--scheme ${params.mlst_scheme}" : ""
    """
    mlst ${scheme_flag} ${chromosome} > ${sample}.mlst.tsv
    """
}

process PARSE_GENOME_EVAL {
    tag "$sample"
    label 'process_low'
    conda "${projectDir}/envs/qc.yml"
    publishDir { "${params.outdir}/${sample}/chromosome" }, mode: 'copy'

    input:
    tuple val(sample),
          path(checkm2_chr_report,  stageAs: 'checkm2_chr.tsv'),
          path(checkm2_full_report, stageAs: 'checkm2_full.tsv'),
          path(busco_summary),
          path(mlst_results),
          path(assembly),
          path(assembly_info)

    output:
    tuple val(sample), path("${sample}.genome_eval.tsv")

    script:
    """
    parse_genome_eval.py \\
        --checkm2       checkm2_chr.tsv \\
        --checkm2_full  checkm2_full.tsv \\
        --busco         ${busco_summary} \\
        --mlst          ${mlst_results} \\
        --assembly      ${assembly} \\
        --assembly_info ${assembly_info} \\
        > ${sample}.genome_eval.tsv
    """
}
