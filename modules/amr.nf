// AMR analysis block: AMRfinder + ResFinder + RGI + CefiderocolFinder
// Run on either chromosome or plasmid fasta; 'source' tag drives publishDir subpath.

process AMRFINDER {
    tag "${sample} (${source})"
    label 'process_medium'
    conda "${projectDir}/envs/amr.yml"
    publishDir { "${params.outdir}/${sample}/${pub_dir}/amrfinder" }, mode: 'copy'

    input:
    tuple val(sample), path(fasta)
    val(source)    // 'chromosome' or 'plasmid'
    val(pub_dir)

    output:
    tuple val(sample), path("${sample}_${source}.amrfinder.tsv"), emit: results

    script:
    def db_flag  = params.amrfinder_db       ? "--database ${params.amrfinder_db}"       : ""
    def org_flag = params.amrfinder_organism ? "--organism ${params.amrfinder_organism}" : ""
    """
    if [ ! -s ${fasta} ]; then
        touch ${sample}_${source}.amrfinder.tsv
    else
        amrfinder --nucleotide ${fasta} \\
                  --threads ${task.cpus} \\
                  ${db_flag} \\
                  ${org_flag} \\
                  --output ${sample}_${source}.amrfinder.tsv
    fi
    """
}

process RESFINDER {
    tag "${sample} (${source})"
    label 'process_medium'
    conda "${projectDir}/envs/amr.yml"
    publishDir { "${params.outdir}/${sample}/${pub_dir}/resfinder" }, mode: 'copy'

    input:
    tuple val(sample), path(fasta)
    val(source)
    val(pub_dir)
    path(resfinder_db)

    output:
    tuple val(sample), path("resfinder_${source}"), emit: results

    script:
    """
    if [ ! -s ${fasta} ]; then
        mkdir -p resfinder_${source}
    else
        python -m resfinder \\
            -ifa ${fasta} \\
            -o resfinder_${source} \\
            -s "Other" \\
            --acquired \\
            -db_res ${resfinder_db} \\
            --min_cov 0.6 \\
            --threshold 0.9
    fi
    """
}

process RGI {
    tag "${sample} (${source})"
    label 'process_medium'
    conda "${projectDir}/envs/amr.yml"
    publishDir { "${params.outdir}/${sample}/${pub_dir}/rgi" }, mode: 'copy'

    input:
    tuple val(sample), path(fasta)
    val(source)
    val(pub_dir)

    output:
    tuple val(sample), path("${sample}_${source}.rgi.txt"), emit: results

    script:
    """
    if [ ! -s ${fasta} ]; then
        touch ${sample}_${source}.rgi.txt
    else
        rgi main --input_sequence ${fasta} \\
                 --output_file ${sample}_${source}.rgi \\
                 --input_type CONTIG \\
                 --alignment_tool BLAST \\
                 --threads ${task.cpus} \\
                 --include_loose \\
                 --clean

        # RGI writes .txt suffix
        mv ${sample}_${source}.rgi.txt ${sample}_${source}.rgi.txt 2>/dev/null || \\
        mv ${sample}_${source}.rgi.txt.txt ${sample}_${source}.rgi.txt 2>/dev/null || true
    fi
    """
}

process CEFIDEROCOLFINDER {
    tag "${sample} (${source})"
    label 'process_low'
    conda "${projectDir}/envs/amr.yml"
    publishDir { "${params.outdir}/${sample}/${pub_dir}/cefiderocolfinder" }, mode: 'copy'

    input:
    tuple val(sample), path(fasta)
    val(source)
    val(pub_dir)

    output:
    tuple val(sample), path("${sample}_${source}.cefiderocol.tsv"), emit: results

    script:
    // CefiderocolFinder CLI — adjust command to match your install
    """
    if [ ! -s ${fasta} ]; then
        touch ${sample}_${source}.cefiderocol.tsv
    elif command -v cefiderocolfinder &>/dev/null; then
        cefiderocolfinder -i ${fasta} -o ${sample}_${source}.cefiderocol.tsv
    elif command -v CefiderocolFinder.py &>/dev/null; then
        CefiderocolFinder.py -i ${fasta} -o ${sample}_${source}.cefiderocol.tsv
    else
        echo "# CefiderocolFinder not found - skipping" > ${sample}_${source}.cefiderocol.tsv
    fi
    """
}

process PARSE_AMR {
    tag "${sample} (${source})"
    label 'process_low'
    conda "${projectDir}/envs/amr.yml"

    input:
    tuple val(sample),
          path(amrfinder_out),
          path(resfinder_out),
          path(rgi_out),
          path(cefi_out)
    val(source)

    output:
    tuple val(sample), path("${sample}_${source}.amr_parsed.tsv"), emit: parsed

    script:
    """
    parse_amr.py \\
        --amrfinder   ${amrfinder_out} \\
        --resfinder   ${resfinder_out} \\
        --rgi         ${rgi_out} \\
        --cefiderocol ${cefi_out} \\
        > ${sample}_${source}.amr_parsed.tsv
    """
}

process COMBINE_AMR {
    tag "$sample"
    label 'process_low'
    conda "${projectDir}/envs/amr.yml"
    publishDir { "${params.outdir}/${sample}" }, mode: 'copy'

    input:
    tuple val(sample),
          path(chr_amr_parsed),
          path(pla_amr_parsed),
          path(contig_map),
          path(qry_amr_parsed)

    output:
    tuple val(sample), path("${sample}.amr.tsv"), emit: combined

    script:
    def qry_arg = (!qry_amr_parsed.name.startsWith('NO_FILE')) ? "--qry_amr ${qry_amr_parsed}" : ''
    """
    combine_amr.py \\
        --chr_amr    ${chr_amr_parsed} \\
        --pla_amr    ${pla_amr_parsed} \\
        --contig_map ${contig_map} \\
        ${qry_arg} \\
        > ${sample}.amr.tsv
    """
}
