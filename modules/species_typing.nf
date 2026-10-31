process KLEBORATE {
    tag "$sample"
    label 'process_medium'
    conda "${projectDir}/envs/amr.yml"
    publishDir { "${params.outdir}/${sample}/chromosome/kleborate" }, mode: 'copy'

    input:
    tuple val(sample), path(chromosome), val(preset)

    output:
    tuple val(sample), path("${sample}.kleborate.txt"), emit: results

    script:
    // Each preset writes its own <module>_output.txt, so copy whichever one
    // appeared rather than hard-coding the KpSC filename -- that mismatch is
    // why non-pneumoniae Klebsiella produced no result.
    """
    kleborate -p ${preset} \\
              -a ${chromosome} \\
              -o kleborate_out
    out=\$(ls kleborate_out/*_output.txt 2>/dev/null | head -1)
    if [ -n "\${out}" ]; then
        cp "\${out}" ${sample}.kleborate.txt
    else
        echo "# kleborate produced no output for preset ${preset}" > ${sample}.kleborate.txt
    fi
    """
}

process PAREX {
    tag "$sample"
    label 'process_medium'
    conda "${projectDir}/envs/amr.yml"
    publishDir { "${params.outdir}/${sample}/chromosome/parex" }, mode: 'copy'

    input:
    tuple val(sample), path(chromosome)

    output:
    tuple val(sample), path("${sample}.parex.txt"), emit: results

    script:
    // ParEx CLI — adjust to match your install path/command name
    """
    if command -v parex &>/dev/null; then
        parex -i ${chromosome} -o ${sample}.parex.txt
    elif command -v ParEx.py &>/dev/null; then
        ParEx.py -i ${chromosome} -o ${sample}.parex.txt
    else
        echo "# ParEx not found - skipping" > ${sample}.parex.txt
    fi
    """
}
