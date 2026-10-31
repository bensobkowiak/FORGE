process FILTLONG {
    tag "$sample"
    label 'process_low'
    conda "${projectDir}/envs/filtlong.yml"

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}.filtered.fastq.gz"), emit: reads

    script:
    def target = params.filtlong_target_bases ? "--target_bases ${params.filtlong_target_bases}" : ""
    """
    filtlong --min_length ${params.filtlong_min_length} ${target} ${reads} \
        | gzip > ${sample}.filtered.fastq.gz
    """
}

// Long-read assembly with Flye. --nano-hq suits R10.4.1 + Dorado sup basecalls.
// assembly_info.txt carries per-contig circularity, used by PARSE_GENOME_EVAL.
process FLYE {
    tag "$sample"
    label 'process_high'
    // flye_bin (a system install) takes precedence and needs no conda env at all
    conda params.flye_bin ? null : (params.flye_env ?: "${projectDir}/envs/flye.yml")
    publishDir { "${params.outdir}/${sample}/assembly/flye" }, mode: 'copy'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("assembly.fasta"),    emit: assembly
    tuple val(sample), path("assembly_info.txt"), emit: info
    path "*"

    script:
    def gsize = params.flye_genome_size ? "--genome-size ${params.flye_genome_size}" : ''
    // Explicit null check, not ?:, so an intentional 0 is preserved. Defaults
    // keep this working if nextflow.config is older than this module.
    def iters = params.flye_iterations != null ? params.flye_iterations : 1
    def exe   = params.flye_bin ?: 'flye'
    // Genome size in bp, for the runtime coverage test below.
    def gs_bp = 0
    if (params.flye_genome_size) {
        def g = (params.flye_genome_size as String).trim().toLowerCase()
        def mult = ['g': 1000000000L, 'm': 1000000L, 'k': 1000L]
        def suf  = mult.keySet().find { g.endsWith(it) }
        gs_bp = suf ? ((g[0..-2] as BigDecimal) * mult[suf]) as long : (g as BigDecimal) as long
    }
    def acov_n = (gs_bp > 0 && params.flye_asm_coverage) ? (params.flye_asm_coverage as int) : 0
    """
    # --asm-coverage caps the disjointig stage so high-coverage samples don't
    # blow up memory, but Flye aborts outright if the data can't reach the
    # requested depth. Apply it only when this sample actually has enough,
    # with a 20% margin — short-read libraries otherwise fail here.
    ACOV=""
    if [ ${acov_n} -gt 0 ]; then
        BASES=\$(zcat ${reads} | awk 'NR%4==2 {n+=length(\$0)} END {print n+0}')
        NEEDED=\$(( ${gs_bp} * ${acov_n} * 12 / 10 ))
        COV=\$(( BASES / ${gs_bp} ))
        if [ "\$BASES" -gt "\$NEEDED" ]; then
            ACOV="--asm-coverage ${acov_n}"
            echo "Coverage ~\${COV}x — using --asm-coverage ${acov_n}" >&2
        else
            echo "Coverage ~\${COV}x below ${acov_n}x threshold — assembling with all reads" >&2
        fi
    fi

    ${exe} --nano-hq ${reads} \\
         --out-dir flye_out \\
         --threads ${task.cpus} \\
         ${gsize} \$ACOV \\
         --iterations ${iters}

    cp flye_out/assembly.fasta    assembly.fasta
    cp flye_out/assembly_info.txt assembly_info.txt
    cp flye_out/flye.log          flye.log
    """
}

process MOB_RECON {
    tag "$sample"
    label 'process_medium'
    conda "${projectDir}/envs/mob_suite.yml"
    publishDir { "${params.outdir}/${sample}/assembly/mob_recon" }, mode: 'copy'

    input:
    tuple val(sample), path(assembly)

    output:
    tuple val(sample), path("mob_recon_out"), emit: mob_dir
    path "mob_recon_out/*"

    script:
    """
    mob_recon --infile ${assembly} \\
              --outdir mob_recon_out \\
              --num_threads ${task.cpus}
    """
}

process CALC_DEPTH {
    tag "$sample"
    label 'process_medium'
    conda "${projectDir}/envs/mapping.yml"
    publishDir { "${params.outdir}/${sample}/assembly" }, mode: 'copy'

    input:
    tuple val(sample), path(reads), path(chromosome), path(plasmids), path(contig_map)

    output:
    tuple val(sample), path("${sample}.depth.tsv"), emit: depth

    script:
    """
    # Map reads to filtered chromosome + plasmids so Chr_Mean_Depth reflects
    # only MALDI-matched contigs and circularity detection is not confused by
    # non-MALDI chromosome contigs removed by FILTER_GENUS_CONTIGS.
    cat ${chromosome} ${plasmids} > combined_ref.fasta
    minimap2 -ax map-ont --secondary=no -t ${task.cpus} combined_ref.fasta ${reads} \
        | samtools sort -@ ${task.cpus} -o aln.bam
    samtools coverage aln.bam > coverage.txt
    parse_depth.py --coverage coverage.txt --contig_map ${contig_map} --bam aln.bam --sample ${sample} > ${sample}.depth.tsv
    """
}

process SEPARATE_CONTIGS {
    tag "$sample"
    label 'process_low'
    conda "${projectDir}/envs/qc.yml"
    publishDir { "${params.outdir}/${sample}/assembly" }, mode: 'copy'

    input:
    tuple val(sample), path(mob_dir), path(assembly), path(asm_info)

    output:
    tuple val(sample), path("contigs/chromosome.fasta"),              emit: chromosome
    tuple val(sample), path("contigs/rescued_mge_contigs.tsv"),       emit: rescued_mge
    tuple val(sample), path("contigs/${sample}_*.fasta"),             emit: plasmids,           optional: true
    tuple val(sample), path("contigs/plasmids_all.fasta"),            emit: plasmids_combined
    tuple val(sample), path("contigs/contig_to_plasmid.tsv"),         emit: contig_plasmid_map
    tuple val(sample), path("contigs/contig_to_mob_annotation.tsv"),  emit: mob_annotation
    tuple val(sample), path("plasmid_summary.txt"),                    emit: pla_summary

    script:
    def rescue = (params.mge_rescue == false) ? '--no_mge_rescue' : ''
    def maxlen = params.mge_max_len ? "--mge_max_len ${params.mge_max_len}" : ''
    """
    mkdir -p contigs
    process_mob_recon.py ${mob_dir} ${assembly} contigs ${sample} \\
        --assembly_info ${asm_info} ${maxlen} ${rescue} > plasmid_summary.txt
    """
}

// Filter assembled contigs to those classified as the MALDI genus.
// Runs Kraken2 per-contig then filter_genus_contigs.py.
// If the assembly is already single-genus (MALDI) no contigs are dropped.
// If multi-species, non-MALDI contigs are excluded before mob_recon.
process FILTER_GENUS_CONTIGS {
    tag "$sample"
    label 'process_medium'
    conda "${projectDir}/envs/qc.yml"
    publishDir { "${params.outdir}/${sample}/chromosome" }, mode: 'copy'

    input:
    tuple val(sample), path(assembly), val(maldi_genus)
    path(kraken_db)

    output:
    tuple val(sample), path("contigs_filtered.fasta"),      emit: assembly
    tuple val(sample), path("maldi_contigs_used.txt"),      emit: flag
    tuple val(sample), path("removed_contig_lengths.txt"),  emit: removed_lengths
    tuple val(sample), path("query_contigs.fasta"),         emit: query
    tuple val(sample), path("query_contig_info.txt"),       emit: query_info
    tuple val(sample), path("chromosome_contig_species.txt"), emit: chr_species

    script:
    def mmap = params.kraken_memory_mapping ? '--memory-mapping' : ''
    """
    _fallback() {
        cp ${assembly} contigs_filtered.fasta
        echo 'No' > maldi_contigs_used.txt
        echo 'NA' > removed_contig_lengths.txt
        : > query_contigs.fasta
        printf '0\\tNA\\n' > query_contig_info.txt
        echo 'NA' > chromosome_contig_species.txt
    }

    if kraken2 --db ${kraken_db} \\
               --fasta-input \\
               --threads ${task.cpus} \\
               --report contig_kraken.report \\
               --output contig_kraken.out \\
               ${mmap} \\
               ${assembly}; then
        filter_genus_contigs.py \\
            --assembly      ${assembly} \\
            --kraken_out    contig_kraken.out \\
            --kraken_report contig_kraken.report \\
            --maldi_genus   "${maldi_genus}" \\
            --out_fasta     contigs_filtered.fasta \\
            --out_flag      maldi_contigs_used.txt \\
            --out_removed   removed_contig_lengths.txt \\
            --out_query     query_contigs.fasta \\
            --out_query_info query_contig_info.txt \\
            --out_species   chromosome_contig_species.txt \\
            --max_len       ${params.mge_max_len ?: 1000000} \\
            --min_genus_len ${params.maldi_min_genus_len ?: 1000000} \\
            --min_genus_frac ${params.maldi_min_genus_frac ?: 0.25} || _fallback
    else
        _fallback
    fi
    """
}
