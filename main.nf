#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { CONCAT_READS; NANOPLOT; PARSE_READ_STATS }    from './modules/read_processing'
include { KRAKEN2 as KRAKEN2_READS; KRAKEN2 as KRAKEN2_ASSEMBLY }   from './modules/kraken2'
include { FILTLONG; FLYE; MOB_RECON; SEPARATE_CONTIGS; CALC_DEPTH;
          FILTER_GENUS_CONTIGS }                                      from './modules/assembly'
include { CHECKM2 as CHECKM2_ALL; CHECKM2 as CHECKM2_CHR;
          BUSCO; MLST; PARSE_GENOME_EVAL }               from './modules/genome_eval'
include { AMRFINDER as AMRFINDER_CHR }                  from './modules/amr'
include { AMRFINDER as AMRFINDER_PLA }                  from './modules/amr'
include { RESFINDER as RESFINDER_CHR }                  from './modules/amr'
include { RESFINDER as RESFINDER_PLA }                  from './modules/amr'
include { RGI as RGI_CHR }                              from './modules/amr'
include { RGI as RGI_PLA }                              from './modules/amr'
include { CEFIDEROCOLFINDER as CEFI_CHR }               from './modules/amr'
include { CEFIDEROCOLFINDER as CEFI_PLA }               from './modules/amr'
include { PARSE_AMR as PARSE_AMR_CHR }                  from './modules/amr'
include { PARSE_AMR as PARSE_AMR_PLA }                  from './modules/amr'
include { AMRFINDER as AMRFINDER_QRY }                  from './modules/amr'
include { RESFINDER as RESFINDER_QRY }                  from './modules/amr'
include { RGI as RGI_QRY }                              from './modules/amr'
include { CEFIDEROCOLFINDER as CEFI_QRY }               from './modules/amr'
include { PARSE_AMR as PARSE_AMR_QRY }                  from './modules/amr'
include { COMBINE_AMR }                                 from './modules/amr'
include { KLEBORATE }                                   from './modules/species_typing'
include { PAREX }                                       from './modules/species_typing'
include { COMPILE_SAMPLE; WRITE_RESULTS }               from './modules/compile_results'

def get_processed_samples(String results_path) {
    def processed = [] as Set
    def f = file(results_path)
    if (f.exists()) {
        f.readLines().drop(1).each { line ->
            def cols = line.tokenize('\t')
            if (cols) processed << cols[0].trim()
        }
    }
    return processed
}

// True when FILTER_GENUS_CONTIGS has produced ALL its outputs. query_contig_info.txt
// is required too, so samples filtered before query contigs existed are re-filtered.
def filter_done(String outdir, String sample) {
    return file("${outdir}/${sample}/chromosome/maldi_contigs_used.txt").exists() &&
           file("${outdir}/${sample}/chromosome/query_contig_info.txt").exists() &&
           file("${outdir}/${sample}/chromosome/chromosome_contig_species.txt").exists()
}

workflow {

    if (!params.input_dir)    error "Provide --input_dir"
    if (!params.results_file) error "Provide --results_file"
    if (!params.kraken_db)    error "Provide --kraken_db"
    if (!params.resfinder_db) error "Provide --resfinder_db"

    def processed = get_processed_samples(params.results_file)
    log.info "Skipping ${processed.size()} already-processed samples."

    // Input discovery — Mode A: barcode subdirs, Mode B: flat fastq files
    def input_dir = file(params.input_dir)
    if (!input_dir.exists()) error "input_dir does not exist: ${params.input_dir}"

    def subdirs = input_dir.listFiles()?.findAll { d ->
        d.isDirectory() && d.listFiles()?.any { f ->
            f.name.endsWith('.fastq.gz') || f.name.endsWith('.fastq')
        }
    }?.sort() ?: []

    def flat_fqs = input_dir.listFiles()?.findAll { f ->
        f.name.endsWith('.fastq.gz') || f.name.endsWith('.fastq')
    }?.sort() ?: []

    def samples_map = [:]

    // Mode A: barcode subdirectories — each subdir name becomes a sample
    if (subdirs) {
        log.info "Mode A: found ${subdirs.size()} barcode subdirectories"
        subdirs.each { dir ->
            samples_map[dir.name] = dir.listFiles().findAll { f ->
                f.name.endsWith('.fastq.gz') || f.name.endsWith('.fastq')
            }.sort()
        }
    }

    // Mode B: flat fastq files — always runs; adds any sample not already found in Mode A
    if (flat_fqs) {
        def n = 0
        flat_fqs.each { f ->
            def prefix = f.name
                .replaceAll(/\.fastq\.gz$/, '').replaceAll(/\.fastq$/, '')
                .replaceAll(/_\d+$/, '')
            if (!samples_map.containsKey(prefix)) {
                samples_map[prefix] = [f]
                n += 1
            }
        }
        if (n > 0) log.info "Mode B: found ${n} flat fastq file(s)"
    }

    if (!samples_map) {
        error "No fastq(.gz) files found in ${input_dir}"
    }

    def new_samples = samples_map.findAll { name, fqs -> !processed.contains(name) }

    if (params.sample_name) {
        // Running as an SGE array task — restrict to one sample
        if (!samples_map.containsKey(params.sample_name)) {
            error "sample_name '${params.sample_name}' not found in input_dir: ${params.input_dir}"
        }
        new_samples = new_samples.findAll { name, fqs -> name == params.sample_name }
        if (!new_samples) {
            log.info "Sample '${params.sample_name}' already processed — skipping."
            return
        }
    } else if (params.sample_list) {
        // Explicit sample list ('#' and blanks ignored). Naming a sample overrides the
        // already-processed check.
        def list_file = file(params.sample_list)
        if (!list_file.exists()) error "sample_list not found: ${params.sample_list}"
        def wanted = list_file.readLines()
                              .collect { it.trim() }
                              .findAll { it && !it.startsWith('#') } as Set
        if (!wanted) error "sample_list is empty: ${params.sample_list}"

        def unknown = wanted.findAll { !samples_map.containsKey(it) }
        if (unknown) {
            log.warn "sample_list: ${unknown.size()} name(s) have no reads in " +
                     "${params.input_dir} and will be skipped: " +
                     "${unknown.sort().take(10).join(', ')}${unknown.size() > 10 ? ' ...' : ''}"
        }
        new_samples = samples_map.findAll { name, fqs -> wanted.contains(name) }
        if (!new_samples) {
            error "None of the ${wanted.size()} samples in ${params.sample_list} " +
                  "were found in ${params.input_dir}"
        }
        def redoing = new_samples.keySet().findAll { processed.contains(it) }
        log.info "sample_list: ${new_samples.size()} of ${wanted.size()} requested " +
                 "(${redoing.size()} already in the results file and will be reprocessed)."
    } else if (!new_samples) {
        error "All samples already in results file — nothing to do."
    }

    log.info "Samples to process: ${new_samples.size()} (${new_samples.keySet().join(', ')})"
    ch_raw = Channel.from(new_samples.collect { name, fqs -> [name, fqs] })

    // MALDI-TOF species map (optional — enables contamination detection)
    // Read the CSV in Groovy directly so Nextflow doesn't try to stage it.
    def ch_maldi_genus   = Channel.empty()
    def ch_maldi_species = Channel.empty()
    if (params.maldi_csv) {
        def maldi_file = file(params.maldi_csv)
        if (!maldi_file.exists()) error "MALDI CSV not found: ${params.maldi_csv}"
        def maldi_entries = []
        maldi_file.readLines().drop(1).each { line ->
            def cols = line.split(',', 2)
            if (cols.size() >= 2) {
                def sname = cols[0].trim()
                def sp    = cols[1].trim()
                def genus = sp ? sp.split()[0] : null
                if (sname && genus) maldi_entries << [sname, sp, genus]
            }
        }
        // Restrict to the samples this run is processing: otherwise the unmatched CSV
        // entries produce remainder tuples that crash the ch_genus map closure.
        def running = new_samples.keySet()
        maldi_entries = maldi_entries.findAll { e -> running.contains(e[0]) }
        if (maldi_entries) {
            ch_maldi_genus   = Channel.from(maldi_entries.collect { e -> tuple(e[0], e[2]) })
            ch_maldi_species = Channel.from(maldi_entries.collect { e -> tuple(e[0], e[1]) })
        }
    }

    // 1. Concatenate reads
    CONCAT_READS(ch_raw)

    // 2. NanoPlot read QC
    ch_needs_nanoplot = CONCAT_READS.out.raw_reads
        .filter { s, fq -> !file("${params.outdir}/${s}/reads/nanoplot/NanoStats.txt").exists() }
    ch_existing_nanoplot = CONCAT_READS.out.raw_reads
        .filter { s, fq -> file("${params.outdir}/${s}/reads/nanoplot/NanoStats.txt").exists() }
        .map    { s, fq -> tuple(s, file("${params.outdir}/${s}/reads/nanoplot/NanoStats.txt")) }

    NANOPLOT(ch_needs_nanoplot)
    ch_nanostats = NANOPLOT.out.stats.mix(ch_existing_nanoplot)

    ch_needs_read_stats    = ch_nanostats
        .filter { s, f -> !file("${params.outdir}/${s}/reads/${s}.read_stats.tsv").exists() }
    ch_existing_read_stats = ch_nanostats
        .filter { s, f -> file("${params.outdir}/${s}/reads/${s}.read_stats.tsv").exists() }
        .map    { s, f -> tuple(s, file("${params.outdir}/${s}/reads/${s}.read_stats.tsv")) }

    PARSE_READ_STATS(ch_needs_read_stats)
    ch_read_stats = PARSE_READ_STATS.out.mix(ch_existing_read_stats)

    // 3. Kraken2 on reads
    ch_needs_kr_reads = CONCAT_READS.out.reads
        .filter { s, fq -> !file("${params.outdir}/${s}/reads/kraken2/${s}_reads.kraken2.parsed").exists() }
    ch_existing_kr_reads = CONCAT_READS.out.reads
        .filter { s, fq -> file("${params.outdir}/${s}/reads/kraken2/${s}_reads.kraken2.parsed").exists() }
        .map    { s, fq -> tuple(s, file("${params.outdir}/${s}/reads/kraken2/${s}_reads.kraken2.parsed")) }

    KRAKEN2_READS(ch_needs_kr_reads, 'reads', 'reads', file(params.kraken_db))
    ch_kraken_reads = KRAKEN2_READS.out.parsed.mix(ch_existing_kr_reads)

    // 4. Filtlong read filtering + Unicycler assembly
    ch_needs_asm = CONCAT_READS.out.reads
        .filter { s, fq -> !file("${params.outdir}/${s}/assembly/flye/assembly.fasta").exists() }
    ch_existing_asm = CONCAT_READS.out.reads
        .filter { s, fq -> file("${params.outdir}/${s}/assembly/flye/assembly.fasta").exists() }

    FILTLONG(ch_needs_asm)
    FLYE(FILTLONG.out.reads)
    ch_assembly = FLYE.out.assembly
        .mix(ch_existing_asm.map { s, fq -> tuple(s, file("${params.outdir}/${s}/assembly/flye/assembly.fasta")) })
    ch_assembly_info = FLYE.out.info
        .mix(ch_existing_asm.map { s, fq -> tuple(s, file("${params.outdir}/${s}/assembly/flye/assembly_info.txt")) })

    // 4b. Contig-level genus filtering (only when MALDI CSV supplied)
    // Genus filtering is applied to the chromosome only, after contig separation.
    // ch_maldi_flag_val / ch_removed_len_val are populated in section 6a below.
    def ch_maldi_flag_val
    def ch_removed_len_val

    // 5. mob_recon (full assembly; genus filtering is chromosome-only, see 6a)
    ch_needs_mob = ch_assembly
        .filter { s, asm -> !file("${params.outdir}/${s}/assembly/mob_recon/mob_recon_out/contig_report.txt").exists() }
    ch_existing_mob = ch_assembly
        .filter { s, asm -> file("${params.outdir}/${s}/assembly/mob_recon/mob_recon_out/contig_report.txt").exists() }
        .map    { s, asm -> tuple(s, file("${params.outdir}/${s}/assembly/mob_recon/mob_recon_out")) }

    MOB_RECON(ch_needs_mob)
    ch_mob_dir = MOB_RECON.out.mob_dir.mix(ch_existing_mob)

    // 6. Separate contigs (chromosome / plasmids)
    // rescued_mge_contigs.tsv is required too: samples separated before the MGE rescue
    // existed lack it and must be re-separated.
    ch_needs_sep = ch_assembly
        .filter { s, asm ->
            !file("${params.outdir}/${s}/assembly/contigs/chromosome.fasta").exists() ||
            !file("${params.outdir}/${s}/assembly/contigs/contig_to_mob_annotation.tsv").exists() ||
            !file("${params.outdir}/${s}/assembly/contigs/rescued_mge_contigs.tsv").exists()
        }
    ch_existing_sep = ch_assembly
        .filter { s, asm ->
            file("${params.outdir}/${s}/assembly/contigs/chromosome.fasta").exists() &&
            file("${params.outdir}/${s}/assembly/contigs/contig_to_mob_annotation.tsv").exists() &&
            file("${params.outdir}/${s}/assembly/contigs/rescued_mge_contigs.tsv").exists()
        }

    // assembly_info carries Flye's per-contig coverage.
    SEPARATE_CONTIGS(ch_mob_dir.join(ch_needs_sep).join(ch_assembly_info))

    ch_chromosome = SEPARATE_CONTIGS.out.chromosome
        .mix(ch_existing_sep.map { s, asm -> tuple(s, file("${params.outdir}/${s}/assembly/contigs/chromosome.fasta")) })

    // Pre-filter chromosome, kept so chromosomal contamination stays visible in the
    // report; ch_chromosome is reassigned to the filtered version in 6a.
    ch_chromosome_raw = ch_chromosome
    ch_pla_merged = SEPARATE_CONTIGS.out.plasmids_combined
        .mix(ch_existing_sep.map { s, asm -> tuple(s, file("${params.outdir}/${s}/assembly/contigs/plasmids_all.fasta")) })
    ch_contig_map = SEPARATE_CONTIGS.out.contig_plasmid_map
        .mix(ch_existing_sep.map { s, asm -> tuple(s, file("${params.outdir}/${s}/assembly/contigs/contig_to_plasmid.tsv")) })
    ch_pla_summary = SEPARATE_CONTIGS.out.pla_summary
        .mix(ch_existing_sep.map { s, asm -> tuple(s, file("${params.outdir}/${s}/assembly/plasmid_summary.txt")) })

    // 6a. Genus contig filtering on chromosome only (when MALDI CSV supplied)
    if (params.maldi_csv) {
        def ch_needs_filter = ch_chromosome
            .join(ch_maldi_genus)
            .filter { s, chr, genus -> !filter_done(params.outdir, s) }

        def ch_existing_filter_chr = ch_chromosome
            .join(ch_maldi_genus)
            .filter { s, chr, genus -> filter_done(params.outdir, s) }
            .map    { s, chr, genus -> tuple(s, file("${params.outdir}/${s}/chromosome/contigs_filtered.fasta")) }

        def ch_existing_flag_val = ch_chromosome
            .join(ch_maldi_genus)
            .filter { s, chr, genus -> filter_done(params.outdir, s) }
            .map    { s, chr, genus ->
                def f = file("${params.outdir}/${s}/chromosome/maldi_contigs_used.txt")
                tuple(s, f.text.trim())
            }

        // Removed-contig lengths from a previous run; older runs predate this file.
        def ch_existing_removed_val = ch_chromosome
            .join(ch_maldi_genus)
            .filter { s, chr, genus -> filter_done(params.outdir, s) }
            .map    { s, chr, genus ->
                def f = file("${params.outdir}/${s}/chromosome/removed_contig_lengths.txt")
                tuple(s, f.exists() ? f.text.trim() : 'NA')
            }

        // Query contigs from a previous run; older runs predate these files.
        def ch_existing_query = ch_chromosome
            .join(ch_maldi_genus)
            .filter { s, chr, genus -> file("${params.outdir}/${s}/chromosome/query_contigs.fasta").exists() }
            .map    { s, chr, genus -> tuple(s, file("${params.outdir}/${s}/chromosome/query_contigs.fasta")) }

        def ch_existing_chr_species_val = ch_chromosome
            .join(ch_maldi_genus)
            .filter { s, chr, genus -> filter_done(params.outdir, s) }
            .map    { s, chr, genus ->
                def f = file("${params.outdir}/${s}/chromosome/chromosome_contig_species.txt")
                tuple(s, f.exists() ? f.text.trim() : 'NA')
            }

        def ch_existing_query_val = ch_chromosome
            .join(ch_maldi_genus)
            .filter { s, chr, genus -> filter_done(params.outdir, s) }
            .map    { s, chr, genus ->
                def f = file("${params.outdir}/${s}/chromosome/query_contig_info.txt")
                tuple(s, f.exists() ? (f.text.trim().split('\t') as List)[0] : '0')
            }

        FILTER_GENUS_CONTIGS(ch_needs_filter, file(params.kraken_db))

        def ch_filtered_chr = FILTER_GENUS_CONTIGS.out.assembly.mix(ch_existing_filter_chr)
        def ch_new_flag_val = FILTER_GENUS_CONTIGS.out.flag.map { s, f -> tuple(s, f.text.trim()) }
        def ch_new_removed_val = FILTER_GENUS_CONTIGS.out.removed_lengths.map { s, f -> tuple(s, f.text.trim()) }
        def ch_new_query_val = FILTER_GENUS_CONTIGS.out.query_info
            .map { s, f -> tuple(s, (f.text.trim().split('\t') as List)[0] ?: '0') }
        def ch_new_chr_species_val = FILTER_GENUS_CONTIGS.out.chr_species
            .map { s, f -> tuple(s, f.text.trim() ?: 'NA') }

        // Samples with no MALDI genus entry keep their chromosome unchanged
        def ch_no_maldi_chr = ch_chromosome
            .join(ch_maldi_genus, remainder: true)
            .filter { s, chr, genus -> genus == null }
            .map    { s, chr, genus -> tuple(s, chr) }

        ch_chromosome     = ch_filtered_chr.mix(ch_no_maldi_chr)
        ch_maldi_flag_val = ch_new_flag_val
            .mix(ch_existing_flag_val)
            .mix(ch_no_maldi_chr.map { s, chr -> tuple(s, 'NA') })
        ch_removed_len_val = ch_new_removed_val
            .mix(ch_existing_removed_val)
            .mix(ch_no_maldi_chr.map { s, chr -> tuple(s, 'NA') })

        // Off-genus contigs below the element-size limit: screened for AMR, excluded from
        // chromosome and plasmid metrics. Empty FASTAs dropped.
        ch_query = FILTER_GENUS_CONTIGS.out.query
            .mix(ch_existing_query)
            .filter { s, fa -> fa.size() > 0 }
        ch_query_count_val = ch_new_query_val
            .mix(ch_existing_query_val)
            .mix(ch_no_maldi_chr.map { s, chr -> tuple(s, '0') })
        ch_chr_species_val = ch_new_chr_species_val
            .mix(ch_existing_chr_species_val)
            .mix(ch_no_maldi_chr.map { s, chr -> tuple(s, 'NA') })
    } else {
        ch_maldi_flag_val  = ch_chromosome.map { s, chr -> tuple(s, 'NA') }
        ch_removed_len_val = ch_chromosome.map { s, chr -> tuple(s, 'NA') }
        // No MALDI genus to filter against, so no query contigs arise.
        ch_query           = Channel.empty()
        ch_query_count_val = ch_chromosome.map { s, chr -> tuple(s, '0') }
        // No genus filtering ran, so no per-contig species summary exists
        ch_chr_species_val = ch_chromosome.map { s, chr -> tuple(s, 'NA') }
    }

    // 7. Depth of coverage (reads mapped to MALDI-filtered chromosome + plasmids)
    ch_needs_depth = CONCAT_READS.out.reads
        .join(ch_chromosome)
        .join(ch_pla_merged)
        .join(ch_contig_map)
        .filter { s, reads, chr, pla, cmap -> !file("${params.outdir}/${s}/assembly/${s}.depth.tsv").exists() }
    ch_existing_depth = CONCAT_READS.out.reads
        .join(ch_chromosome)
        .filter { s, reads, chr -> file("${params.outdir}/${s}/assembly/${s}.depth.tsv").exists() }
        .map    { s, reads, chr -> tuple(s, file("${params.outdir}/${s}/assembly/${s}.depth.tsv")) }

    CALC_DEPTH(ch_needs_depth)
    ch_depth = CALC_DEPTH.out.depth.mix(ch_existing_depth)

    // 8. Kraken2 on the unfiltered chromosome, so chromosomal contamination stays
    //    visible. Plasmids are excluded: they legitimately come from other species.
    ch_needs_kr_asm = ch_chromosome_raw
        .filter { s, asm -> !file("${params.outdir}/${s}/chromosome/kraken2/${s}_assembly.kraken2.parsed").exists() }
    ch_existing_kr_asm = ch_chromosome_raw
        .filter { s, asm -> file("${params.outdir}/${s}/chromosome/kraken2/${s}_assembly.kraken2.parsed").exists() }
        .map    { s, asm -> tuple(s, file("${params.outdir}/${s}/chromosome/kraken2/${s}_assembly.kraken2.parsed")) }

    KRAKEN2_ASSEMBLY(ch_needs_kr_asm, 'assembly', 'chromosome', file(params.kraken_db))
    ch_kraken_asm = KRAKEN2_ASSEMBLY.out.parsed.mix(ch_existing_kr_asm)

    // Parse genus: prefer assembly Kraken result; fall back to reads Kraken, then MALDI genus.
    // The MALDI fallback covers assemblies too fragmented for Kraken2 to classify.
    ch_genus = ch_kraken_asm
        .join(ch_kraken_reads, remainder: true)
        .join(ch_maldi_genus,  remainder: true)
        .map { row ->
            // row is [s, f_asm, f_reads, maldi_g], but collapses to 3 elements when a sample
            // exists only in ch_maldi_genus. Handle both lengths.
            def s       = row[0]
            def f_asm   = row.size() >= 4 ? row[1] : null
            def f_reads = row.size() >= 4 ? row[2] : null
            def maldi_g = row.size() >= 4 ? row[3] : (row.size() == 3 ? row[2] : null)
            def g = f_asm   ? f_asm.text.split('\t')[0].trim()   : 'NA'
            if (g == 'NA' && f_reads) g = f_reads.text.split('\t')[0].trim()
            if (g == 'NA' && maldi_g && maldi_g != 'NA') g = maldi_g
            tuple(s, g ?: 'NA')
        }

    // 9. CheckM2
    def ch_needs_chkm2_all = ch_assembly
        .filter { s, asm -> !file("${params.outdir}/${s}/assembly/checkm2/checkm2_output/quality_report.tsv").exists() }
        .map    { s, asm -> tuple(s, asm, 'assembly') }
    def ch_existing_chkm2_all = ch_assembly
        .filter { s, asm -> file("${params.outdir}/${s}/assembly/checkm2/checkm2_output/quality_report.tsv").exists() }
        .map    { s, asm -> tuple(s, file("${params.outdir}/${s}/assembly/checkm2/checkm2_output/quality_report.tsv")) }

    CHECKM2_ALL(ch_needs_chkm2_all)
    def ch_checkm2_all = CHECKM2_ALL.out.report.mix(ch_existing_chkm2_all)

    def ch_needs_chkm2 = ch_chromosome
        .filter { s, asm -> !file("${params.outdir}/${s}/chromosome/checkm2/checkm2_output/quality_report.tsv").exists() }
        .map    { s, asm -> tuple(s, asm, 'chromosome') }
    def ch_existing_chkm2 = ch_chromosome
        .filter { s, asm -> file("${params.outdir}/${s}/chromosome/checkm2/checkm2_output/quality_report.tsv").exists() }
        .map    { s, asm -> tuple(s, file("${params.outdir}/${s}/chromosome/checkm2/checkm2_output/quality_report.tsv")) }

    CHECKM2_CHR(ch_needs_chkm2)
    def ch_checkm2 = CHECKM2_CHR.out.report.mix(ch_existing_chkm2)

    // 9a. BUSCO (output is a directory)
    ch_needs_busco = ch_chromosome
        .filter { s, asm -> !new File("${params.outdir}/${s}/chromosome/busco/busco_out").isDirectory() }
    ch_existing_busco = ch_chromosome
        .filter { s, asm -> new File("${params.outdir}/${s}/chromosome/busco/busco_out").isDirectory() }
        .map    { s, asm ->
            def d = new File("${params.outdir}/${s}/chromosome/busco/busco_out")
            def f = d.listFiles()?.find { it.name.startsWith('short_summary') && it.name.endsWith('.txt') }
            tuple(s, file(f ? f.path : d.path))
        }

    BUSCO(ch_needs_busco)
    ch_busco = BUSCO.out.summary.mix(ch_existing_busco)

    // 10. MLST
    ch_needs_mlst = ch_chromosome
        .filter { s, asm -> !file("${params.outdir}/${s}/chromosome/mlst/${s}.mlst.tsv").exists() }
    ch_existing_mlst = ch_chromosome
        .filter { s, asm -> file("${params.outdir}/${s}/chromosome/mlst/${s}.mlst.tsv").exists() }
        .map    { s, asm -> tuple(s, file("${params.outdir}/${s}/chromosome/mlst/${s}.mlst.tsv")) }

    MLST(ch_needs_mlst)
    ch_mlst = MLST.out.results.mix(ch_existing_mlst)

    // 11. Chromosome AMR tools
    ch_needs_amrf_chr = ch_chromosome
        .filter { s, asm -> !file("${params.outdir}/${s}/chromosome/amrfinder/${s}_chromosome.amrfinder.tsv").exists() }
    ch_existing_amrf_chr = ch_chromosome
        .filter { s, asm -> file("${params.outdir}/${s}/chromosome/amrfinder/${s}_chromosome.amrfinder.tsv").exists() }
        .map    { s, asm -> tuple(s, file("${params.outdir}/${s}/chromosome/amrfinder/${s}_chromosome.amrfinder.tsv")) }
    AMRFINDER_CHR(ch_needs_amrf_chr, 'chromosome', 'chromosome')
    ch_amrf_chr = AMRFINDER_CHR.out.results.mix(ch_existing_amrf_chr)

    ch_needs_resf_chr = ch_chromosome
        .filter { s, asm -> !new File("${params.outdir}/${s}/chromosome/resfinder/resfinder_chromosome").isDirectory() }
    ch_existing_resf_chr = ch_chromosome
        .filter { s, asm -> new File("${params.outdir}/${s}/chromosome/resfinder/resfinder_chromosome").isDirectory() }
        .map    { s, asm -> tuple(s, file("${params.outdir}/${s}/chromosome/resfinder/resfinder_chromosome")) }
    RESFINDER_CHR(ch_needs_resf_chr, 'chromosome', 'chromosome', file(params.resfinder_db))
    ch_resf_chr = RESFINDER_CHR.out.results.mix(ch_existing_resf_chr)

    ch_needs_rgi_chr = ch_chromosome
        .filter { s, asm -> !file("${params.outdir}/${s}/chromosome/rgi/${s}_chromosome.rgi.txt").exists() }
    ch_existing_rgi_chr = ch_chromosome
        .filter { s, asm -> file("${params.outdir}/${s}/chromosome/rgi/${s}_chromosome.rgi.txt").exists() }
        .map    { s, asm -> tuple(s, file("${params.outdir}/${s}/chromosome/rgi/${s}_chromosome.rgi.txt")) }
    RGI_CHR(ch_needs_rgi_chr, 'chromosome', 'chromosome')
    ch_rgi_chr = RGI_CHR.out.results.mix(ch_existing_rgi_chr)

    ch_needs_cefi_chr = ch_chromosome
        .filter { s, asm -> !file("${params.outdir}/${s}/chromosome/cefiderocolfinder/${s}_chromosome.cefiderocol.tsv").exists() }
    ch_existing_cefi_chr = ch_chromosome
        .filter { s, asm -> file("${params.outdir}/${s}/chromosome/cefiderocolfinder/${s}_chromosome.cefiderocol.tsv").exists() }
        .map    { s, asm -> tuple(s, file("${params.outdir}/${s}/chromosome/cefiderocolfinder/${s}_chromosome.cefiderocol.tsv")) }
    CEFI_CHR(ch_needs_cefi_chr, 'chromosome', 'chromosome')
    ch_cefi_chr = CEFI_CHR.out.results.mix(ch_existing_cefi_chr)

    PARSE_AMR_CHR(ch_amrf_chr.join(ch_resf_chr).join(ch_rgi_chr).join(ch_cefi_chr), 'chr')

    // 12. Species-specific typing (Kleborate / ParEx)
    // Kleborate v3 has one module set per species complex; the KpSC preset on a KoSC
    // genome types it against the wrong scheme and writes a different filename.
    ch_klebs_input = ch_chromosome.join(ch_kraken_asm).join(ch_genus)
        .map { s, asm, kfile, genus ->
            def sp = kfile ? (kfile.text.split('\t').size() > 3
                              ? kfile.text.split('\t')[3].trim().toLowerCase() : '') : ''
            def preset = 'none'
            if (sp =~ /oxytoca|michiganensis|grimontii|pasteurii|huaxiensis|spallanzanii/) {
                preset = 'kosc'
            }
            else if (sp =~ /klebsiella\s+(pneumoniae|variicola|africana|quasi)/) {
                preset = 'kpsc'
            }
            else if (sp =~ /escherichia|shigella/) {
                preset = 'escherichia'
            }
            // Anything else (e.g. K. aerogenes) has no module and stays 'none'.
            tuple(s, asm, preset, genus)
        }
        .filter { s, asm, preset, genus -> preset != 'none' }
        .map { s, asm, preset, genus -> tuple(s, asm, preset) }

    ch_needs_klebs = ch_klebs_input
        .filter { s, asm, preset -> !file("${params.outdir}/${s}/chromosome/kleborate/${s}.kleborate.txt").exists() }
    ch_existing_klebs = ch_klebs_input
        .filter { s, asm, preset -> file("${params.outdir}/${s}/chromosome/kleborate/${s}.kleborate.txt").exists() }
        .map    { s, asm, preset -> tuple(s, file("${params.outdir}/${s}/chromosome/kleborate/${s}.kleborate.txt")) }
    KLEBORATE(ch_needs_klebs)
    ch_kleborate = KLEBORATE.out.results.mix(ch_existing_klebs)

    ch_pseudo_input = ch_chromosome.join(ch_genus)
        .filter { s, asm, genus -> genus.toLowerCase().contains('pseudomonas') }
        .map    { s, asm, genus -> tuple(s, asm) }

    ch_needs_parex = ch_pseudo_input
        .filter { s, asm -> !file("${params.outdir}/${s}/chromosome/parex/${s}.parex.txt").exists() }
    ch_existing_parex = ch_pseudo_input
        .filter { s, asm -> file("${params.outdir}/${s}/chromosome/parex/${s}.parex.txt").exists() }
        .map    { s, asm -> tuple(s, file("${params.outdir}/${s}/chromosome/parex/${s}.parex.txt")) }
    PAREX(ch_needs_parex)
    ch_parex = PAREX.out.results.mix(ch_existing_parex)

    // 13. Plasmid AMR tools
    ch_needs_amrf_pla = ch_pla_merged
        .filter { s, fa -> !file("${params.outdir}/${s}/plasmids/amrfinder/${s}_plasmid.amrfinder.tsv").exists() }
    ch_existing_amrf_pla = ch_pla_merged
        .filter { s, fa -> file("${params.outdir}/${s}/plasmids/amrfinder/${s}_plasmid.amrfinder.tsv").exists() }
        .map    { s, fa -> tuple(s, file("${params.outdir}/${s}/plasmids/amrfinder/${s}_plasmid.amrfinder.tsv")) }
    AMRFINDER_PLA(ch_needs_amrf_pla, 'plasmid', 'plasmids')
    ch_amrf_pla = AMRFINDER_PLA.out.results.mix(ch_existing_amrf_pla)

    ch_needs_resf_pla = ch_pla_merged
        .filter { s, fa -> !new File("${params.outdir}/${s}/plasmids/resfinder/resfinder_plasmid").isDirectory() }
    ch_existing_resf_pla = ch_pla_merged
        .filter { s, fa -> new File("${params.outdir}/${s}/plasmids/resfinder/resfinder_plasmid").isDirectory() }
        .map    { s, fa -> tuple(s, file("${params.outdir}/${s}/plasmids/resfinder/resfinder_plasmid")) }
    RESFINDER_PLA(ch_needs_resf_pla, 'plasmid', 'plasmids', file(params.resfinder_db))
    ch_resf_pla = RESFINDER_PLA.out.results.mix(ch_existing_resf_pla)

    ch_needs_rgi_pla = ch_pla_merged
        .filter { s, fa -> !file("${params.outdir}/${s}/plasmids/rgi/${s}_plasmid.rgi.txt").exists() }
    ch_existing_rgi_pla = ch_pla_merged
        .filter { s, fa -> file("${params.outdir}/${s}/plasmids/rgi/${s}_plasmid.rgi.txt").exists() }
        .map    { s, fa -> tuple(s, file("${params.outdir}/${s}/plasmids/rgi/${s}_plasmid.rgi.txt")) }
    RGI_PLA(ch_needs_rgi_pla, 'plasmid', 'plasmids')
    ch_rgi_pla = RGI_PLA.out.results.mix(ch_existing_rgi_pla)

    ch_needs_cefi_pla = ch_pla_merged
        .filter { s, fa -> !file("${params.outdir}/${s}/plasmids/cefiderocolfinder/${s}_plasmid.cefiderocol.tsv").exists() }
    ch_existing_cefi_pla = ch_pla_merged
        .filter { s, fa -> file("${params.outdir}/${s}/plasmids/cefiderocolfinder/${s}_plasmid.cefiderocol.tsv").exists() }
        .map    { s, fa -> tuple(s, file("${params.outdir}/${s}/plasmids/cefiderocolfinder/${s}_plasmid.cefiderocol.tsv")) }
    CEFI_PLA(ch_needs_cefi_pla, 'plasmid', 'plasmids')
    ch_cefi_pla = CEFI_PLA.out.results.mix(ch_existing_cefi_pla)

    PARSE_AMR_PLA(ch_amrf_pla.join(ch_resf_pla).join(ch_rgi_pla).join(ch_cefi_pla), 'pla')

    // 13a. Query contigs: off-genus, below the element-size limit.  Screened
    // for AMR, but absent from chromosome and plasmid metrics. Not gated on published
    // outputs: query contigs are rare, so Nextflow's own cache suffices for resume.
    AMRFINDER_QRY(ch_query, 'query', 'query_contigs')
    RESFINDER_QRY(ch_query, 'query', 'query_contigs', file(params.resfinder_db))
    RGI_QRY(ch_query,       'query', 'query_contigs')
    CEFI_QRY(ch_query,      'query', 'query_contigs')

    PARSE_AMR_QRY(
        AMRFINDER_QRY.out.results
            .join(RESFINDER_QRY.out.results)
            .join(RGI_QRY.out.results)
            .join(CEFI_QRY.out.results),
        'query'
    )

    // 14. Combine chr + pla (+ query) AMR
    ch_combine_amr_in = PARSE_AMR_CHR.out.parsed
        .join(PARSE_AMR_PLA.out.parsed)
        .join(ch_contig_map)
        .join(PARSE_AMR_QRY.out.parsed, remainder: true)
        .filter { row -> row instanceof List && row.size() == 5 && row[1] != null }
        .map { s, chr_amr, pla_amr, cmap, qry_amr ->
            tuple(s, chr_amr, pla_amr, cmap, qry_amr ?: file('NO_FILE_QRYAMR'))
        }
    COMBINE_AMR(ch_combine_amr_in)

    // 15. Genome evaluation summary
    PARSE_GENOME_EVAL(
        ch_checkm2                // filtered chromosome: completeness etc.
            .join(ch_checkm2_all) // full assembly: contamination
            .join(ch_busco)
            .join(ch_mlst)
            .join(ch_chromosome)
            .join(ch_assembly_info)
    )

    // 16. Per-sample metrics JSON
    // Resolve MALDI species per sample (NA for samples absent from the MALDI CSV)
    def ch_maldi_per_sample = ch_read_stats
        .map   { s, rs -> tuple(s, 'NA') }
        .join  (ch_maldi_species, remainder: true)
        .map   { s, na_val, sp -> tuple(s, sp ?: 'NA') }

    ch_compile = ch_read_stats
        .join(ch_kraken_reads)
        .join(ch_kraken_asm,          remainder: true)
        .join(PARSE_GENOME_EVAL.out,  remainder: true)
        .join(ch_pla_summary)
        .join(ch_kleborate,           remainder: true)
        .join(ch_parex,               remainder: true)
        .join(ch_depth,               remainder: true)
        .join(ch_maldi_flag_val)
        .join(ch_removed_len_val)
        .join(ch_query_count_val)
        .join(ch_chr_species_val)
        .join(ch_maldi_per_sample)
        // A failed PARSE_READ_STATS lets `remainder: true` re-inject a collapsed tuple,
        // which would crash the map closure and abort the whole run. Drop such rows.
        .filter { row -> row instanceof List && row.size() == 14 && row[1] != null }
        .map { s, rs, kr_r, kr_a, geval, pla_sum, klebs, parex, depth, maldi_flag_v, removed_len_v, query_n_v, chr_sp_v, maldi_sp ->
            tuple(s, rs, kr_r,
                  kr_a    ?: file('NO_FILE_KRASM'),
                  geval   ?: file('NO_FILE_GEVAL'),
                  pla_sum,
                  klebs   ?: file('NO_FILE_KLEBS'),
                  parex   ?: file('NO_FILE_PAREX'),
                  depth   ?: file('NO_FILE_DEPTH'),
                  maldi_flag_v  ?: 'NA',
                  removed_len_v ?: 'NA',
                  query_n_v     ?: '0',
                  chr_sp_v      ?: 'NA',
                  maldi_sp)
        }

    COMPILE_SAMPLE(ch_compile)

    // 17. Write wide results TSV (runs once after all samples finish)
    if (!params.skip_write_results) {
        WRITE_RESULTS(
            COMPILE_SAMPLE.out.metrics.collect(),
            COMBINE_AMR.out.combined.map { s, f -> f }.collect(),
            file(params.results_file).toAbsolutePath().toString()
        )
    }
}
