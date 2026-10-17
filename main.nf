nextflow.enable.dsl=2

process PREPARE_GENE_VCFS {
    tag 'chr20 gene subsets'
    publishDir "${params.outdir}/prepared", mode: 'copy', overwrite: false

    input:
    path prepare_script
    path regions
    path truth
    path deepvariant
    path gatk
    val padding

    output:
    path 'subset', emit: subset

    script:
    """
    python3 "${prepare_script}" \
        --regions "${regions}" \
        --giab "${truth}" \
        --deepvariant "${deepvariant}" \
        --gatk "${gatk}" \
        --padding ${padding} \
        --output-dir subset
    """
}

process CHECK_REF_ALLELES {
    tag 'validate REF alleles'
    publishDir "${params.outdir}/reference_validation", mode: 'copy', overwrite: false

    input:
    path check_script
    path subset
    path reference
    path reference_index

    output:
    path 'reference_check.json', emit: report

    script:
    """
    python3 "${check_script}" \
        --reference "${reference}" \
        --input-dir "${subset}" \
        --report reference_check.json
    """
}

process PREPARE_COMPARISON_INPUTS {
    tag 'prepare AD-free copies'
    publishDir "${params.outdir}/comparison_inputs", mode: 'copy', overwrite: false

    input:
    path prepare_script
    path subset
    path reference_report

    output:
    path 'noAD', emit: ready

    script:
    """
    mkdir checked
    cp "${subset}"/*.padded.vcf checked/
    cp "${reference_report}" checked/reference_check.json

    python3 "${prepare_script}" \
        --input-dir checked \
        --output-dir noAD
    """
}


process BENCHMARK_CALLSET {
    tag "${caller}"
    publishDir "${params.outdir}/benchmarks", mode: 'copy', overwrite: false

    // Podman'a bağlanan görev klasöründe gerçek dosyalar bulunsun.
    stageInMode 'copy'

    input:
    tuple val(caller), path(comparison_inputs)
    path reference
    path reference_index
    path confidence_regions
    path image_record

    output:
    tuple val(caller), path("${caller}"), emit: results

    script:
    """
    image_id=\$(python3 -c \
        'import json,sys; print(json.load(open(sys.argv[1]))[0]["Id"])' \
        "${image_record}")

    podman --storage-opt overlay.ignore_chown_errors=true \
        image exists "\$image_id"

    mkdir "${caller}"

    podman --storage-opt overlay.ignore_chown_errors=true \
        run --rm --cpus 1 --memory 1g \
        -e RTG_MEM=512m \
        -v "\$PWD":/work:z -w /work \
        "\$image_id" \
        /opt/hap.py/bin/hap.py \
        "${comparison_inputs}/giab.vcf" \
        "${comparison_inputs}/${caller}.vcf" \
        -r "${reference}" \
        -f "${confidence_regions}" \
        --engine vcfeval \
        --engine-vcfeval-path /opt/hap.py/libexec/rtg-tools-install/rtg \
        --gender male --threads 1 --pass-only \
        -o "${caller}/${caller}" \
        2>&1 | tee "${caller}/${caller}.log"

    test -s "${caller}/${caller}.summary.csv"
    test -s "${caller}/${caller}.vcf.gz"
    """
}


process QUANTIFY_PER_GENE {
    tag "${caller}"
    publishDir "${params.outdir}/per_gene", mode: 'copy', overwrite: false
    stageInMode 'copy'

    input:
    tuple val(caller), path(benchmark_dir), path(comparison_inputs)
    path reference
    path reference_index
    path confidence_regions
    path gene_regions
    path image_record

    output:
    tuple val(caller), path("${caller}_per_gene"), emit: results

    script:
    """
    image_id=\$(python3 -c \
        'import json,sys; print(json.load(open(sys.argv[1]))[0]["Id"])' \
        "${image_record}")

    mkdir "${caller}_per_gene"

    podman --storage-opt overlay.ignore_chown_errors=true \
        run --rm --cpus 1 --memory 1g \
        -v "\$PWD":/work:z -w /work \
        "\$image_id" \
        /opt/hap.py/bin/qfy.py \
        "${benchmark_dir}/${caller}.vcf.gz" \
        --type ga4gh \
        -r "${reference}" \
        -f "${confidence_regions}" \
        --adjust-conf-regions "${comparison_inputs}/giab.vcf" \
        --stratification "${gene_regions}/genes.tsv" \
        --threads 1 --no-roc \
        -o "${caller}_per_gene/${caller}" \
        2>&1 | tee "${caller}_per_gene/${caller}.log"

    test -s "${caller}_per_gene/${caller}.summary.csv"
    test -s "${caller}_per_gene/${caller}.extended.csv"
    """
}


process ANALYZE_GENE_CONTEXT {
    tag "${gene}"
    publishDir "${params.outdir}/context_analysis", mode: 'copy', overwrite: false

    input:
    path analysis_script
    path benchmark_dirs
    path gene_bed
    path repeat_bed
    val gene

    output:
    path "${gene}", emit: context

    script:
    """
    python3 "${analysis_script}" \
        --gene "${gene}" \
        --deepvariant deepvariant/deepvariant.vcf.gz \
        --gatk gatk/gatk.vcf.gz \
        --gene-bed "${gene_bed}" \
        --repeat-bed "${repeat_bed}" \
        --output-dir "${gene}"
    """
}

process PLOT_GENE_CONTEXT {
    tag "${gene}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: false

    input:
    path plot_script
    path context_dir
    val gene

    output:
    path 'figures', emit: figures

    script:
    """
    python3 "${plot_script}" \
        --gene "${gene}" \
        --input-dir "${context_dir}" \
        --output-dir figures
    """
}


process PREPARE_GENE_REGIONS {
    tag 'chr20 CMRG regions'
    publishDir "${params.outdir}", mode: 'copy', overwrite: false

    input:
    path region_script
    path coordinates
    path benchmark_bed
    path reference_index

    output:
    path 'regions', emit: regions

    script:
    """
    python3 "${region_script}" \
        --coordinates "${coordinates}" \
        --benchmark-bed "${benchmark_bed}" \
        --reference-index "${reference_index}" \
        --chromosome chr20 \
        --output-dir regions
    """
}

workflow {
    def manifestFile = file(params.input_manifest, checkIfExists: true)
    def manifest = new groovy.json.JsonSlurper().parse(manifestFile.toFile())

    def inputs = manifest.source_inputs.collectEntries { item ->
        [(item.id): file("${projectDir}/${item.local_path}", checkIfExists: true)]
    }

    def required = [
        'giab_vcf', 'deepvariant', 'gatk', 'reference_chr20',
        'giab_bed', 'cmrg_coordinates', 'tandem_repeats'
    ]
    required.each { id ->
        if (!inputs.containsKey(id)) {
            error "Missing input manifest entry: ${id}"
        }
    }

    PREPARE_GENE_REGIONS(
        file("${projectDir}/scripts/21_prepare_gene_regions.py", checkIfExists: true),
        inputs.cmrg_coordinates,
        inputs.giab_bed,
        file("${inputs.reference_chr20}.fai", checkIfExists: true)
    )

    // Single generated region set is shared by all downstream tasks.
    def generated_regions = PREPARE_GENE_REGIONS.out.regions
    def selected_gene_regions = generated_regions.map { directory ->
        directory.resolve('chr20_cmrg_genes.bed')
    }
    def scoring_regions = generated_regions.map { directory ->
        directory.resolve('chr20_cmrg_evaluable.bed')
    }
    def per_gene_regions = generated_regions.map { directory ->
        directory.resolve('per_gene')
    }
    def context_gene_bed = generated_regions.map { directory ->
        directory.resolve("per_gene/${params.context_gene}.bed")
    }

    PREPARE_GENE_VCFS(
        file("${projectDir}/scripts/04_prepare_gene_vcfs.py", checkIfExists: true),
        selected_gene_regions,
        inputs.giab_vcf,
        inputs.deepvariant,
        inputs.gatk,
        params.padding
    )

    CHECK_REF_ALLELES(
        file("${projectDir}/scripts/05_check_gene_reference.py", checkIfExists: true),
        PREPARE_GENE_VCFS.out.subset,
        inputs.reference_chr20,
        file("${inputs.reference_chr20}.fai", checkIfExists: true)
    )

    PREPARE_COMPARISON_INPUTS(
        file("${projectDir}/scripts/06_prepare_gene_benchmark_inputs.py", checkIfExists: true),
        PREPARE_GENE_VCFS.out.subset,
        CHECK_REF_ALLELES.out.report
    )

    def benchmark_jobs = PREPARE_COMPARISON_INPUTS.out.ready.flatMap { directory ->
        ['deepvariant', 'gatk'].collect { caller ->
            tuple(caller, directory)
        }
    }

    BENCHMARK_CALLSET(
        benchmark_jobs,
        inputs.reference_chr20,
        file("${inputs.reference_chr20}.fai", checkIfExists: true),
        scoring_regions,
        file("${projectDir}/docs/software/happy_image_inspect.json", checkIfExists: true)
    )

    def quantification_jobs = BENCHMARK_CALLSET.out.results.combine(
        PREPARE_COMPARISON_INPUTS.out.ready
    )

    QUANTIFY_PER_GENE(
        quantification_jobs,
        inputs.reference_chr20,
        file("${inputs.reference_chr20}.fai", checkIfExists: true),
        scoring_regions,
        per_gene_regions,
        file("${projectDir}/docs/software/happy_image_inspect.json", checkIfExists: true)
    )

    def context_benchmarks = BENCHMARK_CALLSET.out.results
        .map { caller, directory -> directory }
        .collect()

    ANALYZE_GENE_CONTEXT(
        file("${projectDir}/scripts/09_analyze_gene_context.py", checkIfExists: true),
        context_benchmarks,
        context_gene_bed,
        inputs.tandem_repeats,
        params.context_gene
    )

    PLOT_GENE_CONTEXT(
        file("${projectDir}/scripts/10_plot_gene_context.py", checkIfExists: true),
        ANALYZE_GENE_CONTEXT.out.context,
        params.context_gene
    )
}
