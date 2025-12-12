rule extract_prophage_region:
    input: 
        seq=SEQUENCE,
        bed_file=os.path.join(RESULTS_DIR, "results", "final_coordinates_0-based.bed")
    output:
        extracted_seq=os.path.join(RESULTS_DIR, "results", "prophage_region_sequences.fa")
    threads:
        config["resources"]["small_cpu"]
    resources:  
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "bedtools.yml")
    shell:
        """
        if grep -q "NO PROPHAGE REGIONS FOUND" {input.bed_file}; then
            touch {output.extracted_seq}
        else
            bedtools getfasta -fi {input.seq} \
                -bed {input.bed_file} \
                -name \
                -fo {output.extracted_seq}
        fi   
        """

rule prophage_evaluation:
    """
    Prophage region completeness evaluation using CheckV
    """
    input: 
        extracted_seq=os.path.join(RESULTS_DIR, "results", "prophage_region_sequences.fa")
    output:
        checkv_outdir=directory(os.path.join(RESULTS_DIR, "results", "checkv_evaluation")),
        summary=os.path.join(RESULTS_DIR, "results", "checkv_evaluation", "quality_summary.tsv")
    params:
        database=os.path.join(dir["db"], "checkv-db-v1.5")
    threads:
        config["resources"]["small_cpu"]
    resources:  
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "checkv.yml")
    shell:
        """
        if [[ -s {input.extracted_seq} ]]; then
            checkv end_to_end {input.extracted_seq} {output.checkv_outdir} -t {threads} -d {params.database}
        else
            echo "Sequence file is empty"
            touch {output.summary}
        fi   
        """

