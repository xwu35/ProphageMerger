rule virsorter2:
    """
    run VirSorter2 using default parameters
    """
    input: 
        seq=SEQUENCE,
        done=os.path.join(dir["db"], "virsorter2_db", "Done_all_setup")
    output:
        virsorter2_outdir=directory(os.path.join(RESULTS_DIR, "predictions", "virsorter2_output")),
        fasta_file=os.path.join(RESULTS_DIR, "predictions", "virsorter2_output", "final-viral-combined.fa"),
        boundary=os.path.join(RESULTS_DIR, "predictions", "virsorter2_output", "final-viral-boundary.tsv")
    threads:
        config["resources"]["med_cpu"]
    resources:  
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "virsorter2.yml")
    shell:
        """
        if [[ -s {input.seq} ]]; then
            virsorter run -i {input.seq} -w {output.virsorter2_outdir} -j {threads} all
        else
            echo "Sequence file is empty" 
            touch {output.fasta_file}
            touch {output.boundary}
        fi   
        """

rule checkv_after_virsorter2:
    """
    run CheckV after VirSorter2
    """
    input: 
        fasta_file=os.path.join(RESULTS_DIR, "predictions", "virsorter2_output", "final-viral-combined.fa"),
        done=os.path.join(dir["db"], "checkv-db-v1.5", ".done")
    output:
        checkv_outdir=directory(os.path.join(RESULTS_DIR, "predictions", "checkv_virsorter2")),
        provirus=os.path.join(RESULTS_DIR, "predictions", "checkv_virsorter2", "proviruses.fna")
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
        if [[ -s {input.fasta_file} ]]; then
            checkv end_to_end {input.fasta_file} {output.checkv_outdir} -t {threads} -d {params.database}
        else
            echo "Sequence file is empty"
            touch {output.provirus}
        fi   
        """

rule genomad:
    """
    Prophage prediction using geNomad's find_provirus module
    """
    input: 
        seq=SEQUENCE,
        done=os.path.join(dir["db"], "genomad_db", ".done")
    output:
        genomad_outdir=directory(os.path.join(RESULTS_DIR, "predictions", "genomad_output")),
        provirus=os.path.join(RESULTS_DIR, "predictions", "genomad_output", SEQ_NAME + "_find_proviruses", SEQ_NAME + "_provirus.tsv")
    params:
        database=os.path.join(dir["db"], "genomad_db"),
        sensitivity=config["genomad"]["sensitivity"]
    threads:
        config["resources"]["big_cpu"]
    resources:  
        mem_mb=config["resources"]["big_mem"]
    conda:
        os.path.join(dir["env"], "genomad.yml")
    shell:
        """
        if [[ -s {input.seq} ]]; then
            genomad annotate {input.seq} {output.genomad_outdir} {params.database} --sensitivity {params.sensitivity} -t {threads}
            genomad find-proviruses {input.seq} {output.genomad_outdir} {params.database} --sensitivity {params.sensitivity} -t {threads}
        else
            echo "Sequence file is empty"
            touch {output.provirus} 
        fi   
        """

rule vibrant:
    """
    Prophage prediction using Vibrant with default parameters
    """
    input: 
        seq=SEQUENCE,
        done=os.path.join(dir["db"], "vibrant_db", ".done")
    output:
        vibrant_outdir=directory(os.path.join(RESULTS_DIR, "predictions", "vibrant_output")),
        done=os.path.join(RESULTS_DIR, "predictions", "vibrant_output", ".done")
    params:
        database=os.path.join(dir["db"], "vibrant_db", "databases"),
        prophage=os.path.join(RESULTS_DIR, "predictions", "vibrant_output", "VIBRANT_" + SEQ_NAME, "VIBRANT_results_" + SEQ_NAME, "VIBRANT_integrated_prophage_coordinates_" + SEQ_NAME + ".tsv")
    threads:
        config["resources"]["med_cpu"]
    resources:  
        mem_mb=config["resources"]["med_mem"]
    conda:
        os.path.join(dir["env"], "vibrant.yml")
    shell:
        """
        if [[ -s {input.seq} ]]; then
            VIBRANT_run.py -i {input.seq} -folder {output.vibrant_outdir} -d {params.database} -t {threads} &&
            touch {output.done}
        else
            echo "Sequence file is empty"
            touch {params.prophage} 
            touch {output.done}
        fi   
        """

rule combine_prophage_coordinates:
    """
    Combine all predicted prophage coordinates
    """
    input: 
        vir2_boundary=os.path.join(RESULTS_DIR, "predictions", "virsorter2_output", "final-viral-boundary.tsv"),
        checkv_provirus=os.path.join(RESULTS_DIR, "predictions", "checkv_virsorter2", "proviruses.fna"),
        genomad_provirus=os.path.join(RESULTS_DIR, "predictions", "genomad_output", SEQ_NAME + "_find_proviruses", SEQ_NAME + "_provirus.tsv"),
        # for tracking
        done=os.path.join(RESULTS_DIR, "predictions", "vibrant_output", ".done")
    output:
        all_coordinates=os.path.join(RESULTS_DIR, "results", "all_coordinates.txt"),
        final_coordinates=os.path.join(RESULTS_DIR, "results", "final_coordinates.txt"),
        bed_file=os.path.join(RESULTS_DIR, "results", "final_coordinates_0-based.bed")
    params:
        script=os.path.join(dir["scripts"], "combine_prophage_coordinates.py"),
        vibrant_prophage=os.path.join(RESULTS_DIR, "predictions", "vibrant_output", "VIBRANT_" + SEQ_NAME, "VIBRANT_results_" + SEQ_NAME, "VIBRANT_integrated_prophage_coordinates_" + SEQ_NAME + ".tsv")
    shell:
        """
        python {params.script} \
            --virsorter2 {input.vir2_boundary} \
            --checkv {input.checkv_provirus} \
            --genomad {input.genomad_provirus} \
            --vibrant {params.vibrant_prophage} \
            --all_coordinates {output.all_coordinates} \
            --final_coordinates {output.final_coordinates} \
            --bed_file {output.bed_file}
        """