rule extract_chromosome_seq_for_prediction:
    input: 
        lambda wildcards: SEQUENCE_MAP[wildcards.genome]
    output:
        seq=os.path.join(RESULTS_DIR, "extracted_chromosome", "{genome}.fa")
    params:
        script=os.path.join(dir["scripts"], "extract_longest_sequence.py")
    shell:
        """
        {params.script} -i {input} -n {wildcards.genome} -o {output.seq}
        """

rule get_chromosome_length:
    """
    Get the length of chromosome sequences
    """
    input:
        expand(os.path.join(RESULTS_DIR, "extracted_chromosome", "{genome}.fa"), genome=GENOME)
    output:
        os.path.join(RESULTS_DIR, "extracted_chromosome", "chromosome_length.txt")
    conda:
        os.path.join(dir["env"], "seqkit.yml")
    shell:
        """
        # DON'T NEED TO DEAL WITH EMPTY FASTA, seqkit WILL JUST OUTPUT HEADER
        seqkit fx2tab --length --name --header-line {input} > {output}
        """

rule virsorter2:
    """
    Prophage prediction using VirSorter2 with default parameters
    """
    input: 
        seq=os.path.join(RESULTS_DIR, "extracted_chromosome", "{genome}.fa"),
        done=os.path.join(dir["db"], "virsorter2_db", "Done_all_setup")
    output:
        virsorter2_outdir=directory(os.path.join(RESULTS_DIR, "predictions", "{genome}", "virsorter2_output")),
        fasta_file=os.path.join(RESULTS_DIR, "predictions", "{genome}", "virsorter2_output", "final-viral-combined.fa"),
        boundary=os.path.join(RESULTS_DIR, "predictions", "{genome}", "virsorter2_output", "final-viral-boundary.tsv")
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
        fasta_file=os.path.join(RESULTS_DIR, "predictions", "{genome}", "virsorter2_output", "final-viral-combined.fa"),
        done=os.path.join(dir["db"], "checkv-db-v1.5", ".done")
    output:
        checkv_outdir=directory(os.path.join(RESULTS_DIR, "predictions", "{genome}", "checkv_virsorter2")),
        provirus=os.path.join(RESULTS_DIR, "predictions", "{genome}", "checkv_virsorter2", "proviruses.fna")
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
        seq=os.path.join(RESULTS_DIR, "extracted_chromosome", "{genome}.fa"),
        done=os.path.join(dir["db"], "genomad_db", ".done")
    output:
        genomad_outdir=directory(os.path.join(RESULTS_DIR, "predictions", "{genome}", "genomad_output")),
        provirus=os.path.join(RESULTS_DIR, "predictions", "{genome}", "genomad_output", "{genome}_find_proviruses", "{genome}_provirus.tsv")
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
        seq=os.path.join(RESULTS_DIR, "extracted_chromosome", "{genome}.fa"),
        done=os.path.join(dir["db"], "vibrant_db", ".done")
    output:
        vibrant_outdir=directory(os.path.join(RESULTS_DIR, "predictions", "{genome}", "vibrant_output")),
        done=os.path.join(RESULTS_DIR, "predictions", "{genome}", "vibrant_output", ".done")
    params:
        database=os.path.join(dir["db"], "vibrant_db", "databases"),
        prophage=os.path.join(RESULTS_DIR, "predictions", "{genome}", "vibrant_output", "VIBRANT_{genome}", "VIBRANT_results_{genome}", "VIBRANT_integrated_prophage_coordinates_{genome}.tsv")
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

rule cenote_taker3:
    """
    Prophage prediction using Cenote-Taker3. Dantas lab used this one in their paper, but with v3.2.1
    """
    input: 
        seq=os.path.join(RESULTS_DIR, "extracted_chromosome", "{genome}.fa"),
        db_done=os.path.join(dir["db"], "ct3_db", ".done")
    output:
        # cenote-taker3 doesn't output summary files if no viruses were found, so cannot track summary files
        done=os.path.join(RESULTS_DIR, "predictions", "{genome}", "cenote_taker3", ".done")
    params:
        tmp_dir=directory(os.path.join(RESULTS_DIR, "predictions", "{genome}", "cenote_taker3_tmp")),
        dst_dir=directory(os.path.join(RESULTS_DIR, "predictions", "{genome}", "cenote_taker3")),
        database=os.path.join(dir["db"], "ct3_db"),
        settings=config["cenote_taker3"]["settings"] 
    threads:
        config["resources"]["med_cpu"]
    resources:  
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "cenote.yml")
    shell:
        """
        if [[ -s {input.seq} ]]; then
            cenotetaker3 -c {input.seq} \
                -r cenote_taker3 \
                -t {threads} \
                -wd {params.tmp_dir} \
                --cenote-dbs {params.database} \
                {params.settings}  && 
            mv {params.tmp_dir}/cenote_taker3/* {params.dst_dir} &&
            rm -r {params.tmp_dir} &&
            touch {output.done}
        else
            echo "Sequence file is empty"
            touch {output.done}
        fi   
        """

rule prokka:
    """
    Prokka annotation to produce the .gbk for PhiSpy
    """
    input: 
        seq=os.path.join(RESULTS_DIR, "extracted_chromosome", "{genome}.fa")
    output:
        gbk=os.path.join(RESULTS_DIR, "predictions", "{genome}", "prokka_output", "{genome}.gbk")
    params:
        dir=directory(os.path.join(RESULTS_DIR, "predictions", "{genome}", "prokka_output"))
    threads:
        config["resources"]["med_cpu"]
    resources:  
        mem_mb=config["resources"]["med_mem"]
    conda:
        os.path.join(dir["env"], "prokka.yml")
    shell:
        """
        if [[ -s {input.seq} ]]; then
            prokka \
                --outdir {params.dir} \
                --prefix {wildcards.genome} \
                --cpus {threads} --force \
                {input.seq}
        else
            echo "Sequence file is empty"
            touch {output.gbk}
        fi   
        """

rule phispy:
    """
    Prophage prediction using Vibrant with default parameters
    """
    input: 
        gbk=os.path.join(RESULTS_DIR, "predictions", "{genome}", "prokka_output", "{genome}.gbk")
    output:
        phispy=os.path.join(RESULTS_DIR, "predictions", "{genome}", "phispy_output", "prophage_coordinates.tsv")
    params:
        dir=directory(os.path.join(RESULTS_DIR, "predictions", "{genome}", "phispy_output"))
    threads:
        config["resources"]["small_cpu"]
    resources:  
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "phispy.yml")
    shell:
        """
        if [[ -s {input.gbk} ]]; then
            PhiSpy.py \
                {input.gbk} \
                -o {params.dir} \
                --threads {threads}
        else
            echo "Sequence file is empty"
            touch {output.phispy} 
        fi   
        """

rule combine_prophage_coordinates:
    """
    Combine all predicted prophage coordinates
    """
    input: 
        vir2_boundary=os.path.join(RESULTS_DIR, "predictions", "{genome}", "virsorter2_output", "final-viral-boundary.tsv"),
        checkv_provirus=os.path.join(RESULTS_DIR, "predictions", "{genome}", "checkv_virsorter2", "proviruses.fna"),
        genomad_provirus=os.path.join(RESULTS_DIR, "predictions", "{genome}", "genomad_output", "{genome}_find_proviruses", "{genome}_provirus.tsv"),
        phispy=os.path.join(RESULTS_DIR, "predictions", "{genome}", "phispy_output", "prophage_coordinates.tsv"),
        # for tracking
        vibrant_done=os.path.join(RESULTS_DIR, "predictions", "{genome}", "vibrant_output", ".done"),
        ct3_done=os.path.join(RESULTS_DIR, "predictions", "{genome}", "cenote_taker3", ".done")
    output:
        all_coordinates=os.path.join(RESULTS_DIR, "results", "{genome}_all_coordinates.txt"),
        final_coordinates=os.path.join(RESULTS_DIR, "results", "{genome}_final_coordinates.txt"),
        bed_file=os.path.join(RESULTS_DIR, "results", "{genome}_final_coordinates_0-based.bed")
    params:
        script=os.path.join(dir["scripts"], "combine_prophage_coordinates.py"),
        vibrant_prophage=os.path.join(RESULTS_DIR, "predictions", "{genome}", "vibrant_output", "VIBRANT_{genome}", "VIBRANT_results_{genome}", "VIBRANT_integrated_prophage_coordinates_{genome}.tsv"),
        cenote_prune=os.path.join(RESULTS_DIR, "predictions", "{genome}", "cenote_taker3", "cenote_taker3_prune_summary.tsv"),
        cenote_virus=os.path.join(RESULTS_DIR, "predictions", "{genome}", "cenote_taker3", "cenote_taker3_virus_summary.tsv")
    shell:
        """
        # if cenote_prune file exists
        if [[ -f {params.cenote_prune} ]]; then
            python {params.script} \
                --virsorter2 {input.vir2_boundary} \
                --checkv {input.checkv_provirus} \
                --genomad {input.genomad_provirus} \
                --vibrant {params.vibrant_prophage} \
                --phispy {input.phispy} \
                --cenote_prune {params.cenote_prune} \
                --cenote_virus {params.cenote_virus} \
                --all_coordinates {output.all_coordinates} \
                --final_coordinates {output.final_coordinates} \
                --bed_file {output.bed_file}
        else
            python {params.script} \
                --virsorter2 {input.vir2_boundary} \
                --checkv {input.checkv_provirus} \
                --genomad {input.genomad_provirus} \
                --vibrant {params.vibrant_prophage} \
                --phispy {input.phispy} \
                --all_coordinates {output.all_coordinates} \
                --final_coordinates {output.final_coordinates} \
                --bed_file {output.bed_file}
        fi   
        """

rule include_nonprophage_coordinates:
    """
    use gaps between prophage region as nonprophage coordinates for plotting
    """
    input:
        final_coordinates=os.path.join(RESULTS_DIR, "results", "{genome}_final_coordinates.txt"),
        genome_length=os.path.join(RESULTS_DIR, "extracted_chromosome", "chromosome_length.txt")
    output:
        os.path.join(RESULTS_DIR, "results", "{genome}_whole_genome_coordinates.txt"),
    params:
        script=os.path.join(dir["scripts"], "fill_in_nonprophage_coordinates.py")
    shell:
        """
        if grep -q "NO PROPHAGE REGIONS FOUND" {input.final_coordinates}; then
            cp {input.final_coordinates} {output}
        else
            python {params.script} \
                -i {input.final_coordinates} \
                -g {input.genome_length} \
                -o {output} 
        fi   
        """