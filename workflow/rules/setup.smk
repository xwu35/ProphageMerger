# SET UP DIRECTORY
dir = dict()
dir["env"]     = os.path.join(workflow.basedir, "envs")
dir["scripts"] = os.path.join(workflow.basedir, "scripts")
dir["db"]      = os.path.join(workflow.basedir, "..", "db")

# DOWNLOAD DATABASES
rule download_virsorter2_db:
    output:
        os.path.join(dir["db"], "virsorter2_db", "Done_all_setup")
    params:
        outdir=os.path.join(dir["db"], "virsorter2_db")
    threads:
        config["resources"]["small_cpu"]
    resources:  
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "virsorter2.yml")
    shell:
        """
        virsorter setup -d {params.outdir} -j {threads}
        """

rule download_checkv_db:
    output:
        os.path.join(dir["db"], "checkv-db-v1.5", ".done")
    params:
        outdir=dir["db"]
    threads:
        config["resources"]["small_cpu"]
    resources:  
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "checkv.yml")
    shell:
        """
        checkv download_database {params.outdir} &&
        touch {output}
        """

rule download_genomad_db:
    output:
        os.path.join(dir["db"], "genomad_db", ".done")
    params:
        outdir=dir["db"]
    conda:
        os.path.join(dir["env"], "genomad.yml")
    shell:
        """
        genomad download-database {params.outdir}  &&
        touch {output}
        """

rule download_vibrant_db:
    output:
        os.path.join(dir["db"], "vibrant_db", ".done")
    params:
        outdir=os.path.join(dir["db"], "vibrant_db")
    conda:
        os.path.join(dir["env"], "vibrant.yml")
    shell:
        """
        download-db.sh {params.outdir} &&
        touch {output}
        """


