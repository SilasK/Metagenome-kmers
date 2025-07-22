


include: "sample_table.smk"


rule all_kmc:
    input:
        expand(
            "Intermediate/kmc/samples/{sample}.k{k}",
            sample=SAMPLES,
            k=config["k_mer_size"],
        ),
        histograms=expand(
            "Intermediate/kmc/histograms/{sample}.k{k}.hist",
            sample=SAMPLES,
            k=config["k_mer_size"],
        ),


rule kmc_count:
    input:
        get_qc_reads,
    output:
        db="Intermediate/kmc/samples/{sample}.k{k}",
    params:
        ci=config["min_kmer_count"],
    threads: config["threads"]
    log:
        "logs/kmc/{sample}.{k}.log",
    conda:
        "../envs/kmc.yaml"
    resources:
        mem_mb=config["kmc_mem_mb"],
        mem_gb=config["kmc_mem_mb"] / 1024,
    shell:
        "kmc"
        " -k{wildcards.k} -ci{params.ci} "
        "-t{params.threads} -m{resources.mem_gb} "
        "-fqgz {input} "
        "{output} {resources.tmpdir}"
        " -v &> {log}"


rule kmer_hist:
    input:
        db="Intermediate/kmc/samples/{sample}.k{k}",
    output:
        "Intermediate/kmc/histograms/{sample}.k{k}.hist",
    conda:
        "../envs/kmc.yaml"
    shell:
        """
        kmc_tools hist {input.db} {output}
        """


rule kmc_intersect:
    input:
        db1="Intermediate/kmc/samples/{sample1}.k{k}",
        db2="Intermediate/kmc/samples/{sample2}.k{k}",
    output:
        "Intermediate/kmc/intersections/{sample1}_vs_{sample2}.k{k}",
    conda:
        "../envs/kmc.yaml"
    shell:
        """
        kmc_tools intersect {input.db1} {input.db2} {output.db}
        """


rule kmc_dump:
    input:
        db=f"{SHARED_DIR}/{{sample1}}_vs_{{sample2}}.kmc",
    output:
        txt=f"{SHARED_DIR}/{{sample1}}_vs_{{sample2}}.txt",
    conda:
        "../envs/kmc.yaml"
    shell:
        """
        kmc_tools dump {input.db} {output.txt}
        """
