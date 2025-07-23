


include: "sample_table.smk"


# HACK
# SAMPLES = ["C31519", "C530394", "s6pp", "S13Hp"]

KMC_sample_path = "Intermediate/kmc/samples/{sample}.k{k}"


rule all_kmc:
    input:
        expand(
            KMC_sample_path + ".{ext}",
            sample=SAMPLES,
            k=config["k_mer_size"],
            ext=["kmc_suf", "kmc_pre"],
        ),
        histograms=expand(
            "Intermediate/kmc/histograms/{sample}.k{k}.hist",
            sample=SAMPLES,
            k=config["k_mer_size"],
        ),


rule list_reads:
    input:
        get_qc_reads,
    output:
        temp("Intermediate/kmc/samples/{sample}.list"),
    run:
        with open(output[0], "w") as f:
            for fastq in list(input):
                f.write(f"{fastq}\n")


rule kmc_count:
    input:
        "Intermediate/kmc/samples/{sample}.list",
    output:
        expand("{path}.{ext}", path=KMC_sample_path, ext=["kmc_suf", "kmc_pre"]),
        execution_summary="logs/kmc/{sample}.{k}.json",
    params:
        ci=config["min_kmer_count"],
        cs=config["max_counter_value"],
        output=KMC_sample_path,
    threads: config["threads"]
    log:
        "logs/kmc/{sample}.{k}.log",
    conda:
        "../envs/kmc.yaml"
    resources:
        mem_mb=config["mem_default"],
        mem_gb=config["mem_default"] // 1000,
        time_min=60,
    benchmark:
        "logs/benchmarks/kmc/{sample}.{k}.txt"
    shell:
        "kmc"
        " -k{wildcards.k} -ci{params.ci} "
        "-t{threads} -m{resources.mem_gb} "
        "-j{output.execution_summary} "
        " --opt-out-size -hp "
        " @{input} {params.output} {resources.tmpdir} "
        "&> {log}"


"""
 -v - verbose mode (shows all parameter settings); default: false
  -k<len> - k-mer length (k from 1 to 256; default: 25)
  -m<size> - max amount of RAM in GB (from 1 to 1024); default: 12
  -sm - use strict memory mode (memory limit from -m<n> switch will not be exceeded)
  -hc - count homopolymer compressed k-mers (approximate and experimental)
  -p<par> - signature length (5, 6, 7, 8, 9, 10, 11); default: 9
  -f<a/q/m/bam/kmc> - input in FASTA format (-fa), FASTQ format (-fq), multi FASTA (-fm) or BAM (-fbam) or KMC (-fkmc); default: FASTQ
  -ci<value> - exclude k-mers occurring less than <value> times (default: 2)
  -cs<value> - maximal value of a counter (default: 255)
  -cx<value> - exclude k-mers occurring more of than <value> times (default: 1e9)
  -b - turn off transformation of k-mers into canonical form
  -r - turn on RAM-only mode
  -n<value> - number of bins
  -t<value> - total number of threads (default: no. of CPU cores)
  -sf<value> - number of FASTQ reading threads
  -sp<value> - number of splitting threads
  -sr<value> - number of threads for 2nd stage
  -j<file_name> - file name with execution summary in JSON format
  -w - without output
  -o<kmc/kff> - output in KMC of KFF format; default: KMC
  -hp - hide percentage progress (default: false)
  -e - only estimate histogram of k-mers occurrences instead of exact k-mer counting
  --opt-out-size - optimize output database size (may increase running time)
"""


localrules:
    all_kmc,
    list_reads,
    kmer_hist,


rule kmer_hist:
    input:
        rules.kmc_count.output,
    output:
        "Intermediate/kmc/histograms/{sample}.k{k}.hist",
    conda:
        "../envs/kmc.yaml"
    log:
        "logs/kmc/histograms/{sample}.{k}.log",
    params:
        input_path=KMC_sample_path,
    resources:
        mem_mb=config["mem_simple"],
    benchmark:
        "logs/benchmarks/kmc/histograms/{sample}.{k}.txt"
    shell:
        """
        kmc_tools transform {params.input_path} histogram {output} &> {log}
        """


'''
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
'''
