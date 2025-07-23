import pandas as pd
from pathlib import Path


sample_name_constraint = "[A-Za-z][A-Za-z0-9]+"


wildcard_constraints:
    sample=sample_name_constraint,


if Path("samples.tsv").exists():
    SampleTable = pd.read_table("samples.tsv", index_col=0, sep="\t")
elif Path("samples.csv").exists():
    SampleTable = pd.read_csv("samples.csv", index_col=0)
else:
    raise FileNotFoundError(
        "No sample table found. Please provide a 'samples.tsv' or 'samples.csv' file."
    )


assert SampleTable.index.is_unique, "Sample table index is not unique"
assert SampleTable.index.str.match(
    sample_name_constraint
).all(), "Not all sample names correspond to sample name criteria"


if not SampleTable.columns.str.startswith("Reads_QC").any():
    import warnings

    warnings.warn(
        "QC-Fastq paths not in sample table, use default 'QC/reads/{sample}_{fraction}.fastq.gz'"
    )


SAMPLES = SampleTable.index.tolist()
PAIRED = SampleTable.columns.str.contains("R2").any()

if PAIRED:
    FRACTIONS = ["R1", "R2"]
else:
    FRACTIONS = ["se"]


logger.info(
    f"Found {len(SAMPLES)} {'paired-end' if PAIRED else 'single-end'} samples in the sample table."
)


def get_qc_reads(wildcards):
    headers = ["Reads_QC_" + f for f in FRACTIONS]

    try:
        return SampleTable.loc[wildcards.sample, headers]
    except KeyError:
        return expand("QC/reads/{{sample}}_{fraction}.fastq.gz", fraction=FRACTIONS)


def get_raw_reads(wildcards):
    headers = ["Reads_raw_" + f for f in FRACTIONS]
    fastq_dir = Path(config["fastq_dir"])

    return [fastq_dir / f for f in SampleTable.loc[wildcards.sample, headers]]
