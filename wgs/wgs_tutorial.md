# WGS Sequins Tutorial

This is an example workflow for processing data that have had Sequins spiked
in. It highlights the places where Sequins specific steps must be taken, and it
is not intended as an example of a production workflow. You should adapt each
step to your own needs.

## Install dependencies

### Option 1 self-serve comprehensive

We provide a Docker container that has all dependencies pre-installed, this is
the recommended way to run this tutorial. You can download the latest version
of the image with:

```sh
docker pull ghcr.io/sequinsbio/wgs_tutorial:1.0.0
```

> The Docker container only supports x86_64 architectures. If you’re running
> this tutorial on an ARM64 architecture (e.g., Apple Silicon), you should set
> the DOCKER_DEFAULT_PLATFORM environment variable to linux/amd64 before
> running any Docker commands. This ensures that the container runs in an
> x86_64 emulation mode, which is necessary for compatibility with the tools
> included in the container.
>
> `export DOCKER_DEFAULT_PLATFORM=linux/amd64`

Alternatively, you can add the option `–platform linux/amd64` to every docker
command.

All subsequent commands in this tutorial can be run with the docker image.
Alternatively, if you are unable to run Docker or would prefer to run the tools
natively, you can install each dependency locally.

### Option 2 independent tools

This workflow uses a number of popular bioinformatics tools that need to be
installed.

- bwa (<https://github.com/lh3/bwa>)
- samtools (<https://github.com/samtools/samtools>)
- GATK (<https://github.com/broadinstitute/gatk>)
- RTG Tools (<https://github.com/RealTimeGenomics/rtg-tools>)

`sequintools` is the only Sequins specific tool needed for the workflow. There
are multiple ways you can install it, see details at
<https://github.com/sequinsbio/sequintools>.

## Running the workflow

The following sections step through the workflow in detail, so you can follow
along either running 1) inside the Docker container or 2) on your local
machine. You can start the Docker container with:

```sh
docker run -it --rm -v "$PWD":"$PWD" -w "$PWD" -u "$(id -u)":"$(id -g)" \
  ghcr.io/sequinsbio/wgs_tutorial:1.0.0
```

### 1. Obtain the Sequins resource bundle

The Sequins resource bundle contains Sequins specific files that are needed. To
obtain this resource, please complete the Request Access form
[here](https://sequinsdev.wpenginepowered.com/landing-wgs/).

The following files are provided:

- `sequin_sequences-mirror.fa`: FASTA file containing the sequence of every
  Sequin in the control set.
- `sequin_decoy.chrQ_mirror.fa`: The decoy chromosome that is concatenated to
  the standard reference genome.
- `sequin_regions.chrQ_mirror.bed`: A BED file containing the location of each
  Sequin in the decoy chromosome.
- `sequin_variants.chrQ_mirror.vcf.gz`: A VCF file detailing the variants in
  the Sequins in the decoy chromosome coordinates.
- `sequin_regions.hg38.bed`: A BED file containing the locations in GRCh38 that
  each Sequin is based on.
- `sequin_variants.hg38.vcf.gz`: A VCF file detailing the variants in the
  Sequins converted to their GRCh38 coordinates.

### 2. Build a Sequins augmented reference genome

The first step is to concatenate the decoy chromosome to your reference genome.
Here we are downloading a copy of GRCh38, but you can use your own reference
genome if you prefer.

This step only needs to be done once, and the resulting FASTA file can be used
for all subsequent analyses.

```sh
pushd wgs_core_control_set_v1.1
curl -sSLf -O https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz
gzip -d GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz
cat \
  GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  sequin_decoy.chrQ_mirror.fa \
  >grch38_with_sequins.fasta
samtools faidx grch38_with_sequins.fasta
samtools dict grch38_with_sequins.fasta >grch38_with_sequins.dict
bwa index grch38_with_sequins.fasta
popd
```

### 3. Download the tutorial example data

We provide some example data you can use to test the Sequins workflow. It’s
taken from whole genome sequencing of the reference sample HG002, and has been
subsetted to reduce its size. To obtain this data, please complete the Request
Access form [here](https://sequinsdev.wpenginepowered.com/landing-wgs/).

> The subsetted regions of the human genome are regions that the corresponding
> Sequins were designed to control for.

The following steps will walk you through a basic NGS workflow.

### 4. Align data

We begin by aligning our data with bwa to the Sequins augmented reference
genome we created earlier.

```sh
FASTA=wgs_core_control_set_v1.1/grch38_with_sequins.fasta
R1=wgs_tutorial_data_v1/wgs_tutorial_data_v1_R1.fastq.gz
R2=wgs_tutorial_data_v1/wgs_tutorial_data_v1_R2.fastq.gz
READGROUP="@RG\\tID:A\\tSM:A\\tLB:libA\\tPU:puA\\tPL:ILLUMINA"
bwa mem -M -t 4 -R "$READGROUP" "$FASTA" "$R1" "$R2" |
  samtools fixmate -u -m - - |
  samtools sort -u - |
  samtools markdup -u - - |
  samtools view -b -o wgs_tutorial_data_v1.bam --write-index
```

### 5) Calibrate the Sequins in the aligned BAM

Before calling variants, we need to calibrate the Sequins in the aligned BAM.

We do this to make comparisons between the Sequins and the sample data more
comparable – Sequins will typically be sequenced to a much higher coverage than
the WGS sample, making it more likely that variants can be called. The
calibration step reduces the coverage of each Sequin region to match the
coverage in the sample at the corresponding location. This gives a more
realistic estimation of the true power to call variants in the sample.

```sh
sequintools calibrate \
  --sample-bed wgs_core_control_set_v1.1/sequin_regions.hg38.bed \
  --bed wgs_core_control_set_v1.1/sequin_regions.chrQ_mirror.bed \
  --write-index \
  -o wgs_tutorial_data_v1.calibrated.bam \ wgs_tutorial_data_v1.bam
```

### 6) Call variants

We can now call variants in the calibrated BAM file using GATK’s
HaplotypeCaller.

```sh
gatk HaplotypeCaller \
  -R "$FASTA" \
  -I wgs_tutorial_data_v2.calibrated.bam \
  -O wgs_tutorial_data_v2.calibrated.vcf.gz \
  -L wgs_core_control_set_v1.1/sequin_regions.chrQ_mirror.bed
```

Sequins are compatible with any variant caller, so you can use your preferred
caller at this step; however, your results may vary depending on the caller you
choose.

### 7. Evaluate the results

To evaluate the variants with RTG Tools, we first need to create the required
RTG Sequence Data File (SDF) from the reference FASTA:

```sh
pushd wgs_core_control_set_v1.1
rtg format -o grch38_with_sequins.sdf grch38_with_sequins.fasta
popd
```

We can then perform the evaluation with `vcfeval` command of RTG Tools. To
simplify this tutorial, we will only be evaluating SNVs in this tutorial and
NOT SVs.

First we generate a BED file with the Sequins that represent small variants:

```sh
grep -f wgs_tutorial_data_v1/sv_groups.txt \
  -v wgs_core_control_set_v1.1/sequin_regions.chrQ_mirror.bed \
  >sequin_regions.chrQ_mirror.no_sv.bed
```

then run the evaluation:

```sh
rtg vcfeval \
  -t wgs_core_control_set_v1.1/grch38_with_sequins.sdf \
  -b wgs_core_control_set_v1.1/sequin_variants.chrQ_mirror.vcf.gz \
  -c wgs_tutorial_data_v1.calibrated.vcf.gz \
  -e sequin_regions.chrQ_mirror.no_sv.bed \
  --decompose \
  --no-roc \
  -o vcfeval_calibrated
```

| True Pos | False Pos | False Neg | Precision | Sensitivity | F-measure |
| -------- | --------- | --------- | --------- | ----------- | --------- |
| 62       | 0         | 0         | 1.0000    | 1.0000      | 1.0000    |

Here we can see that Sequins have proven that the sequencing experiment was
successful, and the data is of good quality. If, however, we start to see false
negatives or false positives, we can investigate which types of Sequins are
failing: low/high GC, repeats etc. and make informed decisions about how to
proceed with reporting variants in the sample.
