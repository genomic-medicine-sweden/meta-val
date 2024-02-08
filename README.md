<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-metaval_logo_dark.png">
    <img alt="nf-core/metaval" src="docs/images/nf-core-metaval_logo_light.png">
  </picture>
</h1>

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/metaval)

## Introduction

**nf-core/metaval** is a bioinformatics pipeline that verifies the organisms predicted by the nf-core/taxprofiler pipeline using metagenomic data, including both Illumina short-gun sequencing and Nanopore sequencing data.

At moment, meta-val only checks the classification results from three classifiers `Kraken2`, `Centrifuge` and `diamond`.

## Pipeline summary

1. Extract classified reads for organisms of interest, such as all identified viruses or a predefined list of organisms.

2. Use `BLAST` to identify the closet reference genome for the extracted reads.

3. Map the extracted reads to reference genomes using `Bowtie2` for Illumina reads and `minimap2` for Nanopore reads.

4. Construct consensus maps for the mapped reads.

5. Generate Coverage plots

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,run_accession,instrument_platform,fastq_1,fastq_2,kraken2_report,kraken2_classifiedout,centrifuge_out,centrifuge_result,diamond
sample1,run1,ILLUMINA,sample1.unmapped_1.fastq.gz,sample1.unmapped_2.fastq.gz,sample1.kraken2.kraken2.report.txt,sample1.kraken2.kraken2.classifiedreads.txt,sample1.centrifuge.txt,sample1.centrifuge.results.txt,sample1.diamond.tsv
sample2,run1,ILLUMINA,sample2.unmapped_1.fastq.gz,sample2.unmapped_2.fastq.gz,sample2.kraken2.kraken2.report.txt,sample2.kraken2.kraken2.classifiedreads.txt,sample2.centrifuge.txt,sample2.centrifuge.results.txt,sample2.diamond.tsv
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

Now, you can run the pipeline using:

```bash
nextflow run main.nf \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
   --run_kraken2 --run_centrifuge --run_diamond
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/metaval/usage) and the [parameter documentation](https://nf-co.re/metaval/parameters).

## Test data

There are three test datasets within `assets/data/`, produced by the `nf-core/taxprofiler` pipeline

- `taxprofiler_test_data`: produced by running the `test.config`
- `taxprofiler_test_full_data`: produced by running the `test_full.config`
- `test_data_version2_subset`: produced by running the data downloaded from https://www.nature.com/articles/s41598-021-83812-x

The corresponding input samplesheets are stored in `assets/`

- `samplesheet_v1.csv`:results of taxprofiler test data; no viruses; single-end (`perform_runmerging`).
- `samplesheet_v2.csv`:results of taxprofiler full test data; no viruses; single-end (`perform_runmerging`).
- `samplesheet_v3.csv`: with viruses; subset data from `test_data_version2_subset` (sample 20% of pair-end reads).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/metaval/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/metaval/output).

## Credits

nf-core/metaval was originally written by LilyAnderssonLee.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#metaval` channel](https://nfcore.slack.com/channels/metaval) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/metaval for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
