#!/usr/bin/env python

import os
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat and check the contents of the samplesheet."

    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the structure specified in the provided CSV.
    """

    FQ_EXTENSIONS = (".fq.gz", ".fastq.gz")

    sample_mapping_dict = {}
    with open(file_in, "r") as fin:
        ## Check header
        MIN_COLS = 7
        HEADER = [
            "sample",
            "run_accession",
            "instrument_platform",
            "reads_type",
            "fastq_1",
            "fastq_2",
            "fasta",
            "kraken2_report",
            "kraken2_classifiedout",
            "centrifuge_out",
            "centrifuge_result",
            "diamond",
        ]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]

        ## Check for missing mandatory columns
        missing_columns = list(set(HEADER) - set(header))
        if len(missing_columns) > 0:
            print(
                "ERROR: Missing required column header -> {}. Note some columns can otherwise be empty.".format(
                    ",".join(missing_columns)
                )
            )
            sys.exit(1)

        ## Find locations of mandatory columns
        header_locs = {}
        for i in HEADER:
            header_locs[i] = header.index(i)

        ## Check sample entries
        for line in fin:
            ## Pull out only relevant columns for downstream checking
            line_parsed = [x.strip().strip('"') for x in line.strip().split(",")]

            # Check valid number of columns per row
            if len(line_parsed) < MIN_COLS:
                print_error(
                    f"Invalid number of columns (minimum = {MIN_COLS})!",
                    "Line",
                    line,
                )

            lspl = [line_parsed[i] for i in header_locs.values()]

            ## Check sample name entries

            (
                sample,
                run_accession,
                instrument_platform,
                reads_type,
                fastq_1,
                fastq_2,
                fasta,
                kraken2_report,
                kraken2_classifiedout,
                centrifuge_out,
                centrifuge_result,
                diamond,
            ) = lspl

            # Additional checks specific to your CSV format can be added here

            ## Create sample mapping dictionary = { sample: [ run_accession, instrument_platform, ... ] }
            if sample not in sample_mapping_dict:
                sample_mapping_dict[sample] = [
                    run_accession,
                    instrument_platform,
                    reads_type,
                    fastq_1,
                    fastq_2,
                    fasta,
                    kraken2_report,
                    kraken2_classifiedout,
                    centrifuge_out,
                    centrifuge_result,
                    diamond,
                ]
            else:
                print_error("Samplesheet contains duplicate rows!", "Line", line)

            # Check instrument_platform
            # (You can add more checks here based on your specific requirements)

    ## Write validated samplesheet with appropriate columns
    HEADER_OUT = [
        "sample",
        "run_accession",
        "instrument_platform",
        "reads_type",
        "fastq_1",
        "fastq_2",
        "fasta",
        "kraken2_report",
        "kraken2_classifiedout",
        "centrifuge_out",
        "centrifuge_result",
        "diamond",
    ]
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write("\t".join(HEADER_OUT) + "\n")
            for sample in sorted(sample_mapping_dict.keys()):
                fout.write(f"{sample}\t{','.join(sample_mapping_dict[sample])}\n")
    else:
        print_error("No entries to process!", "Samplesheet: {}".format(file_in))


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())