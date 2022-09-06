# gff_compaction

This script is designed to merge all overlapping ranges in a gff file into single records.

To do this, the gff file is transformed into a GenomicRanges object, and split in subset by each unique element in SequenceID column and by each unique element in Attributes column (specified by "Query="). Then, the records of each subset are sorted and all overlapping ranges are reduce to a single record. Finally, the list of subsets are transformed back into a gff format.

## Usage

gff_compaction.R [options]

## Options

| Options | Description |
| --- | --- |
| -f --file | GFF file |
| -s --source | Method used to obtain the GFF file |
| -h --help | Show this help message and exit |
