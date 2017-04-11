# KAD - K-mer Abundance Database

KAD, is a database build with leveldb to store k-mer abundances of RNA-Seq
samples.

## Installation

KAD depends on [rocksdb]('https://github.com/facebook/rocksdb').

To compile KAD go to `src/` directory and run `make`.

## Usage

By default the kad database is located on the working directy under the ".kad" directory.

Use `kad index [sample_name] counts.tsv` to index the counts from one sample. The counts should be formated as a tabulated file with the kmer sequence in the first column and the count in the second.
