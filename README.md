dhtslib
=======


# Overview

D bindings and convenience wrappers for [htslib](https://github.com/samtools/htslib)

# Installation

Add `dhtslib` as a dependency to `dub.json`:

```
    "dependencies": {
        "dhtslib": "~>0.5.0",
```
(version number 0.5.0 is example but current as of today)

# Usage

## Convenience Wrappers

Object-oriented, idomatic(ish) D wrappers are available for:

* BGZF compressed files (`dhtslib.bgzf`)
* FASTA indexes (`dhtslib.faidx`)
* SAM/BAM/CRAM files and streams (`dhtslib.sam`)
* Tabix-indexed files (`dhtslib.tabix`)
* VCF/BCF files (`dhtslib.vcf`)

For example, this provides access to BGZF files by line as a consumable InputRange.
Or, for BAM files, the ability to query for a range (e.g. "chr1:1000000-2000000") and obtain an InputRange over the BAM records.

## htslib API

Direct bindings to htslib C API are available as submodules under `dhtslib.htslib`. 
Naming remains the same as the original `.h` include files.
For example, `import dhtslib.htslib.faidx` for direct access to the C function calls.
The current compatible versions are 1.7-1.9.

Currently implemented fully or partially:

* bgzf
* faidx
* hts
* hts_log
* kstring
* regidx
* sam
* tbx
* thread_pool (untested)
* vcf

Missing (so far):

* Some CRAM specific functions, although much CRAM functionality works with `sam_` functions
* hfile
* kbitset, kfunc, khash, klist, knetfile, kseq, ksort (mostly used internally anyway)
* synced_bcf_reader
* vcf_sweep
* vcfutils


# Bugs and Warnings

Do not call hts_log hts_log_* from a destructor, as it is potentialy allocating via `toStringz`