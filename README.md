dhtslib
=======

[![unittests](https://github.com/blachlylab/dhtslib/actions/workflows/unittests.yml/badge.svg)](https://github.com/blachlylab/dhtslib/actions/workflows/unittests.yml)
[![codecov.io](https://codecov.io/github/blachlylab/dhtslib/coverage.svg?branch=develop)](https://codecov.io/github/blachlylab/dhtslib?branch=develop)

# Overview

`dhtslib` provides D bindings, high-level abstractions, and additional functionality for [htslib](https://github.com/samtools/htslib), the most widely-used library for manipulation of high-throughput sequencing data. We currently support linux and OSX. Windows support is still in progress (see #38). More extensive documentation can be found at our [gitbook](https://blachlylab.gitbook.io/dhtslib/).

# Installation

Add `dhtslib` as a dependency to `dub.json`:

```
    "dependencies": {
        "dhtslib": "~>0.12.3+htslib-1.10",
    }
```
(version number 0.12.3 is example, `+htslib-1.10` represents the compatible htslib version; see https://dub.pm/package-format-json)

# Requirements

## htslib
A system installation of htslib >= v1.10 is required. You can find detailed install instructions [here](htslib.md).

# Usage
`dhtslib` usage information and examples can be found [here](usage.md).

## Dhtslib API (OOP Wrappers)

Object-oriented, idomatic D wrappers are available for:

* SAM/BAM/CRAM files and streams (`dhtslib.sam`)
* VCF/BCF files (`dhtslib.vcf`)
* BGZF compressed files (`dhtslib.bgzf`)
* FASTA indexes (`dhtslib.faidx`)
* Tabix-indexed files (`dhtslib.tabix`)

Additional functionality is provided for:

* GFF(2|3) files and streams (`dhtslib.gff`)
* BED files and streams (`dhtslib.bed`)
* FASTQ files and streams (`dhtslib.fastq`)
* Compile-time coordinate system (`dhtslib.coordinates`) to avoid off-by-one errors

All htslib bindings can be found under the `htslib` namespace (in prior versions they were under `dhtlsib.htslib`). These can be used directly as you would with `htslib`.


## htslib API

Direct bindings to htslib C API are available as submodules under the `htslib` namespace. Naming remains the same as the original `.h` include files. For example, `import htslib.faidx` for direct access to the C function calls. Where the OOP wrappers manage their own data along the the D garbage collector, these functions use traditional C memory management (or lack thereof). The current compatible htslib versions are 1.10+.

Currently implemented:

* bgzf
* cram (untested)
* faidx
* hfile
* hts\_endian
* hts\_expr (untested)
* hts\_log 
* hts\_os (untested)
* hts
* kbitset (untested)
* kfunc (untested)
* knetfile (untested)
* kroundup
* kstring
* regidx
* sam
* synced\_bcf\_reader (untested)
* tbx
* thread\_pool (untested)
* vcf
* vcf\_sweep (untested)
* vcfutils (untested)

Missing or work-in-progress:

* khash (see [dklib](https://github.com/blachlylab/dklib)), klist, kseq, ksort (mostly used internally anyway)

[dstep](https://github.com/jacob-carlborg/dstep) has matured and is an incredibly powerful tool for machine-assisted C-to-D translation. We've used dstep for the majority of bindings in the since version v0.11.0. After dstep translation, we port inline functions by hand as they are not translated, tweak some macros into templates (done although dstep already does an amazing job on simple `#define` macros translating to D templates!),  and update the documentation comments to ddoc format.

# FAQ

**Q**: Does this work with the latest htslib?

**A**:
Yes

**Q**: Why not use [bioD](https://github.com/biod/BioD)

**A**:
bioD, as a more general bioinformatics framework, is more comparable to bio-python, bio-ruby, bio-rust, etc.
bioD does have some excellent hts file format (BGZF and SAM) handling, and at one time sambamba, which relied on it, was faster than samtools.
However, the development resources poured into `htslib` overall are tremendous, and we with to leverage that rather than writing VCF, tabix, etc. code from scratch.

**Q**: How does this compare to bio-Rust's htslib bindings?

**A**: We love Rust, but dhtslib has way more complete bindings and more and better high level constructs :smile:. We have also implemented a novel compile-time type-safe coordinate system to mostly avoid off-by-one errors.

**Q**: Why am I getting a segfault?

**A**:
It's easy to get a segfault by using the direct C API incorrectly. Or possibly correctly. We have tried to eliminate most of this (use after free, etc.) in the OOP wrappers via refernece counting. If you are getting a segfault you cannot understand when using purely the high-level D API, please post an issue.


# Bugs and Warnings

Do not call `hts_log_*` from a destructor, as it is potentialy allocating via `toStringz`


# Programs made with dhtslib
1. [fade](https://github.com/blachlylab/fade): Fragmentase Artifact Detection and Elimination
2. [recontig](https://github.com/blachlylab/recontig): a program to convert different bioinformatics data types from one reference naming convention to another i.e UCSC to ensembl (chr1 to 1)

# Related projects

1. [gff3d](https://github.com/blachlylab/gff3d): GFF3 record reader/writer
2. [dklib](https://github.com/blachlylab/dklib): Templatized port of attractivechaos' klib, used extensively in htslib
