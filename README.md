dhtslib
=======


# Overview

D bindings and convenience wrappers for [htslib](https://github.com/samtools/htslib),
the most widely-used library for manipulation of high-throughput sequencing data.

# Installation

Add `dhtslib` as a dependency to `dub.json`:

```
    "dependencies": {
        "dhtslib": "~>0.10.0",
```
(version number 0.10.0 is example; see https://dub.pm/package-format-json)

# Requirements

## Dynamically linking to htslib (default)
A system installation of htslib v1.10 or 1.11 is required.

# Statically linking to htslib
`libhts.a` needs to be added to your project's source files.
Remember to link to all dynamic libraries configured when htslib was built. This may
include bz2, lzma, zlib, defalate, crypto, pthreads, curl.
Finally, if statically linking, the `-lhts` flag needs to be removed from compilation
by selecting the dub configuration `source-static` as the dub configuration type for dhtslib
within your own project's dub configuration file:

```
"subConfigurations": {
                "dhtslib": "source-static"
        },
```

# Usage

## D API (OOP Wrappers)

Object-oriented, idomatic D wrappers are available for:

* BGZF compressed files (`dhtslib.bgzf`)
* FASTA indexes (`dhtslib.faidx`)
* SAM/BAM/CRAM files and streams (`dhtslib.sam`)
* Tabix-indexed files (`dhtslib.tabix`)
* VCF/BCF files (`dhtslib.vcf`)

For example, this provides access to BGZF files by line as a consumable InputRange.
Or, for BAM files, the ability to query for a range (e.g. "chr1:1000000-2000000") and obtain an InputRange over the BAM records.
For most file type readers, indexing (`["coordinates"]`) queries return ranges of records. There are multiple options, including
`["chr1", 10_000_000 .. 20_000_000]` and `["chr1:10000000-20000000]`.
See the documentation for more details.

## htslib API

Direct bindings to htslib C API are available as submodules under `dhtslib.htslib`. 
Naming remains the same as the original `.h` include files.
For example, `import dhtslib.htslib.faidx` for direct access to the C function calls.
The current compatible versions are 1.10+

Currently implemented:

* bgzf
* faidx
* hts
* hts\_log
* kstring
* regidx
* sam
* tbx
* thread\_pool (untested)
* vcf

Missing or work-in-progress:

* Some CRAM specific functions, although much CRAM functionality works with `sam_` functions
* hfile
* kbitset, kfunc, khash, klist, knetfile, kseq, ksort (mostly used internally anyway)
* synced\_bcf\_reader
* vcf\_sweep
* vcfutils


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

**A**: We love Rust, but dhtslib has way more complete bindings and more and better high level constructs :smile:

**Q**: Why were htslib bindings ported by hand instead of using a C header/bindings translator as in hts-nim or rust-htslib?

**A**:
Whereas dstep and dpp are incredibly convenient for binding creation, we created these by hand from htslib `.h` files for several reasons.
First, this gave the authors of dhtslib a better familiarity with the htslib API including letting us get to know several lesser-known and internal functions.
Second, some elements (particuarlly `#define` macros) are difficult or impossible in some cases for machines to translate, or translate into efficient code; here we were sometimes able to replace these macros with smarter replacements than a simple macro-expansion-direct-translation. (**see 2020 update below -- dstep translates simple #defines into templates**)
Likewise, we were able to turn certain `#defines` and pseudo-generic functions into D templates, and to `pragma(inline, true)` them.
Finally, instead of dumping all the bindings into an interface file, we left the structure of the file intact to make it easier for the D developer to read the source file as the htslib authors intended the C headers to be read. In addition, this leaves docstring/documentation comments intact, whereas in other projects the direct API has no comments and the developer must refer to the C headers.

**(2020 UPDATE)** [dstep](https://github.com/jacob-carlborg/dstep) has matured and is an incredibly powerful tool for machine-assisted C-to-D translation. We've used dstep for the majority of bindings in the new `htslib-110` branch. After dstep translation, we still need to port inline functions by hand (done), tweak some macros into templates (done although dstep already does an amazing job on simple `#define` macros translating to D templates!), backport some fixes for Windows platforms and update the documentation comments to ddoc format.


**Q**: Why am I getting a segfault?

**A**:
It's easy to get a segfault by using the direct C API incorrectly. Or possibly correctly. We have tried to eliminate most of this (use after free, etc.) in the OOP wrappers. If you are getting a segfault you cannot understand when using purely the high-level D API, please post an issue.


# Bugs and Warnings

Zero-based versus one-based coordinates. Zero-based coordinates are used internally and also by the API for BCF/VCF and SAM/BAM types.
The `fadix` C API expects one-based coordinates; we have built this as a template for the user to specify the coordinate system.
See documentation for more details.

Do not call `hts_log_*` from a destructor, as it is potentialy allocating via `toStringz`


# See Also

1. [gff3d](https://github.com/blachlylab/gff3d) GFF3 record reader/writer
2. [dklib](https://github.com/blachlylab/dklib) Templatized port of attractivechaos' klib, used extensively in htslib
