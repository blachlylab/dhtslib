dhtslib
=======


# Overview

D bindings and convenience wrappers for [htslib](https://github.com/samtools/htslib)

# Installation

Add `dhtslib` as a dependency to `dub.json`:

```
    "dependencies": {
        "dhtslib": "~>0.4.0",
```
(version number 0.4.0 is example but current as of today)

# Usage

## Convenience Wrappers

Object-oriented, idomatic(ish) D wrappers are available for:

* BGZF (`dhtslib.bgzf`)
* FASTA indexes (`dhtslib.faidx`)
* Tabix (`dhtslib.tabix`)

For example, this provides access to BGZF files by line as a consumable InputRange.

## htslib API

Direct bindings to htslib C API are available as submodules under `dhtslib.htslib`. For example, `import dhtslib.htslib.faidx` for direct access to the C function calls. The current compatible versions are 1.7-1.9.

Currently implemented fully or partially:

* bgzf
* faidx
* hts
* regidx
* sam
* tbx
