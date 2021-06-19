For example, this provides access to BGZF files by line as a consumable InputRange.
Or, for BAM files, the ability to query for a range (e.g. "chr1:1000000-2000000") and obtain an InputRange over the BAM records.
For most file type readers, indexing (`["coordinates"]`) queries return ranges of records. There are multiple options, including
`["chr1", 10_000_000 .. 20_000_000]` and `["chr1:10000000-20000000]`.
See the documentation for more details.