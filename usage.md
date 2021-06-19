# dhtslib usage and examples
### Introduction
### Coordinate system
Some of `dhtslib`'s novel functionality is pervasive throughout the library, mostly the compile-time, type-safe coordinate system. `dhtslib`, when dealing with integer-based position or coordinates, requires the use of `dhtslib.coordinates`. This system helps prevent off-by-one errors by asserting at compile-time that the coordinate system must be known for a pair of integer coordinates. To define a `Coordinate`:
```
import dhtslib.coordinates;

auto c1 = Coordinate!(Basis.zero)(0);
auto c2 = Coordinate!(Basis.zero)(1);
c2 = 2;
// c2 = 0; would result in an error as 
// a one-based system cannot have a coordinate zero
```
This defines a singular coordinate as zero or one-based. An easier way of specifing a coordinate is:
```
import dhtslib;

auto c1 = ZeroBased(0);
auto c2 = OneBased(1);
auto c3 = ZB(0);
auto c4 = OB(1);
```
A `Coordinate` can be converted from one `Basis` to another:
```
auto c1 = ZeroBased(0);
assert(c1.to!(Basis.one) == 1);
```

To specify a pair of `Coordinates`:
```
import dhtslib.coordinates;

auto c1 = Coordinates!(CoordSystem.zbho)(0, 1);
auto c2 = Coordinates!(CoordSystem.obho)(1, 2);
auto c3 = Coordinates!(CoordSystem.zbc)(0, 0);
auto c4 = Coordinates!(CoordSystem.obc)(1, 1);

assert(c1.size == 1);
assert(c2.size == 1);
assert(c3.size == 1);
assert(c4.size == 1);
```

This defines a singular coordinate pair as a coordinate system that combines basis and end. `Basis` being zero or one-based. `End` being open or closed (referred to as half-open because the starting coordinate is always closed). The availiable coordinate systems are: zero-based half-open, one-based half-open, zero-based closed, and one-based closed. An easier way of specifing a coordinate pair is:
```
import dhtslib;

auto c1 = ZeroBasedHalfOpen(0, 1);
auto c2 = OneBasedHalfOpen(1, 2);
auto c3 = ZeroBasedClosed(0, 0);
auto c4 = OneBasedClosed(1, 1);

auto c1 = ZBHO(0, 1);
auto c2 = OBHO(1, 2);
auto c3 = ZBC(0, 0);
auto c4 = OBC(1, 1);
```
`Coordinates` can be converted to a different coordinate system
```
auto c1 = ZeroBasedHalfOpen(0, 1);
auto c2 = OneBasedHalfOpen(1, 2);
assert(c1.to!(Coordsystem.obho) == c2);
```
`ChromCoordinates` allows us to specify a contig or chromosome with a pair
of coordinates:
```
auto c1 = ZBHO("chr1:0-1"); // returns a ChromCoordinates type
auto c2 = ChromCoordinates!(Coordsystem.zbho)("chr1",ZBHO(0, 1));
assert(c2.chrom == "chr1");
```

All of the readers, writers, and records for `dhtslib` will return either `Coordinates` or a `Coordinate` and will accept `Coordinates` instead of integer-based coordinates. 
```
// rec.coordinates will return the coordinates of the
// aligned portion of the first read in test.bam
SAMRecord rec = SAMReader("test.bam").allRecords.front;

// we have now filtered the BAM file to only the records that overlap
// the aligned region of the first read
auto filtered_recs = SAMReader("test.bam").query(rec.tid, rec.coordinates);
```
All functions that take coordinates are responsible for converting them to the correct coordsystem.
```
// rec.coordinates will return the coordinates of the
// gff3 record in test.gff3
// These coordinates are One-based closed
GFF3Record rec = GFF3Reader("test.gff3").front;

// we have now filtered the BAM file to only the records that overlap
// with the first region specified in the GFF3 file
// The one-based closed coordinates are automatically converted by
// SAMReader.query to zero-based half-open as is required by htslib
auto filtered_recs = SAMReader("test.bam").query(rec.contig, rec.coordinates);
```

### Readers
dhtslib provides readers for SAM/BAM(/CRAM untested), VCF/BCF, BED, GFF, FASTQ, faidx'd FASTA, generic BGZF or GZIP compressed files, and tabix indexed files. Readers automatically (via `htslib`) have support for reading compressed and remote (https or aws s3) files.
* `dhtslib.sam.reader` : `SAMReader` SAM/BAM(/CRAM untested))
* `dhtslib.vcf.reader` : `VCFReader` VCF/BCF
* `dhtslib.bed.reader` : `BedReader` 
* `dhtslib.gff.reader` : `GFFReader`, `GTFReader`, `GFF2Reader`,`GFF3Reader` 
* `dhtslib.fastq` : `FastqFile` 

The readers generally follow the following format:
* They are structs, but generally must be initialized
* They act as [InputRanges]() that returns the appropriate Record type
* They handle any availiable filtering
* They store headers as the appropriate Header type or as a string
* They own the data that backs the underlying `htslib` `htsFile`
    * They are reference counted
    * They control allocating and freeing that data

Exceptions:
* BGZFile acts as an InputRange via it's `byLine` and `byLineCopy` methods.
### Writers
dhtslib provides writers for SAM/BAM(/CRAM untested), VCF/BCF, BED, and GFF. 
* `dhtslib.sam.writer` : `SAMWriter` SAM/BAM(/CRAM untested))
* `dhtslib.vcf.writer` : `VCFWriter` VCF
* `dhtslib.bed.writer` : `BedWriter` 
* `dhtslib.gff.writer` : `GFF2Writer`,`GFF3Writer` 

The writers generally follow the following format:
* They are structs, but generally must be initialized
* They require the header upon initialization
* They have a write method that accepts their specific Record type
* They own the data that backs the underlying `htslib` `htsFile`
    * They are reference counted
    * They control allocating and freeing that data

Notes:
* SAM/BAM(/CRAM untested) writing accepts an enum that allows it to output SAM, compressed SAM, BAM, uncompressed BAM, and CRAM. By default it will try and deduce this based on the file extension of file it is writing to.
    * VCF/BCF writing currently does not have a similar feature but it will. 
* GFF, VCF, BED writing can only output uncompressed text.

### Records
dhtslib provides record types for SAM/BAM(/CRAM untested), VCF/BCF, BED, and GFF. 
* `dhtslib.sam.record` : `SAMReader` SAM/BAM(/CRAM untested))
* `dhtslib.vcf.record` : `VCFRecord` VCF/BCF
* `dhtslib.bed.record` : `BedRecord` 
* `dhtslib.gff.record` : `GFFRecord`, `GTFRecord`, `GFF2Record`,`GFF3Record` 
* `dhtslib.fastq` : `FastqRecord` 

The records generally follow the following format:
* They are structs, but generally must be initialized
* They can be built from scratch, though are usually generated from a reader
* Some require a header (VCF, SAM(optional)) upon initialization
* They own the data that backs the underlying `htslib` record type or string
    * They are reference counted 
    * They control allocating and freeing that data
    * They have helper methods for mutating the underlying data
* They allow access to the underlying `htslib` datatype pointers

## Examples
#### BAM/SAM manipulation
Loop over records in sam file and do something, then write to bam file. 
```
import dhtslib;

void main(string[] args){

    // open the sam file
    auto bamr = SAMReader("test.sam");
    // open the bam file
    auto bamw = SAMWriter("test.bam", bamr.header);
    foreach(SAMRecord rec; bamr.allRecords){
        // do something ...
        bamw.write(rec);
    }
}
```
#### BAM filtering
For each region in bed file, filter bam file and do something. Must have test.bam.bai (bam must be indexed).
```
import dhtslib;

void main(string[] args){

    // open the sam file
    auto bedr = BedReader("test.bed.gz");
    // open the bam file
    auto bamr = SAMReader("test.bam");
    foreach(auto bedrec; bedr)
    {
        // for each bed region in bed file
        // filter bam file and do something
        // must have test.bam.bai
                                // this is the same as bamr.query 
        foreach(SAMRecord rec; bamr[bedrec.contig, bedrec.coordinates]){
            // do something ...
        }
    }
}
```

#### VCF manipulation
Loop over records in compressed vcf file and do something, then write to vcf file.
```
import dhtslib;

void main(string[] args){

    // open the sam file
    auto vcfr = VCFReader("test.vcf.gz");
    // open the bam file
    auto vcfw = VCFWriter("test.vcf", vcfr.header);
    foreach(VCFRecord rec; vcfr){
        // do something ...
        vcfw.write(rec);
    }
}
```
#### VCF filtering
For each region in gff file, filter vcf file via tabix and do something. Must have test.vcf.gz.tbi (vcf must be bgzipped and tabix'd).
```
import dhtslib;

void main(string[] args){

    // open the sam file
    auto gffr = GFF3Reader("test.gff3.gz");
    // open the bam file
    auto vcfr = VCFReader("test.vcf.gz");
    foreach(auto bedrec; bedr)
    {

        // open tbi
        auto tbi = TabixIndexedFile("test.vcf.gz");
        // create region from bed record
        auto region = ZBHO(bedrec.contig, bedrec.coordinates);

        // filter VCF with region and tbi
        auto vcfr = VCFReader(tbi, region)

        // loop over records
        foreach(SAMRecord rec; vcfr){
            // do something ...
        }
    }
}
```