module dhtslib.bed;

import std.range : inputRangeObject, InputRangeObject;
import std.algorithm.iteration: splitter;
import std.range: drop, enumerate;
import std.conv: to;

import dhtslib.coordinates;
import dhtslib.bgzf;
import dhtslib.tabix;

struct RGB
{
    ubyte red;
    ubyte green;
    ubyte blue;
}

/// Represents a record in a bed file.
/// Based on UCSC format and methods are derived from
/// UCSC's specifications.
struct BedRecord
{
    string line;
    
    /// column 1: The name of the chromosome or scaffold.
    @property contig() const { return this.line.splitter('\t').front; }

    /// Columns 2 & 3 as coordinate pair, Zero-based half-open.
    /// column 2: The starting position of the feature in the chromosome or scaffold.
    /// column 3: The ending position of the feature in the chromosome or scaffold. 
    @property coordinates() const
    {
        auto start = (cast(string)this.line.splitter('\t').drop(1).front).to!long;
        auto end = (cast(string)this.line.splitter('\t').drop(2).front).to!long;
        return ZBHO(start, end);
    }
    
    /// column 2: The starting position of the feature in the chromosome or scaffold.
    @property start() const { return this.coordinates.start; }
    
    /// column 3: The ending position of the feature in the chromosome or scaffold. 
    @property end() const { return this.coordinates.end; }

    /// column 4: Defines the name of the BED line.
    @property name() const { return this.line.splitter('\t').drop(3).front; }

    /// column 5: A score between 0 and 1000.
    @property score() const { return this.line.splitter('\t').drop(4).front.to!int; }
    
    /// column 6: Defines the strand. Either "." (=no strand) or "+" or "-".
    @property strand() const {
        return cast(char)this.line.splitter('\t').drop(5).front[0];
    }

    /// column 7: The starting position at which the feature is drawn thickly;
    @property thickStart() const { return this.line.splitter('\t').drop(6).front.to!int; }

    /// column 8: The ending position at which the feature is drawn thickly
    @property thickEnd() const { return this.line.splitter('\t').drop(7).front.to!int; }

    /// column 9: An RGB value of the form R,G,B (e.g. 255,0,0).
    @property itemRGB() const 
    {
        auto str = this.line.splitter('\t').drop(8).front;
        if(str == "0") return RGB.init;
        return RGB(
            str.splitter(',').front.to!ubyte,
            str.splitter(',').drop(1).front.to!ubyte,
            str.splitter(',').drop(2).front.to!ubyte
            );
    }

    /// column 10: The number of blocks (exons) in the BED line.
    @property blockCount() const { return this.line.splitter('\t').drop(9).front.to!int; }
    
    /// column 11: A comma-separated list of the block sizes. 
    /// The number of items in this list should correspond to blockCount.
    @property blockSizes() const 
    {
        auto str = this.line.splitter('\t').drop(10).front;
        int[] arr = new int[this.blockCount];
        foreach (i, key; str.splitter(',').enumerate)
        {
            arr[i] = key.to!int;
        }
        return arr;
    }

    /// column 12: A comma-separated list of block starts. 
    /// All of the blockStart positions should be calculated relative to chromStart. 
    /// The number of items in this list should correspond to blockCount.
    @property blockStarts() const 
    {
        auto str = this.line.splitter('\t').drop(11).front;
        int[] arr = new int[this.blockCount];
        foreach (i, key; str.splitter(',').enumerate)
        {
            arr[i] = key.to!int;
        }
        return arr;
    }

    /// get column idx as a string
    /// helps if your bed file isn't a ucsc bed file
    auto opIndex(ulong idx)
    {
        return this.line.splitter('\t').drop(idx).front;
    }
}

auto BedReader(string fn)
{
    return RecordReader!BedRecord(fn);
}

auto BedReader(CoordSystem cs)(string fn, ChromCoordinates!cs region, string fnIdx = "")
{
    return RecordReaderRegion!(BedRecord, cs)(fn, region, fnIdx);
}

auto BedReader(CoordSystem cs)(string fn, string chrom, Coordinates!cs coords, string fnIdx = "")
{
    return RecordReaderRegion!(BedRecord, cs)(fn, chrom, coords, fnIdx);
}

debug(dhtslib_unittest) unittest
{
    import std.stdio;
    import htslib.hts_log;
    import std.algorithm : map;
    import std.array : array;
    import std.path : buildPath, dirName;
    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
    hts_log_info(__FUNCTION__, "Testing BedReader");
    hts_log_info(__FUNCTION__, "Loading test file");

    auto bed = BedReader(buildPath(dirName(dirName(dirName(__FILE__))),"htslib","test","tabix","bed_file.bed"));
    auto rec = bed.front;
    assert(rec.contig == "X");
    assert(rec.coordinates == ZBHO(1000, 1100));
    assert(rec.name == "X1");
    assert(rec.score == 500);
    assert(rec.strand == '+');
    assert(rec.thickStart == 1000);
    assert(rec.thickEnd == 1100);
    assert(rec.itemRGB == RGB(255,0,0));
    bed.popFront;

    rec = bed.front;
    assert(rec.contig == "X");
    assert(rec.coordinates == ZBHO(1200, 1300));
    assert(rec.name == "X2");
    assert(rec.score == 500);
    assert(rec.strand == '+');
    assert(rec.thickStart == 1200);
    assert(rec.thickEnd == 1300);
    assert(rec.itemRGB == RGB(255,0,0));

}

debug(dhtslib_unittest) unittest
{
    import std.stdio;
    import htslib.hts_log;
    import htslib.tbx : tbx_index_build2, tbx_conf_gff;
    import std.algorithm : map;
    import std.array : array;
    import std.path : buildPath, dirName;
    import std.utf : toUTFz;
    import std.array : array;

    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
    hts_log_info(__FUNCTION__, "Testing BedReader (Tabix)");
    hts_log_info(__FUNCTION__, "Loading test file");
    
    auto bed = BedReader(
        buildPath(dirName(dirName(dirName(__FILE__))),"htslib","test","tabix","bed_file.bed.gz"),
        ZBHO("X:1000-1400")
        );

    assert(bed.array.length == 2);

}