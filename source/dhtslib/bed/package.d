module dhtslib.bed;
/**

BED file reading and writing

This module provides a readable, writeable abstraction of BED records and files.

Authors: Thomas Gregory <charles.gregory@osumc.edu>

Standards: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
*/

public import dhtslib.bed.record;
public import dhtslib.bed.reader;
public import dhtslib.bed.writer;

debug(dhtslib_unittest) unittest
{
    import dhtslib.coordinates;
    import std.stdio;
    import htslib.hts_log;
    import std.algorithm : map;
    import std.array : array;
    import std.path : buildPath, dirName;
    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
    hts_log_info(__FUNCTION__, "Testing BedReader");
    hts_log_info(__FUNCTION__, "Loading test file");

    auto bed = BedReader(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","bed_file.bed"));
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

    rec.contig = "X1";
    rec.coordinates = ZBHO(1201, 1301);
    rec.name = "X21";
    rec.score = 501;
    rec.strand = '-';
    rec.thickStart = 1201;
    rec.thickEnd = 1301;
    rec.itemRGB = RGB(255,0,1);

    assert(rec.contig == "X1");
    assert(rec.coordinates == ZBHO(1201, 1301));
    assert(rec.name == "X21");
    assert(rec.score == 501);
    assert(rec.strand == '-');
    assert(rec.thickStart == 1201);
    assert(rec.thickEnd == 1301);
    assert(rec.itemRGB == RGB(255,0,1));

    assert(rec.toString == "X1\t1201\t1301\tX21\t501\t-\t1201\t1301\t255,0,1");

}

debug(dhtslib_unittest) unittest
{
    import dhtslib.coordinates;
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
        buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","bed_file.bed.gz"),
        ChromZBHO("X:1000-1400")
        );

    assert(bed.array.length == 2);

}