module dhtslib.gff.reader;

import std.range : inputRangeObject, InputRangeObject;

import dhtslib.bgzf;
import dhtslib.tabix;
import dhtslib.gff;
import dhtslib.coordinates;

/// Returns a RecordReader range for a GFF3 file
auto GFF3Reader(string fn)
{
    return RecordReader!GFF3Record(fn);
}

/// ditto
auto GFF3Reader(CoordSystem cs)(string fn, string chrom, Interval!cs coords, int extra_threads = -1, string fnIdx = "")
{
    return RecordReaderRegion!(GFF3Record, cs)(fn, chrom, coords, extra_threads, fnIdx);
}

/// Returns a RecordReader range for a GTF file
auto GTFReader(string fn)
{
    return RecordReader!GTFRecord(fn);
}

/// ditto
auto GTFReader(CoordSystem cs)(string fn, string chrom, Interval!cs coords, int extra_threads = -1, string fnIdx = "")
{
    return RecordReaderRegion!(GTFRecord, cs)(fn, chrom, coords, extra_threads, fnIdx);
}

/// Returns a RecordReader range for a GFF2 file
auto GFF2Reader(string fn)
{
    return RecordReader!GTFRecord(fn);
}

/// ditto
auto GFF2Reader(CoordSystem cs)(string fn, string chrom, Interval!cs coords, string fnIdx = "")
{
    return RecordReaderRegion!(GTFRecord, cs)(fn, chrom, coords, fnIdx);
}

debug(dhtslib_unittest) unittest
{
    import std.stdio;
    import htslib.hts_log;
    import std.algorithm : map;
    import std.array : array;
    import std.path : buildPath, dirName;
    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
    hts_log_info(__FUNCTION__, "Testing GFF3Reader");
    hts_log_info(__FUNCTION__, "Loading test file");

    auto gff = GFF3Reader(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","gff_file.gff"));
    auto rec = gff.front;
    assert(rec.contig == "X");
    assert(rec.source == "Vega");
    assert(rec.type == "exon");
    assert(rec.coordinates == OBC(2934816, 2935190));
    assert(rec.score == -1);
    assert(rec.strand == '-');
    assert(rec.phase == -1);
    assert(rec["Name"] == "OTTHUME00001604789");
    assert(rec["Parent"] == "OTTHUMT00000055643");
    gff.popFront;

    rec = gff.front;
    assert(rec.contig == "X");
    assert(rec.source == "Vega");
    assert(rec.type == "gene");
    assert(rec.coordinates == OBC(2934816, 2964270));
    assert(rec.score == -1);
    assert(rec.strand == '-');
    assert(rec.phase == -1);
    assert(rec["Name"] == "OTTHUMG00000137358");

}

debug(dhtslib_unittest) unittest
{
    import dhtslib.util;
    import std.stdio;
    import htslib.hts_log;
    import htslib.tbx : tbx_index_build2, tbx_conf_gff;
    import std.algorithm : map;
    import std.array : array;
    import std.path : buildPath, dirName;
    import std.utf : toUTFz;
    import std.array : array;

    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
    hts_log_info(__FUNCTION__, "Testing GFF3Reader");
    hts_log_info(__FUNCTION__, "building test idx file");
    // auto err = tbx_index_build2(
    //     toUTFz!(char *)(buildPath(dirName(dirName(dirName(__FILE__))),"htslib","test","tabix","gff_file.gff.gz")),
    //     toUTFz!(char *)("test.tbi"),
    //     0,
    //     &tbx_conf_gff
    // );
    // writeln(err);
    hts_log_info(__FUNCTION__, "Loading test file");
    
    auto reg = getIntervalFromString("X:2934832-2935190");
    auto gff = GFF3Reader(
        buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","gff_file.gff.gz"),
        reg.contig, reg.interval
        );

    assert(gff.array.length == 4);

}

debug(dhtslib_unittest) unittest
{
    import std.stdio;
    import htslib.hts_log;
    import std.algorithm : map;
    import std.array : array;
    import std.path : buildPath, dirName;
    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
    hts_log_info(__FUNCTION__, "Testing GFF3Reader save");
    hts_log_info(__FUNCTION__, "Loading test file");

    auto gff = GFF3Reader(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","gff_file.gff"));
    assert(gff.array.length == 62);

}

debug(dhtslib_unittest) unittest
{
    import dhtslib.util;
    import std.stdio;
    import htslib.hts_log;
    import htslib.tbx : tbx_index_build2, tbx_conf_gff;
    import std.algorithm : map;
    import std.array : array;
    import std.path : buildPath, dirName;
    import std.utf : toUTFz;
    import std.array : array;

    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
    hts_log_info(__FUNCTION__, "Testing GFF3Reader save (tabix)");
    hts_log_info(__FUNCTION__, "Loading test file");
    
    auto reg = getIntervalFromString("X:2934832-2935190");
    auto gff = GFF3Reader(
        buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","gff_file.gff.gz"),
        reg.contig, reg.interval
        );
    auto range1 = gff.save;
    gff.popFront;

    auto range2 = gff.save;
    gff.popFront;

    auto range3 = gff.save;
    gff.popFront;

    auto range4 = gff.save;
    gff.popFront;

    auto range5 = gff.save;
    assert(range1.array.length == 4);
    assert(range2.array.length == 3);
    assert(range3.array.length == 2);
    assert(range4.array.length == 1);
    assert(range5.array.length == 0);

}