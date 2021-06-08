module dhtslib.gff.gffreader;

import std.range : inputRangeObject, InputRangeObject;

import dhtslib.bgzf;
import dhtslib.tabix;
import dhtslib.gff3d;
import dhtslib.coordinates;

auto GFF3Reader(string fn)
{
    return RecordReader!GFF3Record(fn);
}

auto GFF3Reader(CoordSystem cs)(string fn, ChromCoordinates!cs region, string fnIdx = "")
{
    return RecordReaderRegion!(GFF3Record, cs)(fn, region, fnIdx);
}

auto GFF3Reader(CoordSystem cs)(string fn, string chrom, Coordinates!cs coords, string fnIdx = "")
{
    return RecordReaderRegion!(GFF3Record, cs)(fn, chrom, coords, fnIdx);
}

auto GTFReader(string fn)
{
    return RecordReader!GTFRecord(fn);
}

auto GTFReader(CoordSystem cs)(string fn, ChromCoordinates!cs region, string fnIdx = "")
{
    return RecordReaderRegion!(GTFRecord, cs)(fn, region, fnIdx);
}

auto GTFReader(CoordSystem cs)(string fn, string chrom, Coordinates!cs coords, string fnIdx = "")
{
    return RecordReaderRegion!(GTFRecord, cs)(fn, chrom, coords, fnIdx);
}

auto GFF2Reader(string fn)
{
    return RecordReader!GTFRecord(fn);
}

auto GFF2Reader(CoordSystem cs)(string fn, ChromCoordinates!cs region, string fnIdx = "")
{
    return RecordReaderRegion!(GTFRecord, cs)(fn, region, fnIdx);
}

auto GFF2Reader(CoordSystem cs)(string fn, string chrom, Coordinates!cs coords, string fnIdx = "")
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

    auto gff = GFF3Reader(buildPath(dirName(dirName(dirName(__FILE__))),"htslib","test","tabix","gff_file.gff"));
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

// debug(dhtslib_unittest) unittest
// {
//     import std.stdio;
//     import htslib.hts_log;
//     import htslib.tbx : tbx_index_build2, tbx_conf_gff;
//     import std.algorithm : map;
//     import std.array : array;
//     import std.path : buildPath, dirName;
//     import std.utf : toUTFz;
//     import std.array : array;

//     hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
//     hts_log_info(__FUNCTION__, "Testing GFF3Reader");
//     hts_log_info(__FUNCTION__, "building test idx file");
//     auto err = tbx_index_build2(
//         toUTFz!(char *)(buildPath(dirName(dirName(dirName(__FILE__))),"htslib","test","tabix","gff_file.gff.gz")),
//         toUTFz!(char *)("test.tbi"),
//         0,
//         &tbx_conf_gff
//     );
//     writeln(err);
//     hts_log_info(__FUNCTION__, "Loading test file");
    
//     auto gff = GFF3Reader(
//         buildPath(dirName(dirName(dirName(__FILE__))),"htslib","test","tabix","gff_file.gff.gz"),
//         OBC("X:2934832-2935190")
//         );

//     assert(gff.array.length == 4);

// }