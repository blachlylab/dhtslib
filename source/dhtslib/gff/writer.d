module dhtslib.gff.writer;

import dhtslib.gff.record;
import dhtslib.recordwriter;

/// Record writer for GFF3Records
alias GFF2Writer = RecordWriter!GFF2Record;

/// Record writer for GFF3Records
alias GFF3Writer = RecordWriter!GFF3Record;

debug(dhtslib_unittest) unittest
{
    import std.stdio;
    import htslib.hts_log;
    import std.range : take;
    import std.algorithm : map;
    import std.array : array;
    import std.path : buildPath, dirName;
    import dhtslib.gff.reader;
    import dhtslib.coordinates;
    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
    hts_log_info(__FUNCTION__, "Testing GFF3Reader");
    hts_log_info(__FUNCTION__, "Loading test file");

    auto gffr = GFF3Reader(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","gff_file.gff"));
    auto gffw = GFF3Writer("/tmp/test.gff", gffr.header);
    auto rec = gffr.front;
    gffw.write(rec);
    gffr.popFront;
    auto moreRecs = gffr.take(3);
    gffw.writeRecords(moreRecs);
    destroy(gffw);
    gffr = GFF3Reader("/tmp/test.gff");
    rec = gffr.front;

    assert(rec.contig == "X");
    assert(rec.source == "Vega");
    assert(rec.type == "exon");
    assert(rec.coordinates == OBC(2934816, 2935190));
    assert(rec.score == -1);
    assert(rec.strand == '-');
    assert(rec.phase == -1);
    assert(rec["Name"] == "OTTHUME00001604789");
    assert(rec["Parent"] == "OTTHUMT00000055643");
    gffr.popFront;

    rec = gffr.front;
    assert(rec.contig == "X");
    assert(rec.source == "Vega");
    assert(rec.type == "gene");
    assert(rec.coordinates == OBC(2934816, 2964270));
    assert(rec.score == -1);
    assert(rec.strand == '-');
    assert(rec.phase == -1);
    assert(rec["Name"] == "OTTHUMG00000137358");

}