module dhtslib.bed.writer;

import dhtslib.bed.record;
import dhtslib.recordwriter;

/// Record writer for BedRecords
alias BedWriter = RecordWriter!BedRecord;


debug(dhtslib_unittest) unittest
{
    import std.stdio;
    import htslib.hts_log;
    import std.range : take;
    import std.algorithm : map;
    import std.array : array;
    import std.path : buildPath, dirName;
    import dhtslib.bed.reader;
    import dhtslib.coordinates;
    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
    hts_log_info(__FUNCTION__, "Testing GFF3Reader");
    hts_log_info(__FUNCTION__, "Loading test file");

    auto bedr = BedReader(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","bed_file.bed"));
    auto bedw = BedWriter("/tmp/test.bed", bedr.header);
    auto rec = bedr.front;
    bedw.write(rec);
    bedr.popFront;
    auto moreRecs = bedr.take(3);
    bedw.writeRecords(moreRecs);
    destroy(bedw);

    bedr = BedReader("/tmp/test.bed");
    rec = bedr.front;
    assert(rec.contig == "X");
    assert(rec.coordinates == ZBHO(1000, 1100));
    assert(rec.name == "X1");
    assert(rec.score == 500);
    assert(rec.strand == '+');
    assert(rec.thickStart == 1000);
    assert(rec.thickEnd == 1100);
    assert(rec.itemRGB == RGB(255,0,0));
    bedr.popFront;

    rec = bedr.front;
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