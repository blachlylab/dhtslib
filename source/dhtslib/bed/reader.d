module dhtslib.bed.reader;

import dhtslib.coordinates;
import dhtslib.bed.record;
import dhtslib.bgzf;
import dhtslib.tabix;

/// returns a BedRecord range from a RecordReader
auto BedReader(string fn)
{
    return RecordReader!BedRecord(fn);
}

/// ditto
auto BedReader(CoordSystem cs)(string fn, string chrom, Interval!cs coords, int extra_threads = -1, string fnIdx = "")
{
    return RecordReaderRegion!(BedRecord, cs)(fn, chrom, coords, extra_threads, fnIdx);
}

