module dhtslib.bed.bedreader;

import dhtslib.coordinates;
import dhtslib.bed.bedrecord;
import dhtslib.bgzf;
import dhtslib.tabix;

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

