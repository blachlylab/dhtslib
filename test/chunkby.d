module test.chunkby;

import std.stdio;
import dhtslib.sam;
import htslib.hts_log;
import std.path : buildPath, dirName;
import std.string : fromStringz;
import std.array : array; 
import std.algorithm : chunkBy, map;

void main()
{
    hts_set_log_level(htsLogLevel.HTS_LOG_WARNING);
    hts_log_info(__FUNCTION__, "Testing chunkBy");
    hts_log_info(__FUNCTION__, "Loading test file");
    

    auto sam = SAMFile(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","range.bam"), 0);
    /// group by contig, then convert groups to arrays
    /// print info from groups
    writeln("1");
    auto range = sam.allRecords.chunkBy!((a, b) => a.tid == b.tid).map!(x=> x.array);
    foreach(recs; range){
        foreach(rec; recs){
            writeln(rec.queryName);
            writeln(rec.references);
            writeln(rec.pos);
            writeln(rec.cigar.toString);
            writeln(rec.sequence);
            writeln();
            break;
        }
        break;
    }

    /// group by contig, then convert groups to arrays
    /// print info from groups
    writeln("2");
    auto range2 = sam.allRecords.chunkBy!((a, b) => a.tid == b.tid);//.map!(x=> x.array);
    foreach(group; range2){
        auto recs = group.array;
        foreach(rec; recs){
            writeln(rec.queryName);
            writeln(rec.references);
            writeln(rec.pos);
            writeln(rec.cigar.toString);
            writeln(rec.sequence);
            writeln();
            break;
        }
        break;
    }
    /// group by contig, then convert groups to arrays
    /// print info from groups
    writeln("3");
    auto range3 = sam.allRecords.chunkBy!((a, b) => a.tid == b.tid);//.map!(x=> x.array);
    foreach(group; range3){
        SAMRecord[] recs;
        foreach(rec; group){
            recs ~= rec;
        }
        foreach(rec; recs){
            writeln(rec.queryName);
            writeln(rec.references);
            writeln(rec.pos);
            writeln(rec.cigar.toString);
            writeln(rec.sequence);
            writeln();
            break;
        }
        break;
    }

    /// group by contig, DON'T convert groups to arrays
    /// print info from groups
    /// Only this one is correct
    writeln("4");
    auto range4 = sam.allRecords.chunkBy!((a, b) => a.tid == b.tid);//.map!(x=> x.array);
    foreach(recs; range4){
        foreach(rec; recs){
            writeln(rec.queryName);
            writeln(rec.references);
            writeln(rec.pos);
            writeln(rec.cigar.toString);
            writeln(rec.sequence);
            writeln();
            break;
        }
        break;
    }
}