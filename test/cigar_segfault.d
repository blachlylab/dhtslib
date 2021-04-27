module test.cigar_segfault;

import std.stdio;
import std.exception : enforce;
import dhtslib;
import htslib.hts_log;
import core.thread;
import core.memory: GC;

void main(string[] args)
{
    //hts_set_log_level(htsLogLevel.HTS_LOG_DEBUG);
    string cigar_string;
    Cigar cigar;
    {
        auto bam = SAMReader("htslib/test/range.bam", 0);
        foreach(rec;bam.allRecords){
            cigar = rec.cigar;
            cigar_string = cigar.toString;
            assert(rec.cigar.references == 1);
        }
        GC.collect;
    }
    enforce(cigar_string == cigar.toString);
}