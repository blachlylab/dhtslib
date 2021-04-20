module test.cigar_segfault;

import std.stdio;
import std.exception : enforce;
import dhtslib;
import core.thread;
import core.memory: GC;

void main(string[] args)
{
    string cigar_string;
    Cigar cigar;
    {
        auto bam = SAMReader(args[1], 0);
        foreach(rec;bam.allRecords){
            cigar = rec.cigar;
        }
        writeln(cigar);
        cigar_string = cigar.toString;
        GC.collect;
    }
    writeln(cigar);
    enforce(cigar_string == cigar.toString);
}