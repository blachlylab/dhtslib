module test.samreader;

import std.stdio;
import std.string;

import dhtslib.sam;
import dhtslib.htslib.sam;
import dhtslib.htslib.hts_log;

int main()
{
    debug(dhtslib_debug)
    {
        writeln("enabling debug logging");
        hts_set_log_level(htsLogLevel.HTS_LOG_DEBUG);
    }

    auto sf = SAMFile("/Users/james/Documents/Development/blachlylab/funmap/wgEncodeUwRepliSeqBg02esG1bAlnRep1.bam");

    writeln("Basic SAM/BAM/CRAM data");
    writefln("N targets: %d", sf.n_targets);
    writefln("Length of contig 0: %d", sf.target_len(0));
    writefln("Length of all targets: %s", sf.target_lens);
    writefln("Target names array: %s", sf.target_names);
    writefln("Lookup by name: \"chr5\" is the %d'th target (0-indexed).", sf.target_id("chr5"));

    int i;
    auto x = sf.all_records;
    foreach(r; x) {
        i++;
        writeln(i);
        writeln(fromStringz(r.sequence));
    }

 

    return 0;
}