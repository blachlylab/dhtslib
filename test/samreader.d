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
        hts_set_log_level(htsLogLevel.HTS_LOG_TRACE);

        // Test log levels
        hts_log_info(__FUNCTION__, "Testing log levels: expect I(nfo) D(ebug) W(arning) E(rror) T(race)");
        hts_log_debug(__FUNCTION__, "Test debug");
        hts_log_warning(__FUNCTION__, "Test warning");
        hts_log_error(__FUNCTION__, "Test error");
        hts_log_trace(__FUNCTION__, "Test trace (after which will reset to Debug level)");

        hts_set_log_level(htsLogLevel.HTS_LOG_DEBUG);
    }


    //auto sf = SAMFile("/Users/james/Documents/Development/blachlylab/funmap/wgEncodeUwRepliSeqBg02esG1bAlnRep1.bam");
    auto sf = SAMFile("/Users/james/Documents/Development/blachlylab/funmap/ENCFF399AWI.bam");

    writeln("Basic SAM/BAM/CRAM data\n-----------------------");
    writefln("Input file: %s", sf.filename);
    writefln("N targets: %d", sf.n_targets);
    writefln("Length of contig 0: %d", sf.target_len(0));
    writefln("Length of all targets: %s", sf.target_lens);
    writefln("Target names array: %s", sf.target_names);
    writefln("Lookup by name: \"chr5\" is the %d'th target (0-indexed).", sf.target_id("chr5"));

    writeln("Now testing default class constructor");
    for(int j; j<100_000; j++)
    {
        auto x = new SAMRecord();
    }

    writeln("Query c raw htslib");
    {
        int j = 0;
        auto fn = toStringz("/Users/james/Documents/Development/blachlylab/funmap/ENCFF399AWI.bam");
        import dhtslib.htslib.hts;
        auto fp = hts_open(fn, cast(immutable(char)*)"r".ptr);
        auto idx= sam_index_load(fp, fn);
        bam1_t *b = bam_init1();
        hts_itr_t *iter;
        int r;
        if ((iter = sam_itr_queryi(idx, 0, 1_000_000, 2_000_000)) == null) {
            hts_log_error(__FUNCTION__, "Failed to parse region");
            return 1;
        }
        writefln("iter == %x", iter);
        
        while ((r = sam_itr_next(fp, iter, b)) >= 0) {
            j++;
        }

        writefln("Processed %d records with raw iter", j);

        hts_itr_destroy(iter);
        bam_destroy1(b);
        hts_close(fp);
    }

    writeln("Testing query with D wrapper");
    int j;
    //auto qr = sf.query("chr1:1000000-2000000");
    auto qr = sf.query(0, 1_000_000, 2_000_000);
    foreach(r; qr) {
        j++;
    }
    writefln("%d records", j);

    writeln("Testing query with expected no results");
    j = 0;
    qr = sf.query(0, 1, 2);
    foreach(r; qr) {
        j++;
    }
    writefln("%d records", j);

    writeln("Now testing AllRecordsRange");
    int i;
    auto x = sf.all_records;
    foreach(r; x) {
        i++;
        //writeln(i);
        //writeln(fromStringz(r.sequence));
    }
    writefln("%d records", i);

    writeln("SAMFile going out of scope?");
 
    return 0;
}