module test.samreader;

import std.stdio;
import std.string;

import dhtslib.sam;
import htslib.sam;
import htslib.hts_log;

int main()
{
    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
    auto sf = SAMFile("./htslib/test/range.bam");

    writeln("Basic SAM/BAM/CRAM data\n-----------------------");
    writefln("Input file: %s", sf.filename);
    writefln("N targets: %d", sf.n_targets);
    writefln("Length of contig 0: %d", sf.target_len(0));
    writefln("Length of all targets: %s", sf.target_lens);
    writefln("Target names array: %s", sf.target_names);
    writefln("Lookup by name: \"CHROMOSOME_I\" is the %d'th target (0-indexed).", sf.target_id("CHROMOSOME_I"));

    hts_log_info(__FUNCTION__, "Default class constructor");
    for(int j; j<100_000; j++)
    {
        auto x = new SAMRecord();
    }

    hts_log_info(__FUNCTION__, "Query via htslib C interface");
    {
        int j = 0;
        auto fn = toStringz("./htslib/test/range.bam");

        import htslib.hts;
        auto fp = hts_open(fn, cast(immutable(char)*)"r".ptr);
        auto idx= sam_index_load(fp, fn);
        bam1_t *b = bam_init1();
        hts_itr_t *iter;
        int r;
        if ((iter = sam_itr_queryi(idx, 0, 1_000, 2_000)) == null) {
            hts_log_error(__FUNCTION__, "Failed to parse region");
            return 1;
        }
        writefln("iter == %x", iter);
        
        while ((r = sam_itr_next(fp, iter, b)) >= 0) {
            j++;
        }

        writefln("Processed %d records with raw iter", j);
        assert(j == 14);

        hts_itr_destroy(iter);
        bam_destroy1(b);
        hts_close(fp);
    }

    hts_log_info(__FUNCTION__, "Query via D wrapper (SAMFile.query())");
    {
        int j;
        //auto qr = sf.query("chr1:1000000-2000000");
        //auto qr = sf.query(0, 1_000_000, 2_000_000);
        auto qr = sf.query(0, 1_000, 2_000);
        foreach(r; qr) {
            j++;
        }
        writefln("%d records", j);
        assert(j == 14);
    }

    hts_log_info(__FUNCTION__, "Test query with zero expected results");
    int j = 0;
    auto qr = sf.query(0, 1, 2);
    foreach(r; qr) {
        j++;
    }
    writefln("%d records", j);
    assert(j == 0);

    hts_log_info(__FUNCTION__, "Test AllRecordsRange iterator");
    int i;
    // TODO: remove below once RecordRange is fixed
    sf = SAMFile("./htslib/test/range.bam");
    auto x = sf.all_records;
    foreach(r; x) {
        i++;
        //writeln(r.queryName);
    }
    writefln("%d records", i);
    assert(i == 112);   // confirmed by samtools flagstat

    return 0;
}
