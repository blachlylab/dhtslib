module dhtslib.sam;

import core.stdc.stdlib: calloc, free;
import std.format;
import std.parallelism: totalCPUs;
import std.stdio: writeln, writefln;
import std.string: fromStringz, toStringz;

import dhtslib.htslib.hts: htsFile, hts_open, hts_close;
import dhtslib.htslib.hts: hts_itr_t;
import dhtslib.htslib.hts: seq_nt16_str;
import dhtslib.htslib.hts: hts_set_threads;

import dhtslib.htslib.hts_log;
import dhtslib.htslib.kstring;
import dhtslib.htslib.sam;

/**
Encapsulates a SAM/BAM/CRAM record,
using the bam1_t type for memory efficiency,
and the htslib helper functions for speed.
**/
class SAMRecord {
    ///
    bam1_t *b;

    ///
    this()
    {
        hts_log_debug(__FUNCTION__, "ctor()"); /// This line triggers memory error when __FUNCTION__, but not when "Other string"
        //writeln(__FUNCTION__);    // this will not trigger the memory error
        this.b = bam_init1();
        //assert(0);                // This will elide(?) the memory error
        //assert(1 == 2);           // This will elide(?) the memory error
    }
    ///
    this(bam1_t *b)
    {
        hts_log_debug(__FUNCTION__, "ctor(bam1_t *b)");
        this.b = b;
    }
    ~this()
    {
        writeln("Record dtor");
        hts_log_debug(__FUNCTION__, "dtor");
        //bam_destroy1(this.b); // we don't own it!
    }
}

int main()
{
    hts_set_log_level(htsLogLevel.HTS_LOG_DEBUG);

    htsFile *fp = hts_open("/Users/james/Documents/Development/blachlylab/funmap/ENCFF399AWI.bam", cast(immutable(char)*)"r");

    bam_hdr_t *header = sam_hdr_read(fp);

    bam1_t *b = bam_init1();

    for(int i;i<100_000_000;i++)
    {
        sam_read1(fp, header, b);
        auto x = new SAMRecord(bam_dup1(b));
        writeln(i);
    }

    
    bam_hdr_destroy(header);
    hts_close(fp);

    return 0;
}
