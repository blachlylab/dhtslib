module test.sam_iter;

import std.stdio;
import std.string;

import dhtslib.htslib.hts;
import dhtslib.htslib.hts_log;
import dhtslib.htslib.sam;

int main()
{
    writeln("Query c raw htslib");
    int j = 0;
    auto fn = toStringz("/Users/james/Documents/Development/blachlylab/funmap/ENCFF399AWI.bam");

    htsFile* fp = hts_open(fn, cast(immutable(char)*)"r".ptr);
    writefln("fp    : %x", fp);
    writefln("fp->fn: %s", fromStringz(fp.fn));
    writefln("fp->fp: %x", fp.fp.bgzf);
    writefln("fp->fmt %s", fp.format);

    auto idx= sam_index_load(fp, fn);
    bam1_t *b = bam_init1();
    hts_itr_t *iter;
    int r;
    if ((iter = sam_itr_queryi(idx, 0, 1_000_000, 2_000_000)) == null) {
        hts_log_error(__FUNCTION__, "Failed to parse region");
        return 1;
    }
    writefln("iter == %x", iter);
    writefln("fp->fp: %x", fp.fp.bgzf);
    while ((r = sam_itr_next(fp, iter, b)) >= 0) {
        j++;
    }

    writefln("Processed %d records with raw iter", j);

    hts_itr_destroy(iter);
    bam_destroy1(b);
    hts_close(fp);

    return 0;
}