module dhtslib.sam;

import core.stdc.stdlib: malloc, free;
import std.format;
import std.parallelism: totalCPUs;
import std.stdio: writeln, writefln;
import std.string: fromStringz, toStringz;

import dhtslib.htslib.hts: htsFile, hts_itr_t;
import dhtslib.htslib.sam;

/**
Encapsulates a SAM/BAM file.
Implements InputRange interface using htslib calls to ().
*/
struct SAMFile {

    /// filename; reference needed to avoid GC reaping result of toStringz when ctor goes out of scope
    private immutable(char)* fn;

    /// SAM/BAM index 
    private hts_idx_t* idx;

    /// htslib data structure representing the BGZF compressed file/stream fp
    //private BGZF* bgzf;

    private kstring_t line;

    // ref counting to prevent closing file multiple times
    // (free is instead now in popFront instead of dtor)
    private int rc = 1;

    // postblit ref counting
    this(this)
    {
        this.rc++;
    }

    ///
    this(string fn)
    {
        debug(dhtslib_debug) { writeln("SAMFile ctor"); }

        // open file
        this.fn = toStringz(fn);
        //this.bgzf = bgzf_open(this.fn, "r");


        // Do not prime the range with popFront(),
        // because otherwise attempting to iterate again will yield the first row (only)

    }
    ~this()
    {
        debug(dhtslib_debug) { writefln("SAMFile dtor | rc=%d", this.rc); }

        if(!--rc) {
            debug(dhtslib_debug) { 
                writefln("SAMFile closing file (rc=%d)", rc);
            }
            // free(this.line.s) not necessary as should be taken care of in popFront
            // (or front() if using pre-primed range and fetching each row in popFront)
            // on top of this, it should never have been malloc'd in this refcount=0 copy
            //if (bgzf_close(this.bgzf) != 0) writefln("hts_close returned non-zero status: %s\n", fromStringz(this.fn));
        }
    }

    /** Query a region and return matching alignments as an InputRange */
    /// Query by chr, start, end
    void query(string chrom, int start, int end)
    {
        string q = format("%s:%d-%d", chrom, start, end);
        return query(q);
    }
    /// Query by string chr:start-end
    void query(string q)
    {

    }
    /// Query by contig id, start, end
    void query(int tid, int start, int end)
    {
        // hts_itr_t *sam_itr_queryi(const hts_idx_t *idx, int tid, int beg, int end);
        //hts_itr_t* qiter = sam_itr_queryi(this.idx, tid, beg, end);

    }

    /// InputRange interface; returned by query()
    struct RecordRange
    {
        private kstring_t line;     // shut up the compiler
        private htsFile *fp;        // this is shared/will point to a copy possibly used by other iterators
        private hts_itr_t *iter;

        private bam1_t *b;          /// This is the alignment object

        @property bool empty()
        {
            // equivalent to htslib ks_release
            this.line.l = 0;
            this.line.m = 0;
            this.line.s = null;
            
            // int bgzf_getline(BGZF *fp, int delim, kstring_t *str);
            //immutable int res = bgzf_getline(this.bgzf, cast(int)'\n', &this.line);
            //return (res < 0 ? true : false);
            return true;
        }
        /// ditto
        void popFront()
        {

            free(this.line.s);

            // equivalent to htslib ks_release
            this.line.l = 0;
            this.line.m = 0;
            this.line.s = null;
            
        }
        /// ditto
        string front()
        {
            auto ret = fromStringz(this.line.s).idup;
            return ret;
        }
    }
}