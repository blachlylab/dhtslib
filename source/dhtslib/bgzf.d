module dhtslib.bgzf;

import core.stdc.stdlib: malloc, free;
import std.parallelism: totalCPUs;
import std.stdio: writeln, writefln;
import std.string: fromStringz, toStringz;

import htslib.bgzf;
import htslib.kstring;

/**
Encapsulates a bgzipped (block gzipped) file.
Implements InputRange interface using htslib calls to bgzf_getline().
*/
struct BGZFile {

    /// filename; reference needed to avoid GC reaping result of toStringz when ctor goes out of scope
    private immutable(char)* fn;

    /// htslib data structure representing the BGZF compressed file/stream fp
    private BGZF* bgzf;

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
        debug(dhtslib_debug) { writeln("BGZFile ctor"); }

        // open file
        this.fn = toStringz(fn);
        this.bgzf = bgzf_open(this.fn, "r");

        // enable multi-threading
        // (only effective if library was compiled with -DBGZF_MT)
        // int bgzf_mt(BGZF *fp, int n_threads, int n_sub_blks);
        // n_sub_blks : blocks per thread; 64-256 recommended
        if(totalCPUs > 1) {
            immutable int ret = bgzf_mt(this.bgzf, totalCPUs, 64);
            debug(dhtslib_debug) {
                writefln("Total CPUs: %d", totalCPUs);
                writefln("bgzf_mt() -> %d", ret);
            }
        }

        // Do not prime the range with popFront(),
        // because otherwise attempting to iterate again will yield the first row (only)

    }
    ~this()
    {
        debug(dhtslib_debug) { writefln("BGZFile dtor | rc=%d", this.rc); }

        if(!--rc) {
            debug(dhtslib_debug) { 
                writefln("BGZFile closing file (rc=%d)", rc);
            }
            // free(this.line.s) not necessary as should be taken care of in popFront
            // (or front() if using pre-primed range and fetching each row in popFront)
            // on top of this, it should never have been malloc'd in this refcount=0 copy
            if (bgzf_close(this.bgzf) != 0) writefln("hts_close returned non-zero status: %s\n", fromStringz(this.fn));
        }
    }

    /// InputRange interface
    @property bool empty()
    {
        // equivalent to htslib ks_release
        this.line.l = 0;
        this.line.m = 0;
        this.line.s = null;
        
        // int bgzf_getline(BGZF *fp, int delim, kstring_t *str);
        immutable int res = bgzf_getline(this.bgzf, cast(int)'\n', &this.line);
        return (res < 0 ? true : false);
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
