module dhtslib.bgzf;

import core.stdc.stdlib: malloc, free;
import std.parallelism: totalCPUs;
import std.stdio: writeln, writefln;
import std.string: fromStringz, toStringz;

import dhtslib.htslib.bgzf;

/**
Encapsulates a bgzipped (block gzipped) file.
Implements InputRange interface using htslib calls.
*/
struct BGZFile {

    /// filename; reference needed to avoid GC reaping result of toStringz when ctor goes out of scope
    private immutable(char)* fn;

    /// htslib data structure representing the BGZF compressed file/stream fp
    private BGZF* bgzf;

    private bool EOF = true;
    private kstring_t line;

    // ref counting to prevent free() errors
    private int rc = 1;

    // postblit ref counting
    this(this)
    {
        this.rc++;
    }

    ///
    this(string fn)
    {
        debug{ writeln("BGZFFile ctor"); }

        // open file
        this.fn = toStringz(fn);
        this.bgzf = bgzf_open(this.fn, "r");

        // enable multi-threading
        // (only effective if library was compiled with -DBGZF_MT)
        // int bgzf_mt(BGZF *fp, int n_threads, int n_sub_blks);
        // n_sub_blks : blocks per thread; 64-256 recommended
        if(totalCPUs > 1) {
            immutable int ret = bgzf_mt(this.bgzf, totalCPUs, 64);
            debug {
                writefln("Total CPUs: %d", totalCPUs);
                writefln("bgzf_mt() -> %d", ret);
            }
        }

        // Prime range
        popFront();

    }
    ~this()
    {
        debug{ writefln("BGZFFile dtor | rc=%d", this.rc); }

        //if (this.line.s != null) {writeln("freeing"); free(this.line.s); }

/+
        debug{ 
            writefln("(rc=%d)", rc);
            writefln("EOF? %s", this.EOF);
            writefln("this.line.s : %s", fromStringz(this.line.s));
            writefln("this.line.s : %x", this.line.s);
            writefln("&this.line.s: %x", &this.line.s);
        }
+/
        if(!--rc) {
            debug{ 
                writefln("closing file (rc=%d)", rc);
            }
            // free(this.line.s) not necessary as should be taken care of in front()
            // on top of this, it should never have been malloc'd in this refcount=0 copy
            if (bgzf_close(this.bgzf) != 0) writefln("hts_close returned non-zero status: %s\n", fromStringz(this.fn));
        }
    }

    /// InputRange interface
    @property bool empty()
    {
        debug writeln("empty()");
        return this.EOF;
    }
    /// ditto
    void popFront()
    {
        debug writeln("popFront()");
        // equivalent to htslib ks_release
        this.line.l = 0;
        this.line.m = 0;
        this.line.s = null;
        
        // int bgzf_getline(BGZF *fp, int delim, kstring_t *str);
        immutable int res = bgzf_getline(this.bgzf, cast(int)'\n', &this.line);
        if (res < 0) {
            // we are done
            this.EOF = true;
        } else {
            this.EOF = false;
        }
    }
    /// ditto
    string front()
    {
        debug writeln("front()");
        auto ret = fromStringz(this.line.s).idup;
        free(this.line.s); // Now calling front() again will result in undefined behavior
        return ret;
    }
}