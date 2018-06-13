module dhtslib.tabix;

import std.stdio;
import std.string;
import core.stdc.stdlib : malloc, free;

import dhtslib.htslib.hts;
//import dhtslib.htslib.kstring;
import dhtslib.htslib.tbx;
//import dhtslib.htslib.regidx;

struct TabixIndexedFile {

    htsFile *fp;
    tbx_t   *tbx;

    string header;

    this(string fn)
    {
        debug{ writeln("TabixIndexedFile ctor"); }
        this.fp = hts_open( toStringz(fn), "r");
        if ( !this.fp ) {
            writefln("Could not read %s\n", fn);
            throw new Exception("Couldn't read file");
        }
        //enum htsExactFormat format = hts_get_format(fp)->format;
        
        this.tbx = tbx_index_load( toStringz(fn) );
        if (!this.tbx) { 
            writefln("Could not load .tbi index of %s\n", fn );
            throw new Exception("Couldn't load tabix index file");
        }

        loadHeader();
    }
    ~this()
    {
        debug{ writeln("TabixIndexedFile dtor"); }
        tbx_destroy(this.tbx);

        if ( hts_close(this.fp) ) writefln("hts_close returned non-zero status: %s\n", fromStringz(this.fp.fn));
    }

    private void loadHeader()
    {
        kstring_t str;
        scope(exit) { free(str.s); }

        while ( hts_getline(this.fp, '\n', &str) >= 0 )
        {
            if ( !str.l || str.s[0] != this.tbx.conf.meta_char ) break;
            this.header ~= fromStringz(str.s);
        }

        debug{ writeln("end loadHeader"); }
    }

    /** region(r)
     *  returns an InputRange that iterates through lines from the input file within the given range
     */
    auto region(string r)
    {
        struct Region {

            /** TODO: determine how thread-(un)safe this is */
            private htsFile *fp;
            private tbx_t   *tbx;

            private hts_itr_t *itr;
            private string current;
            private string next;
            private kstring_t kstr;
            
            this(htsFile *fp, tbx_t *tbx, string r)
            {
                this.fp = fp;
                this.tbx= tbx;

                this.itr = tbx_itr_querys(tbx, toStringz(r) );
                writeln("Region ctor");
                writeln("this.itr: ", this.itr);
                if (this.itr) {
                    // Load the first record
                    this.popFront();
                }
                else {
                    // TODO handle error
                    throw new Exception("could not allocate this.itr");
                }
            }
            ~this()
            {
                writeln("Region dtor");
                writeln("this.itr: ", this.itr);
                //tbx_itr_destroy(itr);
                //free(this.kstr.s);
            }

            @property bool empty() const {
                if (!this.next) return true;
                return false;
            }

            @property string front() {
                return this.current;
            }

            void popFront() {
                // closure over fp and tbx? (i.e. potentially unsafe?)

                // Move next into current (note that the ctor runs popFront() once automatically)
                this.current = this.next;

                // Get next entry
                auto res = tbx_itr_next(this.fp, this.tbx, this.itr, &this.kstr);
                if (res < 0) {
                    // we are done
                    this.next = null;
                } else {
                // Otherwise load into next
                this.next = fromStringz(this.kstr.s).idup;
                }
            }
        }

        return Region(this.fp, this.tbx, r);
    }

}
