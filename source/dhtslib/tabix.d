module dhtslib.tabix;

import std.stdio;
import std.string;
import core.stdc.stdlib : malloc, free;

import dhtslib.htslib.hts;
import dhtslib.htslib.kstring;
import dhtslib.htslib.tbx;
//import dhtslib.htslib.regidx;

/** Encapsulates a position-sorted record-oriented NGS flat file,
 *  indexed with Tabix, including BED, GFF3, VCF.
 *
 *  region(string r) returns an InputRange that iterates through rows of the file intersecting or contained within the requested range 
 */
struct TabixIndexedFile {

    htsFile *fp;    /// pointer to htsFile struct
    tbx_t   *tbx;   /// pointer to tabix handle

    string header;  /// NGS flat file's header (if any; e.g. BED may not have one)

    /// Initialize with a complete file path name to the tabix-indexed file
    /// The tabix index (.tbi) must already exist alongside
    this(string fn,string fntbi="")
    {
        debug(dhtslib_debug) { writeln("TabixIndexedFile ctor"); }
        this.fp = hts_open( toStringz(fn), "r");
        if ( !this.fp ) {
            writefln("Could not read %s\n", fn);
            throw new Exception("Couldn't read file");
        }
        //enum htsExactFormat format = hts_get_format(fp)->format;
        if(fntbi!="") this.tbx = tbx_index_load2( toStringz(fn),toStringz(fntbi) );
        else this.tbx = tbx_index_load( toStringz(fn) );
        if (!this.tbx) { 
            writefln("Could not load .tbi index of %s\n", fn );
            throw new Exception("Couldn't load tabix index file");
        }

        loadHeader();
    }
    ~this()
    {
        debug(dhtslib_debug) { writeln("TabixIndexedFile dtor"); }
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
            this.header ~= fromStringz(str.s) ~ '\n';
        }

        debug(dhtslib_debug) { writeln("end loadHeader"); }
    }

    /// tbx.d: const(char **) tbx_seqnames(tbx_t *tbx, int *n);  // free the array but not the values
    @property string[] sequenceNames()
    {
        // TODO: Check for memory leaks (free the array after copy into sequence_names)
        int nseqs;

        string[] sequence_names;

        const(char **) seqnames = tbx_seqnames(this.tbx, &nseqs);

        for(int i; i<nseqs; i++) {
            sequence_names ~= cast(string) fromStringz(seqnames[i]);
        }

        return sequence_names;
    }

    /** region(r)
     *  returns an InputRange that iterates through rows of the file intersecting or contained within the requested range 
     */
    auto region(string r)
    {
        struct Region {

            /** TODO: determine how thread-(un)safe this is (i.e., using a potentially shared *fp and *tbx) */
            private htsFile *fp;
            private tbx_t   *tbx;

            private hts_itr_t *itr;
            private string next;

            // necessary because the alternative strategy of preloading the first row
            // leads to problems when Struct inst is blitted ->
            // re-iterating always returns first row only (since *itr is expended 
            // but first row was preloaded in this.next)
            private bool active;

            this(htsFile *fp, tbx_t *tbx, string r)
            {
                this.fp = fp;
                this.tbx= tbx;

                this.itr = tbx_itr_querys(tbx, toStringz(r) );
                debug(dhtslib_debug) { writeln("Region ctor // this.itr: ", this.itr); }
                if (this.itr) {
                    // Load the first record
                    //this.popFront(); // correction, do not load the first record
                }
                else {
                    // TODO handle error
                    throw new Exception("could not allocate this.itr");
                }
            }
            ~this()
            {
                debug(dhtslib_debug) { writeln("Region dtor // this.itr: ", this.itr); }
                //tbx_itr_destroy(itr);
                //free(this.kstr.s);
            }

            // had to remove "const" property from empty() due to manipulation of this.active
            @property bool empty() {

                if (!this.active) {
                    // this is the first call to empty() (and use of the range)
                    // Let's make it active and attempt to load the first record, if one exists
                    this.active = true;
                    this.popFront();
                }

                if (!this.next) return true;
                else return false;
            }

            @property string front() const {
                return this.next;
            }

            void popFront() {
                // closure over fp and tbx? (i.e. potentially unsafe?)

                // Get next entry
                kstring_t kstr;
                immutable res = tbx_itr_next(this.fp, this.tbx, this.itr, &kstr);
                if (res < 0) {
                    // we are done
                    this.next = null;
                } else {
                    // Otherwise load into next
                    this.next = fromStringz(kstr.s).idup;
                    free(kstr.s);
                }
            }
        }

        return Region(this.fp, this.tbx, r);
    }

}
