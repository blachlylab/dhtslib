/**
Module providing a wrapper, `TabixIndexedFile` over a line-oriented NGS flat file,
such as BED, GFF3, VCF that has been indexed with tabix.

The wrapper provides a list of reference sequence names, as well as iterator over
all rows falling within a sequence range, e.g. "chr1:1000-2000"
*/

module dhtslib.tabix;

import std.stdio;
import std.string;
import std.traits : ReturnType;
import std.range : inputRangeObject, InputRangeObject;
import std.parallelism : totalCPUs;
import core.stdc.stdlib : malloc, free;

import htslib.hts;
import htslib.hts_log;
import htslib.kstring;
import htslib.tbx;
import dhtslib.coordinates;
import dhtslib.memory;
import dhtslib.util;
//import htslib.regidx;

/** Encapsulates a position-sorted record-oriented NGS flat file,
 *  indexed with Tabix, including BED, GFF3, VCF.
 *
 *  region(string r) returns an InputRange that iterates through rows of the file intersecting or contained within the requested range 
 */
struct TabixIndexedFile {

    HtsFile fp;    /// rc wrapper around htsFile ptr
    Tbx   tbx;   /// rc wrapper around tabix ptr

    string header;  /// NGS flat file's header (if any; e.g. BED may not have one)

    /// Initialize with a complete file path name to the tabix-indexed file
    /// The tabix index (.tbi) must already exist alongside
    this(const(char)[] fn, int extra_threads = -1, const(char)[] fntbi = "")
    {
        debug(dhtslib_debug) { writeln("TabixIndexedFile ctor"); }
        this.fp = HtsFile(hts_open( toStringz(fn), "r"));
        if ( !this.fp ) {
            writefln("Could not read %s\n", fn);
            throw new Exception("Couldn't read file");
        }

        if (extra_threads == -1)
        {
            if ( totalCPUs > 1)
            {
                hts_log_info(__FUNCTION__,
                        format("%d CPU cores detected; enabling multithreading", totalCPUs));
                // hts_set_threads adds N _EXTRA_ threads, so totalCPUs - 1 seemed reasonable,
                // but overcomitting by 1 thread (i.e., passing totalCPUs) buys an extra 3% on my 2-core 2013 Mac
                hts_set_threads(this.fp, totalCPUs);
            }
        } else if (extra_threads > 0)
        {
            if ((extra_threads + 1) > totalCPUs)
                hts_log_warning(__FUNCTION__, "More threads requested than CPU cores detected");
            hts_set_threads(this.fp, extra_threads);
        }
        else if (extra_threads == 0)
        {
            hts_log_debug(__FUNCTION__, "Zero extra threads requested");
        }

        //enum htsExactFormat format = hts_get_format(fp)->format;
        if(fntbi!="") this.tbx = Tbx(tbx_index_load2( toStringz(fn), toStringz(fntbi) ));
        else this.tbx = Tbx(tbx_index_load( toStringz(fn) ));
        if (!this.tbx) { 
            writefln("Could not load .tbi index of %s\n", fn );
            throw new Exception("Couldn't load tabix index file");
        }

        loadHeader();
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
    auto region(CoordSystem cs)(string chrom, Interval!cs coords)
    {
        auto newCoords = coords.to!(CoordSystem.zbho);
        struct Region
        {

            /** TODO: determine how thread-(un)safe this is (i.e., using a potentially shared *fp and *tbx) */
            private HtsFile fp;
            private Tbx tbx;

            private HtsItr itr;
            private string next;

            // necessary because the alternative strategy of preloading the first row
            // leads to problems when Struct inst is blitted ->
            // re-iterating always returns first row only (since *itr is expended 
            // but first row was preloaded in this.next)
            private bool active;
            private string chrom;

            this(HtsFile fp, Tbx tbx,  string chrom, Interval!(CoordSystem.zbho) coords)
            {
                this.fp = fp;
                this.tbx = tbx;
                this.chrom = chrom;
                this.itr = HtsItr(tbx_itr_queryi(tbx, tbx_name2id(tbx, toStringz(this.chrom)), coords.start, coords.end));
                this.empty;
                debug(dhtslib_debug) { writeln("Region ctor // this.itr: ", this.itr); }
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
            Region save()
            {
                Region newRange;
                newRange.fp = HtsFile(copyHtsFile(fp));
                newRange.itr = HtsItr(copyHtsItr(itr));
                newRange.tbx = tbx;
                newRange.next = next;
                newRange.active = active;
                newRange.chrom = chrom;
                return newRange;
            }
        }

        return Region(this.fp, this.tbx, chrom, newCoords);
    }

}

// TODO: figure out how to make this unittest with just htslib files

// debug(dhtslib_unittest) unittest
// {
//     import htslib.hts_log;
//     import dhtslib.vcf;
//     import std.path : buildPath, dirName;

//     hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
//     hts_log_info(__FUNCTION__, "Testing TabixIndexedFile");
//     hts_log_info(__FUNCTION__, "Loading test file");
//     auto vcf = VCFReader(buildPath(dirName(dirName(dirName(__FILE__))),"htslib","test","index.vcf"));
//     auto vcfw = VCFWriter()

// }

/**
    Range that allows reading a record based format via tabix.
    Needs a record type that encompasses only one line of text
    and a ChromInterval region to use for tabix filtering.
    Rectype could be GFF3Record, BedRecord ...
    This is a sister struct to dhtslib.bgzf.RecordReader.
*/
struct RecordReaderRegion(RecType, CoordSystem cs)
{
    /// file reader
    TabixIndexedFile file;
    /// file reader range
    ReturnType!(this.initializeRange) range;
    /// chrom of region
    string chrom;
    /// coordinates of region
    Interval!cs coords;
    /// keep the header
    string header;
    
    bool emptyLine = false;
    
    /// string chrom and Interval ctor
    this(string fn, string chrom, Interval!cs coords, int extra_threads = -1, string fnIdx = "")
    {
        this.file = TabixIndexedFile(fn, extra_threads, fnIdx);
        this.chrom = chrom;
        this.coords = coords;
        this.header = this.file.header;
        this.range = this.initializeRange;
        this.range.empty;
    }

    /// copy the TabixIndexedFile.region range
    auto initializeRange()
    {
        return this.file.region(this.chrom, this.coords);
    }

    /// returns RecType
    RecType front()
    {
        return RecType(this.range.front);
    }

    /// move the range
    void popFront()
    {
        this.range.popFront;
        if(this.range.front == "") this.emptyLine = true;
    }

    /// is the range done
    auto empty()
    {
        return this.emptyLine || this.range.empty;
    }

    typeof(this) save()
    {
        typeof(this) newRange = this;
        newRange.range = this.range.save;
        return newRange;
    }
}