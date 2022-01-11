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
import std.algorithm : map;
import core.stdc.stdlib : malloc, free;

import htslib.hts;
import htslib.kstring;
import htslib.tbx;
import dhtslib.coordinates;
import dhtslib.memory;
import dhtslib.util;
import dhtslib.file;
//import htslib.regidx;

/** Encapsulates a position-sorted record-oriented NGS flat file,
 *  indexed with Tabix, including BED, GFF3, VCF.
 *
 *  region(string r) returns an InputRange that iterates through rows of the file intersecting or contained within the requested range 
 */
struct TabixIndexedFile {

    HtslibFile f;    /// HtslibFile 
    Tbx   tbx;   /// rc wrapper around tabix ptr
    string header;

    /// Initialize with a complete file path name to the tabix-indexed file
    /// The tabix index (.tbi) must already exist alongside
    this(const(char)[] fn, const(char)[] fntbi = "")
    {
        debug(dhtslib_debug) { writeln("TabixIndexedFile ctor"); }
        this.f = HtslibFile(fn);
        if ( !this.f.fp ) {
            stderr.writefln("Could not read %s\n", fn);
            throw new Exception("Couldn't read file");
        }
        //enum htsExactFormat format = hts_get_format(fp)->format;
        if(fntbi!="") this.f.loadTabixIndex(fntbi);
        else this.f.loadTabixIndex;
        this.tbx = this.f.tbx;
        if (!this.tbx) { 
            stderr.writefln("Could not load .tbi index of %s\n", fn );
            throw new Exception("Couldn't load tabix index file");
        }

        this.f.loadHeader!Kstring('#');
        this.header = fromStringz(this.f.textHdr.s).idup;
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
        auto tid = tbx_name2id(tbx, toStringz(chrom));
        auto newF = this.f.dup;
        newF.resetToFirstRecord;
        return newF.query!Kstring(tid,newCoords.start, newCoords.end)
                    .map!(x => fromStringz(x.s).idup);
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
    this(string fn, string chrom, Interval!cs coords, string fnIdx = "")
    {
        this.file = TabixIndexedFile(fn, fnIdx);
        this.chrom = chrom;
        this.coords = coords;
        this.header = this.file.header;
        this.range = this.initializeRange;
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
        if(!this.range.empty && this.range.front == "") this.emptyLine = true;
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