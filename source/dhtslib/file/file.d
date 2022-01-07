module dhtslib.file.file;

import core.stdc.stdio : SEEK_SET;

import std.stdio;
import std.string : toStringz, fromStringz;
import std.traits : isSomeString;

import dhtslib.memory;
import dhtslib : SAMHeader, VCFHeader, initKstring;
import htslib;

struct HtslibFileMode
{
    char[] openMode;
    char[] formatMode;

    this(string mode, HtslibFileFormatMode[] formats ...)
    {
        openMode = mode.dup;
        formatMode = cast(char[])formats;
    }

    char[] toMode()
    {
        return openMode ~ formatMode ~ '\0';
    }
}

enum HtslibFileFormatMode
{
    Binary = 'b',
    Cram = 'c',
    Gzip = 'g',
    Uncompressed = 'u',
    Bgzf = 'z',
    ZlibCompress0 = '0',
    ZlibCompress1 = '1',
    ZlibCompress2 = '2',
    ZlibCompress3 = '3',
    ZlibCompress4 = '4',
    ZlibCompress5 = '5',
    ZlibCompress6 = '6',
    ZlibCompress7 = '7',
    ZlibCompress8 = '8',
    ZlibCompress9 = '9',
}



struct HtslibFile
{
    /// dhtslib.memory htsFile* rc wrapper
    /// File pointer
    HtsFile fp; 

    /// D File if used
    File f;

    /// dhtslib.memory kstring_t rc wrapper
    /// Since we need the filename as a null terminated string anyway
    /// just use kstring
    Kstring fn;

    const char[] mode;

    /// SAM/BAM/CRAM index 
    HtsIdx idx;

    /// Tabix index 
    Tbx tbx;

    bool eof;

    /// allow HtslibFile to be used as 
    /// underlying ptr type
    alias getFilePtr this;

    /// get underlying file pointer wrapper
    @property nothrow pure @nogc
    ref inout(HtsFile) getFilePtr() inout return
    {
        return this.fp;
    }

    this(T)(T f)
    if ((is(T == string) || is(T == File)))
    {
        this(f, "r");
    }

    this(T1, T2)(T1 f, T2 mode)
    if ((is(T1 == string) || is(T1 == File)) && 
        (isSomeString!T2 || is(T2 == HtslibFileMode)))
    {
        static if(is(T2 == HtslibFileMode)){
            this.mode = mode.toMode();
        } else {
            this.mode = mode.dup;
        }
        // open file
        static if (is(T1 == string))
        {
            this.fn = Kstring(initKstring());
            ks_initialize(this.fn);
            kputsn(f.ptr, f.length, this.fn);
            this.fp = HtsFile(hts_open(this.fn.s, this.mode.ptr));
        }
        else static if (is(T1 == File))
        {
            this.fn = Kstring(initKstring());
            ks_initialize(this.fn);
            kputsn(f.name.ptr, f.name.length, this.fn);

            auto hf = hdopen(f.fileno, m.ptr);
            this.fp = HtsFile(hts_hopen(hf, this.fn.s, m.ptr));
        }
        else static assert(0);

        this.idx = HtsIdx(null);
    }

    HtslibFile dup()
    {
        auto filename = fromStringz(this.fn.s).idup;
        auto newFile = HtslibFile(filename, this.mode.dup);
        return newFile;
    }

    void setExtraThreads(int extra)
    {
        hts_set_threads(this.fp, extra);
    }

    off_t tell()
    {
        if(this.fp.is_bgzf) return bgzf_tell(this.fp.fp.bgzf);
        else return htell(this.fp.fp.hfile);
    }

    void seek(long loc)
    {
        long err;
        if(this.fp.is_bgzf) err = bgzf_seek(this.fp.fp.bgzf, loc, SEEK_SET);
        else err = hseek(this.fp.fp.hfile, loc, SEEK_SET);
        if(err < 0) hts_log_error(__FUNCTION__, "Error seeking htsFile");
    }

    auto readHeader(T)()
    {
        static if(is(T == BamHdr)){
            return BamHdr(sam_hdr_read(this.fp));
        }else static if(is(T == BcfHdr)){
            return BcfHdr(bcf_hdr_read(this.fp));
        }else static assert(0);
        
    }

    HtsIdx loadHtsIndex()
    {
        this.idx = HtsIdx(sam_index_load(this.fp, this.fn.s));
        return this.idx;
    }

    HtsIdx loadHtsIndex(string idxFile)
    {
        this.idx = HtsIdx(sam_index_load2(this.fp, this.fn.s, toStringz(idxFile)));
        return this.idx;
    }

    Tbx loadTabixIndex()
    {
        this.tbx = Tbx(tbx_index_load(this.fn.s));
        return this.tbx;
    }

    Tbx loadTabixIndex(string idxFile)
    {
        this.tbx = Tbx(tbx_index_load2(this.fn.s, toStringz(idxFile)));
        return this.tbx;
    }

    HtsItr getItr(int tid, hts_pos_t beg, hts_pos_t end)
    {
        return HtsItr(sam_itr_queryi(this.idx, tid, beg, end));
    }

    void write(T1, T2)(T1 rec, T2 header)
    {
        long err;
        static if(is(T1 == Bam1) && is(T2 == BamHdr)){
            err = sam_write1(this.fp, header, rec);
        }
        else static if(is(T1 == Bcf1) && is(T2 == BcfHdr)){
            err = bcf_write(this.fp, header, rec);
        }else static assert(0);
        if(err < 0) hts_log_error(__FUNCTION__, "Error writing BAM/BCF record");
    }

    void write(T)(T rec)
    {
        long err;
        static if(is(T == BamHdr)){
            err = sam_hdr_write(this.fp, rec);
        }
        else static if(is(T == BcfHdr)){
            err = bcf_hdr_write(this.fp, rec);
        }else static if(isSomeString!T || is(T == ubyte[])){
            if(this.fp.is_bgzf) err = bgzf_write(this.fp.fp.bgzf, rec.ptr, rec.length);
            else err = hwrite(this.fp.fp.hfile, rec.ptr, rec.length);
        }else static assert(0);
        if(err < 0) hts_log_error(__FUNCTION__, "Error writing BamHdr/BcfHdr or string record");
    }

    void writeln(T)(T rec)
    if(isSomeString!T)
    {
        rec = rec ~ '\n';
        write(rec);
    }

    auto readRecord(T)(T header)
    {
        long err;
        static if(is(T == BamHdr)){
            auto b = bam_init1;
            err = sam_read1(this.fp, header, b);
            auto rec = Bam1(b);

        }
        else static if(is(T == BcfHdr)){
            auto b = bcf_init;
            err = bcf_read(this.fp, header, b);
            auto rec = Bcf1(b);
        }
        if(err < -1) hts_log_error(__FUNCTION__, "Error reading Bam/Bcf record");
        else if(err == -1) eof = true;
        return rec;
    }

    string readln()
    {
        long err;
        auto ks = Kstring(initKstring);
        ks_initialize(ks);
        err = hts_getline(this.fp, cast(int)'\n', ks);
        if(err < -1) hts_log_error(__FUNCTION__, "Error reading Bam/Bcf record");
        else if(err == -1) eof = true;
        auto s = fromStringz(ks.s).idup;
        return s;
    }

}

unittest
{
    {
        auto f = HtslibFile("/tmp/test.txt", "wu");
        f.writeln("hello");
        f.writeln("test");
    }
    {
        auto f = HtslibFile("/tmp/test.txt");
        assert(f.readln == "hello");
        assert(f.readln == "test");
    }
}

unittest
{
    {
        auto f = HtslibFile("/tmp/test.txt.gz", "wg");
        f.writeln("hello");
        f.writeln("test");
    }
    {
        auto f = HtslibFile("/tmp/test.txt.gz");
        assert(f.readln == "hello");
        assert(f.readln == "test");
    }
}

debug(dhtslib_unittest) unittest
{
    import std.path:buildPath,dirName;
    auto fn = buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","range.bam");
    auto f = HtslibFile(fn);
    
    auto header = f.readHeader!BamHdr();
    auto read = f.readRecord(header);
    assert(fromStringz(bam_get_qname(read)) == "HS18_09653:4:1315:19857:61712");
}