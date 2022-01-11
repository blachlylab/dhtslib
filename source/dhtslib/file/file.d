module dhtslib.file.file;

import core.stdc.stdio : SEEK_SET;

import std.stdio;
import std.string : toStringz, fromStringz;
import std.utf : toUTFz;
import std.traits : isSomeString;
import std.algorithm : map;
import std.array : array;

import dhtslib.memory;
import dhtslib.file.iterator;
import dhtslib : initKstring;
import htslib;

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


/** 
 * Some shortcut mode strings for different file types
 * As a reminder:
 *      b  binary format (BAM, BCF, etc) rather than text (SAM, VCF, etc)
        c  CRAM format
        g  gzip compressed
        u  uncompressed
        z  bgzf compressed
        [0-9]  zlib compression level
 */
enum HtslibFileWriteMode
{
    Bam             = "wb",
    Cram            = "wc",
    Sam             = "w",
    UncompressedBam = "wb0",
    GzippedSam      = "wg",
    BgzippedSam     = "wz",

    Bcf             = "wb",
    UncompressedBcf = "wb0",
    Vcf             = "w",
    GzippedVcf      = "wg",
    BgzippedVcf     = "wz",

    Text            = "w",
    GzippedText     = "wg",
    BgzippedText    = "wz",
    
}

/** 
 * HtslibFile is an abstraction for htslib's htsFile using dhtslib.memory for reference counting.
 * 
 * Ideally this struct and its methods can replace file io across dhtslib through a 
 * standard interface. It can also centralize the ability to duplicate a file so that 
 * it can be iterated several times. It pairs with HtslibIterator.
 */
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

    char[] mode;

    /// SAM/BAM/CRAM index 
    HtsIdx idx;

    /// SAM/BAM/CRAM header
    BamHdr bamHdr;
    /// BCF/VCF header
    BcfHdr bcfHdr;
    /// text header
    Kstring textHdr;

    /// Tabix index 
    Tbx tbx;

    /// is file end reached?
    bool eof;

    /// offset of file after header
    /// if it has been loaded 
    off_t headerOffset = -1;

    /// allow HtslibFile to be used as 
    /// underlying ptr type
    alias getFilePtr this;

    /// get underlying file pointer wrapper
    @property nothrow pure @nogc
    ref inout(HtsFile) getFilePtr() inout return
    {
        return this.fp;
    }

    /// File or string ctor
    this(T)(T f)
    if ((is(T == string) || is(T == File)))
    {
        this(f, "r");
    }

    /// File or string ctor with mode
    this(T1, T2)(T1 f, T2 mode)
    if ((is(T1 == string) || is(T1 == File)) && (isSomeString!T2 || is(T2 == HtslibFileWriteMode)))
    {
        this.mode = mode.dup ~ '\0';
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

    /// Duplicate the file with the same offset
    HtslibFile dup()
    {
        // open a new file
        auto filename = fromStringz(this.fn.s).idup;
        auto newFile = HtslibFile(filename, this.mode);

        // seek to the current offset
        newFile.seek(this.tell);

        // copy internal fields
        newFile.idx = this.idx;
        newFile.bamHdr = this.bamHdr;
        newFile.bcfHdr = this.bcfHdr;
        newFile.textHdr = this.textHdr;
        newFile.tbx = this.tbx;
        newFile.eof = this.eof;
        newFile.headerOffset = this.headerOffset;
        return newFile;
    }

    /// set extra multithreading
    void setExtraThreads(int extra)
    {
        hts_set_threads(this.fp, extra);
    }

    /// get file offset
    off_t tell()
    {
        if(this.fp.is_bgzf) return bgzf_tell(this.fp.fp.bgzf);
        else return htell(this.fp.fp.hfile);
    }

    /// seek to offset in file
    void seek(long loc)
    {
        long err;
        if(this.fp.is_bgzf) err = bgzf_seek(this.fp.fp.bgzf, loc, SEEK_SET);
        else err = hseek(this.fp.fp.hfile, loc, SEEK_SET);
        if(err < 0) hts_log_error(__FUNCTION__, "Error seeking htsFile");
    }

    /// reset file reading
    void reset()
    {
        this.seek(0);
    }

    /// reset file reading
    void resetToFirstRecord()
    {
        if(headerOffset != -1)
            this.seek(headerOffset);
        else {
            hts_log_error(__FUNCTION__, "Cannot resetToFirstRecord, header offset unknown");
            hts_log_error(__FUNCTION__, "Has the header been loaded?");
        }
    }

    /// Read a SAM/BAM header, VCF/BCF header, or text header and internally store it
    void loadHeader(char commentChar = '#')
    {
        switch(this.fp.format.format){
            case htsExactFormat.sam:
            case htsExactFormat.bam:
            case htsExactFormat.cram:
                this.bamHdr = BamHdr(sam_hdr_read(this.fp));
                break;
            case htsExactFormat.vcf:
            case htsExactFormat.bcf:
                this.bcfHdr = BcfHdr(bcf_hdr_read(this.fp));
                break;
            default:
                // parse header from a text file
                this.textHdr = Kstring(initKstring);
                ks_initialize(this.textHdr);

                auto ks = Kstring(initKstring);
                ks_initialize(ks);

                // peek at first char, if commentChar (#)
                // save the line and repeat
                if(this.fp.is_bgzf){
                    
                    while(true){
                        if(cast(char) bgzf_peek(this.fp.fp.bgzf) == commentChar){
                            hts_getline(this.fp, cast(int)'\n', ks);
                            kputs(ks.s, this.textHdr);
                            kputc(cast(int)'\n', this.textHdr);
                        } else break;
                    }
                }else{
                    while(true){
                        char c;
                        hpeek(this.fp.fp.hfile, &c, 1);
                        if(c  == commentChar){
                            hts_getline(this.fp, cast(int)'\n', ks);
                            kputs(ks.s, this.textHdr);
                            kputc(cast(int)'\n', this.textHdr);
                        } else break;
                    }
                }
        }
        // set header offset
        this.headerOffset = this.tell;
    }

    /// set header for a HtlsibFile 
    void setHeader(T)(T hdr)
    if(is(T == BamHdr) || is(T == BcfHdr))
    {
        static if(is(T == BamHdr))
            this.bamHdr = hdr;
        else static if(is(T == BcfHdr))
            this.bcfHdr = hdr;
        else static assert(0);
    }

    /// write internally stored header
    void writeHeader()
    {
        assert(this.bamHdr.initialized || this.bcfHdr.initialized);
        assert(!(this.bamHdr.initialized && this.bcfHdr.initialized));
        assert((this.bamHdr.initialized && this.bamHdr.getRef != null) 
            || (this.bcfHdr.initialized && this.bcfHdr.getRef != null));
        int err;
        if(this.bamHdr.initialized)
            err = sam_hdr_write(this.fp, this.bamHdr);
        else
            err = bcf_hdr_write(this.fp, this.bcfHdr);
        if(err < 0) hts_log_error(__FUNCTION__, "Error writing SAM/BAM/VCF/BCF header");
    }

    /// load BAM/BCF index
    HtsIdx loadHtsIndex()
    {
        if(this.fp.format.format == htsExactFormat.bam || this.fp.format.format == htsExactFormat.cram)
            this.idx = HtsIdx(sam_index_load(this.fp, cast(const(char)*)this.fn.s));
        else if(this.fp.format.format == htsExactFormat.bcf)
            this.idx = HtsIdx(bcf_index_load(cast(const(char)*)this.fn.s));
        return this.idx;
    }

    /// load SAM/BAM index from filename
    HtsIdx loadHtsIndex(string idxFile)
    {
        if(this.fp.format.format == htsExactFormat.bam || this.fp.format.format == htsExactFormat.cram)
            this.idx = HtsIdx(sam_index_load2(this.fp, this.fn.s, toStringz(idxFile)));
        else if(this.fp.format.format == htsExactFormat.bcf)
            this.idx = HtsIdx(bcf_index_load2(this.fn.s, toStringz(idxFile)));
        return this.idx;
    }

    /// load tabix index
    Tbx loadTabixIndex()
    {
        this.tbx = Tbx(tbx_index_load(this.fn.s));
        return this.tbx;
    }

    /// load tabix index from file
    Tbx loadTabixIndex(string idxFile)
    {
        this.tbx = Tbx(tbx_index_load2(this.fn.s, toStringz(idxFile)));
        return this.tbx;
    }

    /// Query htsFile with tid, start, and end
    /// returns an HtslibIterator that has type T
    /// requires index be loaded first
    auto query(T)(int tid, long beg, long end)
    if(is(T == Bam1) || is(T == Bcf1))
    {
        assert(this.idx.getRef != null);
        static if(is(T == Bam1)){
            auto itr = HtsItr(sam_itr_queryi(this.idx, tid, beg, end));
            return HtslibIterator!(Bam1, false)(this, itr);
        }else static if(is(T == Bcf1)){
            auto itr = HtsItr(bcf_itr_queryi(this.idx, tid, beg, end));
            return HtslibIterator!(Bcf1, false)(this, itr);
        }else static assert(0);
    }

    /// Query htsFile with string region
    /// returns an HtslibIterator that has type T
    /// requires index and header be loaded first
    auto query(T)(string region)
    if(is(T == Bam1) || is(T == Bcf1))
    {
        assert(this.idx.getRef != null);
        static if(is(T == Bam1)){
            auto itr = HtsItr(sam_itr_querys(this.idx, this.bamHdr, toStringz(region)));
            return HtslibIterator!(Bam1, false)(this, itr);
        }else static if(is(T == Bcf1)){
            auto itr = HtsItr(bcf_itr_querys(this.idx, this.bcfHdr, toStringz(region)));
            return HtslibIterator!(Bcf1, false)(this, itr);
        }else static assert(0);
    }

    /// iterate over all records in a file
    /// returns a HtslibIterator
    auto byRecord(T)()
    if(is(T == Bam1) || is(T == Bcf1) || is(T == Kstring))
    {
        return HtslibIterator!T(this);
    }

    /// Query htsFile with array of regions
    /// returns an HtslibIterator that has type T
    /// requires index and header be loaded first
    auto query(string[] regions)
    {
        auto cQueries = regions.map!(toUTFz!(char *)).array;

        assert(this.idx.getRef != null);
        auto itr = HtsItr(sam_itr_regarray(this.idx, this.bamHdr, cQueries.ptr, cast(int)regions.length));
        return HtslibIterator!(Bam1, false)(this, itr);
    }

    /// Query tabix'd htsFile with tid, start, and end
    /// returns an HtslibIterator that has type T
    /// requires tabix be loaded first
    auto queryTabix(T)(int tid, long beg, long end)
    if(is(T == Bam1) || is(T == Bcf1) || is(T == Kstring))
    {
        assert(this.tbx.initialized && this.tbx.getRef);
        static if(is(T == Bam1)){
            auto itr = HtsItr(tbx_itr_queryi(this.tbx, tid, beg, end));
            return HtslibIterator!(Bam1, true)(this, itr);
        }else static if(is(T == Bcf1)){
            auto itr = HtsItr(tbx_itr_queryi(this.tbx, tid, beg, end));
            return HtslibIterator!(Bcf1, true)(this, itr);
        }else static if(is(T == Kstring)){
            auto itr = HtsItr(tbx_itr_queryi(this.tbx, tid, beg, end));
            return HtslibIterator!(Kstring, true)(this, itr);
        } else static assert(0);
    }

    /// Query tabix'd htsFile with tid, start, and end
    /// returns an HtslibIterator that has type T
    /// requires tabix be loaded first
    auto queryTabix(T)(string region)
    if(is(T == Bam1) || is(T == Bcf1) || is(T == Kstring))
    {
        assert(this.tbx.initialized && this.tbx.getRef);
        static if(is(T == Bam1)){
            auto itr = HtsItr(tbx_itr_querys(this.tbx, toStringz(region)));
            return HtslibIterator!(Bam1, true)(this, itr);
        }else static if(is(T == Bcf1)){
            auto itr = HtsItr(tbx_itr_querys(this.tbx, toStringz(region)));
            return HtslibIterator!(Bcf1, true)(this, itr);
        }else static if(is(T == Kstring)){
            auto itr = HtsItr(tbx_itr_querys(this.tbx, toStringz(region)));
            return HtslibIterator!(Kstring, true)(this, itr);
        } else static assert(0);
    }

    /// write SAM/BAM/VCF/BCF record, string, or ubyte data
    /// requires the header be loaded if writing SAM/BAM/VCF/BCF record
    void write(T)(T rec)
    if(isSomeString!T || is(T == Bam1) || is(T == Bcf1) || is(T == Kstring) || is(T == ubyte[]))
    {
        long err;
        static if(is(T == Bam1)){
            err = sam_write1(this.fp, this.bamHdr, rec);
        }else static if(is(T == Bcf1)){
            err = bcf_write(this.fp, this.bcfHdr, rec);
        }else static if(isSomeString!T || is(T == ubyte[])){
            if(this.fp.is_bgzf) err = bgzf_write(this.fp.fp.bgzf, rec.ptr, rec.length);
            else err = hwrite(this.fp.fp.hfile, rec.ptr, rec.length);
        }else static if(is(T == Kstring)){
            if(this.fp.is_bgzf) err = bgzf_write(this.fp.fp.bgzf, rec.s, rec.l);
            else err = hwrite(this.fp.fp.hfile, rec.s, rec.l);
        }else static assert(0);
        if(err < 0) hts_log_error(__FUNCTION__, "Error writing SAM/BAM/VCF/BCF record, string, or ubyte data");
    }

    /// write a string with a newline appended
    void writeln(T)(T rec)
    if(isSomeString!T)
    {
        rec = rec ~ '\n';
        write(rec);
    }

    /// read a BAM/SAM/BCF/VCF record
    auto readRecord(T)()
    {
        long err;
        static if(is(T == Bam1)){
            assert(this.bamHdr.initialized && this.bamHdr.getRef != null);
            auto b = bam_init1;
            err = sam_read1(this.fp, this.bamHdr, b);
            auto rec = Bam1(b);
        }
        else static if(is(T == Bcf1)){
            assert(this.bcfHdr.initialized && this.bcfHdr.getRef != null);
            auto b = bcf_init;
            err = bcf_read(this.fp, this.bcfHdr, b);
            auto rec = Bcf1(b);
        }else static if(is(T == Kstring)){
            auto rec = Kstring(initKstring);
            ks_initialize(rec);
            err = hts_getline(this.fp, cast(int)'\n', rec);
        }
        if(err < -1) hts_log_error(__FUNCTION__, "Error reading Bam/Bcf record");
        else if(err == -1) eof = true;
        return rec;
    }

}

debug(dhtslib_unittest) unittest
{
    hts_log_info(__FUNCTION__, "Testing HtslibFile text reading and writing");
    {
        auto f = HtslibFile("/tmp/htslibfile.test.txt", HtslibFileWriteMode.Text);
        f.writeln("hello");
        f.writeln("test");

        f = HtslibFile("/tmp/htslibfile.test.txt.gz", HtslibFileWriteMode.GzippedText);
        f.writeln("hello");
        f.writeln("test");

        f = HtslibFile("/tmp/htslibfile.test.txt.bgz", HtslibFileWriteMode.BgzippedText);
        f.writeln("hello");
        f.writeln("test");

        f = HtslibFile("/tmp/htslibfile.test.txt.9gz", "w9");
        f.writeln("hello");
        f.writeln("test");
    }
    {
        auto f = HtslibFile("/tmp/htslibfile.test.txt");
        assert(fromStringz(f.readRecord!Kstring.s) == "hello");
        assert(fromStringz(f.readRecord!Kstring.s) == "test");

        f = HtslibFile("/tmp/htslibfile.test.txt.gz");
        assert(fromStringz(f.readRecord!Kstring.s) == "hello");
        assert(fromStringz(f.readRecord!Kstring.s) == "test");

        f = HtslibFile("/tmp/htslibfile.test.txt.bgz");
        assert(fromStringz(f.readRecord!Kstring.s) == "hello");
        assert(fromStringz(f.readRecord!Kstring.s) == "test");

        f = HtslibFile("/tmp/htslibfile.test.txt.9gz");
        assert(fromStringz(f.readRecord!Kstring.s) == "hello");
        assert(fromStringz(f.readRecord!Kstring.s) == "test");

    }
    {
        auto f = HtslibFile("/tmp/htslibfile.test.txt.gz");
        assert(fromStringz(f.readRecord!Kstring.s) == "hello");
        assert(fromStringz(f.readRecord!Kstring.s) == "test");
    }
    {
        auto f = HtslibFile("/tmp/htslibfile.test.txt");
        assert(f.byRecord!Kstring.map!(x=> fromStringz(x.s).idup).array == ["hello", "test"] );
    }
}

debug(dhtslib_unittest) unittest
{
    hts_log_info(__FUNCTION__, "Testing HtslibFile SAM/BAM reading and writing");
    import std.path:buildPath,dirName;
    auto fn = buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","range.bam");
    {
        auto f = HtslibFile(fn);
        f.loadHeader;
        auto h = f.bamHdr;
        auto read = f.readRecord!Bam1();

        f = HtslibFile("/tmp/htslibfile.test.sam", HtslibFileWriteMode.Sam);
        f.setHeader(h);
        f.writeHeader;
        f.write(read);
        
        f = HtslibFile("/tmp/htslibfile.test.bam", HtslibFileWriteMode.Bam);
        f.setHeader(h);
        f.writeHeader;
        f.write(read);

        f = HtslibFile("/tmp/htslibfile.test.ubam", HtslibFileWriteMode.UncompressedBam);
        f.setHeader(h);
        f.writeHeader;
        f.write(read);

        f = HtslibFile("/tmp/htslibfile.test.sam.gz", HtslibFileWriteMode.GzippedSam);
        f.setHeader(h);
        f.writeHeader;
        f.write(read);

        f = HtslibFile("/tmp/htslibfile.test.sam.bgz", HtslibFileWriteMode.BgzippedSam);
        f.setHeader(h);
        f.writeHeader;
        f.write(read);

    }
    {

        auto f = HtslibFile("/tmp/htslibfile.test.sam");
        f.loadHeader;
        auto h = f.bamHdr;
        auto read = f.readRecord!Bam1();
        assert(fromStringz(bam_get_qname(read)) == "HS18_09653:4:1315:19857:61712");
        
        f = HtslibFile("/tmp/htslibfile.test.bam");
        f.loadHeader;
        h = f.bamHdr;
        read = f.readRecord!Bam1();
        assert(fromStringz(bam_get_qname(read)) == "HS18_09653:4:1315:19857:61712");

        f = HtslibFile("/tmp/htslibfile.test.ubam");
        f.loadHeader;
        h = f.bamHdr;
        read = f.readRecord!Bam1();
        assert(fromStringz(bam_get_qname(read)) == "HS18_09653:4:1315:19857:61712");

        f = HtslibFile("/tmp/htslibfile.test.sam.gz");
        f.loadHeader;
        h = f.bamHdr;
        read = f.readRecord!Bam1();
        assert(fromStringz(bam_get_qname(read)) == "HS18_09653:4:1315:19857:61712");

        f = HtslibFile("/tmp/htslibfile.test.sam.bgz");
        f.loadHeader;
        h = f.bamHdr;
        read = f.readRecord!Bam1();
        assert(fromStringz(bam_get_qname(read)) == "HS18_09653:4:1315:19857:61712");
    }
}

debug(dhtslib_unittest) unittest
{
    hts_log_info(__FUNCTION__, "Testing HtslibFile BCF reading and writing");
    import std.path:buildPath,dirName;
    auto fn = buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","vcf_file.vcf.gz");
    {
        auto f = HtslibFile(fn);
        f.loadHeader;
        auto h = f.bcfHdr;
        auto read = f.readRecord!Bcf1();

        f = HtslibFile("/tmp/htslibfile.test.vcf", HtslibFileWriteMode.Vcf);
        f.setHeader(h);
        f.writeHeader;
        f.write(read);
        
        f = HtslibFile("/tmp/htslibfile.test.bcf", HtslibFileWriteMode.Bcf);
        f.setHeader(h);
        f.writeHeader;
        f.write(read);

        f = HtslibFile("/tmp/htslibfile.test.ubcf", HtslibFileWriteMode.UncompressedBcf);
        f.setHeader(h);
        f.writeHeader;
        f.write(read);

        f = HtslibFile("/tmp/htslibfile.test.vcf.gz", HtslibFileWriteMode.GzippedVcf);
        f.setHeader(h);
        f.writeHeader;
        f.write(read);

        f = HtslibFile("/tmp/htslibfile.test.vcf.bgz", HtslibFileWriteMode.BgzippedVcf);
        f.setHeader(h);
        f.writeHeader;
        f.write(read);

    }
    {

        auto f = HtslibFile("/tmp/htslibfile.test.vcf");
        f.loadHeader;
        auto h = f.bcfHdr;
        auto read = f.readRecord!Bcf1();
        assert(read.pos == 3000149);
        
        f = HtslibFile("/tmp/htslibfile.test.bcf");
        f.loadHeader;
        h = f.bcfHdr;
        read = f.readRecord!Bcf1();
        assert(read.pos == 3000149);

        f = HtslibFile("/tmp/htslibfile.test.ubcf");
        f.loadHeader;
        h = f.bcfHdr;
        read = f.readRecord!Bcf1();
        assert(read.pos == 3000149);

        f = HtslibFile("/tmp/htslibfile.test.vcf.gz");
        f.loadHeader;
        h = f.bcfHdr;
        read = f.readRecord!Bcf1();
        assert(read.pos == 3000149);

        f = HtslibFile("/tmp/htslibfile.test.vcf.bgz");
        f.loadHeader;
        h = f.bcfHdr;
        read = f.readRecord!Bcf1();
        assert(read.pos == 3000149);
    }
}