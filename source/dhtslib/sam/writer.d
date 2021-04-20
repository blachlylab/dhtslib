module dhtslib.sam.writer;

import std.stdio : stderr, writeln, File;
import std.string : fromStringz, toStringz;
import std.parallelism : totalCPUs;
import std.format;

import htslib.hts : htsFile, hts_open, hts_close, hts_hopen, hts_set_threads;
import htslib.hts : hts_idx_t, hts_itr_t, hts_itr_multi_t, hts_reglist_t, hts_pair32_t;
import htslib.hfile : hdopen, hclose, hFILE;

import htslib.kstring;
import htslib.sam;
import htslib.hts_log;
import dhtslib.sam.record;
import dhtslib.sam.header;
import dhtslib.sam : parseSam;


/// SAM/BAM/CRAM on-disk format.
/// `DEDUCE` will attempt to auto-detect from filename or other means
enum SAMWriterTypes
{
    BAM,
    UBAM,
    SAM,
    CRAM,
    DEDUCE
}

/// Encapsulates a SAM/BAM/CRAM, but as write-only
struct SAMWriter
{
    /// filename; as usable from D
    string filename;

    /// filename \0-terminated C string; reference needed to avoid GC reaping result of toStringz when ctor goes out of scope
    private const(char)* fn;

    /// htsFile
    private htsFile* fp;

    /// hFILE if required
    private hFILE* f;

    /// header struct
    SAMHeader header;

    private kstring_t line;

    /// disallow copying
    @disable this(this);

    /** Create a representation of SAM/BAM/CRAM file from given filename or File

    Params:
        fn =            string filename (complete path passed to htslib; may support S3:// https:// etc.)
        extra_threads = extra threads for compression/decompression
                        -1 use all available cores (default)
                        0  use no extra threads
                        >1 add indicated number of threads (to a default of 1)
    */
    this(T)(T f,SAMHeader header, SAMWriterTypes t=SAMWriterTypes.DEDUCE,int extra_threads = -1)
    if (is(T == string) || is(T == File))
    {
        import std.parallelism : totalCPUs;
        char[] mode;
        if(t == SAMWriterTypes.BAM) mode=['w','b','\0'];
        else if(t == SAMWriterTypes.UBAM) mode=['w','b','u','\0'];
        else if(t == SAMWriterTypes.SAM) mode=['w','\0'];
        else if(t == SAMWriterTypes.CRAM) mode=['w','c','\0'];
        // open file
        static if (is(T == string))
        {
            if(t == SAMWriterTypes.DEDUCE){
                import std.path:extension;
                auto ext=extension(f);
                if(ext==".bam") mode=['w','b','\0'];
                else if(ext==".sam") mode=['w','\0'];
                else if(ext==".cram") mode=['w','c','\0'];
                else {
                    hts_log_error(__FUNCTION__,"extension "~ext~" not valid");
                    throw new Exception("DEDUCE SAMWriterType used with non-valid extension");
                }
            }
            this.filename = f;
            this.fn = toStringz(f);
            this.fp = hts_open(this.fn, mode.ptr);
        }
        else static if (is(T == File))
        {
            assert(t!=SAMWriterTypes.DEDUCE);
            this.filename = f.name();
            this.fn = toStringz(f.name);
            this.f = hdopen(f.fileno, cast(immutable(char)*) "w");
            this.fp = hts_hopen(this.f, this.fn, mode.ptr);
        }
        else assert(0);

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
        }
        else if (extra_threads > 0)
        {
            if ((extra_threads + 1) > totalCPUs)
                hts_log_warning(__FUNCTION__, "More threads requested than CPU cores detected");
            hts_set_threads(this.fp, extra_threads);
        }
        else if (extra_threads == 0)
        {
            hts_log_debug(__FUNCTION__, "Zero extra threads requested");
        }
        else
        {
            hts_log_warning(__FUNCTION__, "Invalid negative number of extra threads requested");
        }

        // read header
        this.header = header.dup;
        sam_hdr_write(this.fp,this.header.h);
    }

    ~this()
    {
        debug (dhtslib_debug)
        {
            writeln("SAMWriter dtor");
        }

    }

    /// close file
    void close(){
        const auto ret = hts_close(this.fp);
        if (ret < 0)
            stderr.writefln("There was an error closing %s", fromStringz(this.fn));
    }

    /// Write a SAMRecord to disk
    void write(SAMRecord rec){
        const auto ret = sam_write1(this.fp, this.header.h, rec.b);
        assert(ret>=0);
    }
}
///
debug(dhtslib_unittest) unittest
{
    import dhtslib.sam;
    import htslib.hts_log : hts_log_info;
    import std.path : buildPath,dirName;
    import std.string : fromStringz;
    import std.array : array;

    hts_set_log_level(htsLogLevel.HTS_LOG_TRACE);
    hts_log_info(__FUNCTION__, "Testing SAMWriter");
    hts_log_info(__FUNCTION__, "Loading test file");
    auto sam = SAMFile(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","auxf#values.sam"), 0);
    auto sam2 = SAMWriter("/tmp/test.bam",sam.header);
    auto sam3 = SAMWriter("/tmp/test2.bam",sam.header, SAMWriterTypes.DEDUCE, 4);
    auto sam4 = SAMWriter("/tmp/test3.bam",sam.header, SAMWriterTypes.DEDUCE, 0);

    hts_log_info(__FUNCTION__, "Getting read 1");
    auto readrange = sam.allRecords;
    assert(readrange.empty == false);
    auto read = readrange.front();

    sam2.write(read);
    sam3.write(read);
    sam4.write(read);
    sam2.close;
    destroy(sam2);
    sam3.close;
    destroy(sam3);
    sam4.close;
    destroy(sam4);
    auto sam5 = SAMFile("/tmp/test.bam");
    auto sam6 = SAMFile("/tmp/test2.bam");
    auto sam7 = SAMFile("/tmp/test3.bam");

    assert(sam5.allRecords.array.length == 1);
    assert(sam5.allRecords.array.length == sam6.allRecords.array.length);
    assert(sam5.allRecords.array.length == sam7.allRecords.array.length);
    assert(sam7.allRecords.array.length == sam6.allRecords.array.length);

}



