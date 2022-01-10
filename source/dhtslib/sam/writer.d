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
import dhtslib.memory;
import dhtslib.util;
import dhtslib.file;


/// SAM/BAM/CRAM on-disk format.
/// `DEDUCE` will attempt to auto-detect from filename or other means
enum SAMWriterTypes
{
    BAM = HtslibFileWriteMode.Bam,
    UBAM = HtslibFileWriteMode.UncompressedBam,
    SAM = HtslibFileWriteMode.Sam,
    CRAM = HtslibFileWriteMode.Cram,
    DEDUCE = []
}

/// Encapsulates a SAM/BAM/CRAM, but as write-only
struct SAMWriter
{
    /// htsFile
    private HtslibFile f;

    /// header struct
    SAMHeader header;

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
        char[] mode = cast(char[])t;
        // open file
        static if (is(T == string))
        {
            if(t == SAMWriterTypes.DEDUCE){
                import std.path:extension;
                auto ext=extension(f);
                if(ext==".bam") mode=cast(char[])HtslibFileWriteMode.Bam;
                else if(ext==".sam") mode=cast(char[])HtslibFileWriteMode.Sam;
                else if(ext==".cram") mode=cast(char[])HtslibFileWriteMode.Cram;
                else {
                    hts_log_error(__FUNCTION__,"extension "~ext~" not valid");
                    throw new Exception("DEDUCE SAMWriterType used with non-valid extension");
                }
            }
            this.f = HtslibFile(f, mode);
        }
        else static if (is(T == File))
        {
            assert(t!=SAMWriterTypes.DEDUCE);
            this.f = HtslibFile(f, mode);
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
                this.f.setExtraThreads(totalCPUs);
            }
        }
        else if (extra_threads > 0)
        {
            if ((extra_threads + 1) > totalCPUs)
                hts_log_warning(__FUNCTION__, "More threads requested than CPU cores detected");
            this.f.setExtraThreads(extra_threads);
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
        this.f.setHeader(this.header.h);
        this.f.writeHeader;
    }


    /// Write a SAMRecord to disk
    void write(SAMRecord rec){
        this.f.write(rec.b);
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
    destroy(sam2);
    destroy(sam3);
    destroy(sam4);
    auto sam5 = SAMFile("/tmp/test.bam");
    auto sam6 = SAMFile("/tmp/test2.bam");
    auto sam7 = SAMFile("/tmp/test3.bam");

    assert(sam5.allRecords.array.length == 1);
    assert(sam5.allRecords.array.length == sam6.allRecords.array.length);
    assert(sam5.allRecords.array.length == sam7.allRecords.array.length);
    assert(sam7.allRecords.array.length == sam6.allRecords.array.length);

}



