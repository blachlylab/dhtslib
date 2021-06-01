/// @file htslib/bgzf.h
/// Low-level routines for direct BGZF operations.
/*
   Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology
                 2011, 2012 Attractive Chaos <attractor@live.co.uk>
   Copyright (C) 2009, 2013, 2014, 2017, 2018-2019 Genome Research Ltd

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

/* The BGZF library was originally written by Bob Handsaker from the Broad
 * Institute. It was later improved by the SAMtools developers. */

import core.sys.posix.sys.types;

extern (C):

enum BGZF_BLOCK_SIZE = 0xff00; // make sure compressBound(BGZF_BLOCK_SIZE) < BGZF_MAX_BLOCK_SIZE
enum BGZF_MAX_BLOCK_SIZE = 0x10000;

enum BGZF_ERR_ZLIB = 1;
enum BGZF_ERR_HEADER = 2;
enum BGZF_ERR_IO = 4;
enum BGZF_ERR_MISUSE = 8;
enum BGZF_ERR_MT = 16; // stream cannot be multi-threaded
enum BGZF_ERR_CRC = 32;

struct hFILE;
struct hts_tpool;
struct kstring_t;
struct bgzf_mtaux_t;
struct bgzidx_t;
struct bgzf_cache_t;
struct z_stream_s;

struct BGZF
{
    import std.bitmanip : bitfields;

    mixin(bitfields!(
        uint, "errcode", 16,
        uint, "reserved", 1,
        uint, "is_write", 1,
        uint, "no_eof_block", 1,
        uint, "is_be", 1,
        int, "compress_level", 9,
        uint, "last_block_eof", 1,
        uint, "is_compressed", 1,
        uint, "is_gzip", 1));

    // Reserved bits should be written as 0; read as "don't care"

    int cache_size;
    int block_length;
    int block_clength;
    int block_offset;
    long block_address;
    long uncompressed_address;
    void* uncompressed_block;
    void* compressed_block;
    bgzf_cache_t* cache;
    hFILE* fp; // actual file handle
    bgzf_mtaux_t* mt; // only used for multi-threading
    bgzidx_t* idx; // BGZF index
    int idx_build_otf; // build index on the fly, set by bgzf_index_build_init()
    z_stream_s* gz_stream; // for gzip-compressed files
    long seeked; // virtual offset of last seek
}

/******************
 * Basic routines *
 ******************/

/**
 * Open an existing file descriptor for reading or writing.
 *
 * @param fd    file descriptor
 *              Note that the file must be opened in binary mode, or else
 *              there will be problems on platforms that make a difference
 *              between text and binary mode.
 * @param mode  mode matching /[rwag][u0-9]+/: 'r' for reading, 'w' for
 *              writing, 'a' for appending, 'g' for gzip rather than BGZF
 *              compression (with 'w' only), and digit specifies the zlib
 *              compression level.
 *              Note that there is a distinction between 'u' and '0': the
 *              first yields plain uncompressed output whereas the latter
 *              outputs uncompressed data wrapped in the zlib format.
 * @return      BGZF file handler; 0 on error
 */
BGZF* bgzf_dopen(int fd, const(char)* mode);

alias bgzf_fdopen = bgzf_dopen; // for backward compatibility

/**
 * Open the specified file for reading or writing.
 */
BGZF* bgzf_open(const(char)* path, const(char)* mode);

/**
 * Open an existing hFILE stream for reading or writing.
 */
BGZF* bgzf_hopen(hFILE* fp, const(char)* mode);

/**
 * Close the BGZF and free all associated resources.
 *
 * @param fp    BGZF file handler
 * @return      0 on success and -1 on error
 */
int bgzf_close(BGZF* fp);

/**
 * Read up to _length_ bytes from the file storing into _data_.
 *
 * @param fp     BGZF file handler
 * @param data   data array to read into
 * @param length size of data to read
 * @return       number of bytes actually read; 0 on end-of-file and -1 on error
 */
ssize_t bgzf_read(BGZF* fp, void* data, size_t length);

/**
 * Write _length_ bytes from _data_ to the file.  If no I/O errors occur,
 * the complete _length_ bytes will be written (or queued for writing).
 *
 * @param fp     BGZF file handler
 * @param data   data array to write
 * @param length size of data to write
 * @return       number of bytes written (i.e., _length_); negative on error
 */
ssize_t bgzf_write(BGZF* fp, const(void)* data, size_t length);

/**
 * Write _length_ bytes from _data_ to the file, the index will be used to
 * decide the amount of uncompressed data to be written to each bgzip block.
 * If no I/O errors occur, the complete _length_ bytes will be written (or
 * queued for writing).
 * @param fp     BGZF file handler
 * @param data   data array to write
 * @param length size of data to write
 * @return       number of bytes written (i.e., _length_); negative on error
 */
ssize_t bgzf_block_write(BGZF* fp, const(void)* data, size_t length);

/**
 * Returns the next byte in the file without consuming it.
 * @param fp     BGZF file handler
 * @return       -1 on EOF,
 *               -2 on error,
 *               otherwise the unsigned byte value.
 */
int bgzf_peek(BGZF* fp);

/**
 * Read up to _length_ bytes directly from the underlying stream without
 * decompressing.  Bypasses BGZF blocking, so must be used with care in
 * specialised circumstances only.
 *
 * @param fp     BGZF file handler
 * @param data   data array to read into
 * @param length number of raw bytes to read
 * @return       number of bytes actually read; 0 on end-of-file and -1 on error
 */
ssize_t bgzf_raw_read(BGZF* fp, void* data, size_t length);

/**
 * Write _length_ bytes directly to the underlying stream without
 * compressing.  Bypasses BGZF blocking, so must be used with care
 * in specialised circumstances only.
 *
 * @param fp     BGZF file handler
 * @param data   data array to write
 * @param length number of raw bytes to write
 * @return       number of bytes actually written; -1 on error
 */
ssize_t bgzf_raw_write(BGZF* fp, const(void)* data, size_t length);

/**
 * Write the data in the buffer to the file.
 *
 * @param fp     BGZF file handle
 * @return       0 on success and -1 on error
 */
int bgzf_flush(BGZF* fp);

/**
 * Return a virtual file pointer to the current location in the file.
 * No interpretation of the value should be made, other than a subsequent
 * call to bgzf_seek can be used to position the file at the same point.
 * Return value is non-negative on success.
 */
extern (D) auto bgzf_tell(T)(auto ref T fp)
{
    return (fp.block_address << 16) | (fp.block_offset & 0xFFFF);
}

/**
 * Set the file to read from the location specified by _pos_.
 *
 * @param fp     BGZF file handler
 * @param pos    virtual file offset returned by bgzf_tell()
 * @param whence must be SEEK_SET
 * @return       0 on success and -1 on error
 *
 * @note It is not permitted to seek on files open for writing,
 * or files compressed with gzip (as opposed to bgzip).
 */
long bgzf_seek(BGZF* fp, long pos, int whence);

/**
 * Check if the BGZF end-of-file (EOF) marker is present
 *
 * @param fp    BGZF file handler opened for reading
 * @return      1 if the EOF marker is present and correct;
 *              2 if it can't be checked, e.g., because fp isn't seekable;
 *              0 if the EOF marker is absent;
 *              -1 (with errno set) on error
 */
int bgzf_check_EOF(BGZF* fp);

/** Return the file's compression format
 *
 * @param fp  BGZF file handle
 * @return    A small integer matching the corresponding
 *            `enum htsCompression` value:
 *   - 0 / `no_compression` if the file is uncompressed
 *   - 1 / `gzip` if the file is plain GZIP-compressed
 *   - 2 / `bgzf` if the file is BGZF-compressed
 * @since 1.4
 */
int bgzf_compression(BGZF* fp);

/**
 * Check if a file is in the BGZF format
 *
 * @param fn    file name
 * @return      1 if _fn_ is BGZF; 0 if not or on I/O error
 */
int bgzf_is_bgzf(const(char)* fn);

/*********************
 * Advanced routines *
 *********************/

/**
 * Set the cache size. Only effective when compiled with -DBGZF_CACHE.
 *
 * @param fp    BGZF file handler
 * @param size  size of cache in bytes; 0 to disable caching (default)
 */
void bgzf_set_cache_size(BGZF* fp, int size);

/**
 * Flush the file if the remaining buffer size is smaller than _size_
 * @return      0 if flushing succeeded or was not needed; negative on error
 */
int bgzf_flush_try(BGZF* fp, ssize_t size);

/**
 * Read one byte from a BGZF file. It is faster than bgzf_read()
 * @param fp     BGZF file handler
 * @return       byte read; -1 on end-of-file or error
 */
int bgzf_getc(BGZF* fp);

/**
 * Read one line from a BGZF file. It is faster than bgzf_getc()
 *
 * @param fp     BGZF file handler
 * @param delim  delimiter
 * @param str    string to write to; must be initialized
 * @return       length of the string; -1 on end-of-file; <= -2 on error
 */
int bgzf_getline(BGZF* fp, int delim, kstring_t* str);

/**
 * Read the next BGZF block.
 */
int bgzf_read_block(BGZF* fp);

/**
 * Enable multi-threading (when compiled with -DBGZF_MT) via a shared
 * thread pool.  This means both encoder and decoder can balance
 * usage across a single pool of worker jobs.
 *
 * @param fp          BGZF file handler; must be opened for writing
 * @param pool        The thread pool (see hts_create_threads)
 */
int bgzf_thread_pool(BGZF* fp, hts_tpool* pool, int qsize);

/**
 * Enable multi-threading (only effective when the library was compiled
 * with -DBGZF_MT)
 *
 * @param fp          BGZF file handler; must be opened for writing
 * @param n_threads   #threads used for writing
 * @param n_sub_blks  #blocks processed by each thread; a value 64-256 is recommended
 */
int bgzf_mt(BGZF* fp, int n_threads, int n_sub_blks);

/**
 * Compress a single BGZF block.
 *
 * @param dst    output buffer (must have size >= BGZF_MAX_BLOCK_SIZE)
 * @param dlen   size of output buffer; updated on return to the number
 *               of bytes actually written to dst
 * @param src    buffer to be compressed
 * @param slen   size of data to compress (must be <= BGZF_BLOCK_SIZE)
 * @param level  compression level
 * @return       0 on success and negative on error
 */
int bgzf_compress(
    void* dst,
    size_t* dlen,
    const(void)* src,
    size_t slen,
    int level);

/*******************
 * bgzidx routines *
 *******************/

/**
 *  Position BGZF at the uncompressed offset
 *
 *  @param fp           BGZF file handler; must be opened for reading
 *  @param uoffset      file offset in the uncompressed data
 *  @param where        must be SEEK_SET
 *
 *  Returns 0 on success and -1 on error.
 *
 *  @note It is not permitted to seek on files open for writing,
 *  or files compressed with gzip (as opposed to bgzip).
 */
int bgzf_useek(BGZF* fp, off_t uoffset, int where);

/**
 *  Position in uncompressed BGZF
 *
 *  @param fp           BGZF file handler; must be opened for reading
 *
 *  Returns the current offset on success and -1 on error.
 */
off_t bgzf_utell(BGZF* fp);

/**
 * Tell BGZF to build index while compressing.
 *
 * @param fp          BGZF file handler; can be opened for reading or writing.
 *
 * Returns 0 on success and -1 on error.
 *
 * @note This function must be called before any data has been read or
 * written, and in particular before calling bgzf_mt() on the same
 * file handle (as threads may start reading data before the index
 * has been set up).
 */
int bgzf_index_build_init(BGZF* fp);

/// Load BGZF index
/**
 * @param fp          BGZF file handler
 * @param bname       base name
 * @param suffix      suffix to add to bname (can be NULL)
 * @return 0 on success and -1 on error.
 */
int bgzf_index_load(BGZF* fp, const(char)* bname, const(char)* suffix);

/// Load BGZF index from an hFILE
/**
 * @param fp   BGZF file handle
 * @param idx  hFILE to read from
 * @param name file name (for error reporting only; can be NULL)
 * @return 0 on success and -1 on error.
 *
 * Populates @p fp with index data read from the hFILE handle @p idx.
 * The file pointer to @idx should point to the start of the index
 * data when this function is called.
 *
 * The file name can optionally be passed in the @p name parameter.  This
 * is only used for printing error messages; if NULL the word "index" is
 * used instead.
 */
int bgzf_index_load_hfile(BGZF* fp, hFILE* idx, const(char)* name);

/// Save BGZF index
/**
 * @param fp          BGZF file handler
 * @param bname       base name
 * @param suffix      suffix to add to bname (can be NULL)
 * @return 0 on success and -1 on error.
 */
int bgzf_index_dump(BGZF* fp, const(char)* bname, const(char)* suffix);

/// Write a BGZF index to an hFILE
/**
 * @param fp     BGZF file handle
 * @param idx    hFILE to write to
 * @param name   file name (for error reporting only, can be NULL)
 * @return 0 on success and -1 on error.
 *
 * Write index data from @p fp to the file @p idx.
 *
 * The file name can optionally be passed in the @p name parameter.  This
 * is only used for printing error messages; if NULL the word "index" is
 * used instead.
 */

int bgzf_index_dump_hfile(BGZF* fp, hFILE* idx, const(char)* name);

