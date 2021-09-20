/// @file htslib/hfile.h
/// Buffered low-level input/output streams.
/*
    Copyright (C) 2013-2020 Genome Research Ltd.

    Author: John Marshall <jm18@sanger.ac.uk>

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */
module htslib.hfile;

import core.sys.posix.sys.types;
import core.stdc.string : memcpy, strlen;

import htslib.kstring : kstring_t;

@system:
nothrow:
@nogc:

extern (C):

//#include <sys/types.h>
alias off_t = size_t;
alias ssize_t = size_t;

/// internal
struct hFILE_backend; // @suppress(dscanner.style.phobos_naming_convention)

/// Low-level input/output stream handle
/** The fields of this structure are declared here solely for the benefit
of the hFILE-related inline functions.  They may change in future releases.
User code should not use them directly; you should imagine that hFILE is an
opaque incomplete type.
*/
struct hFILE
{
    import std.bitmanip : bitfields;

    // @cond internal
    char* buffer;
    char* begin;
    char* end;
    char* limit;
    const(hFILE_backend)* backend;
    off_t offset;

    mixin(bitfields!(
        uint, "at_eof", 1,
        uint, "mobile", 1,
        uint, "readonly", 1,
        uint, "", 5));

    int has_errno;
    // @endcond
}

/// Open the named file or URL as a stream
/** @return An hFILE pointer, or `NULL` (with _errno_ set) if an error occurred.

The usual `fopen(3)` _mode_ letters are supported: one of
`r` (read), `w` (write), `a` (append), optionally followed by any of
`+` (update), `e` (close on `exec(2)`), `x` (create exclusively),
`:` (indicates scheme-specific variable arguments follow).
*/
hFILE* hopen(const(char)* filename, const(char)* mode, ...);

/// Associate a stream with an existing open file descriptor
/** @return An hFILE pointer, or `NULL` (with _errno_ set) if an error occurred.

Note that the file must be opened in binary mode, or else
there will be problems on platforms that make a difference
between text and binary mode.

For socket descriptors (on Windows), _mode_ should contain `s`.
*/
hFILE* hdopen(int fd, const(char)* mode);

/// Report whether the file name or URL denotes remote storage
/** @return  0 if local, 1 if remote.

"Remote" means involving e.g. explicit network access, with the implication
that callers may wish to cache such files' contents locally.
*/
int hisremote(const(char)* filename);

/// Append an extension or replace an existing extension
/** @param buffer     The kstring to be used to store the modified filename
    @param filename   The filename to be (copied and) adjusted
    @param replace    If non-zero, one extension (if any) is removed first
    @param extension  The extension to be added (e.g. ".csi")
    @return  The modified filename (i.e., `buffer.s`), or NULL on error.
    @since   1.10

If _filename_ is an URL, alters extensions at the end of the `hier-part`,
leaving any trailing `?query` or `#fragment` unchanged.
*/
char* haddextension(
    kstring_t* buffer,
    const(char)* filename,
    int replace,
    const(char)* extension);

/// Flush (for output streams) and close the stream
/** @return  0 if successful, or `EOF` (with _errno_ set) if an error occurred.
*/
int hclose(hFILE* fp);

/// Close the stream, without flushing or propagating errors
/** For use while cleaning up after an error only.  Preserves _errno_.
*/
void hclose_abruptly(hFILE* fp);

/// Return the stream's error indicator
/** @return  Non-zero (in fact, an _errno_ value) if an error has occurred.

This would be called `herror()` and return true/false to parallel `ferror(3)`,
but a networking-related `herror(3)` function already exists.
*/
pragma(inline, true)
int herrno(hFILE* fp)
{
    return fp.has_errno;
}

/// Clear the stream's error indicator
pragma(inline, true)
void hclearerr(hFILE* fp)
{
    fp.has_errno = 0;
}

/// Reposition the read/write stream offset
/** @return  The resulting offset within the stream (as per `lseek(2)`),
    or negative if an error occurred.
*/
off_t hseek(hFILE* fp, off_t offset, int whence);

/// Report the current stream offset
/** @return  The offset within the stream, starting from zero.
*/
pragma(inline, true)
off_t htell(hFILE* fp)
{
    return fp.offset + (fp.begin - fp.buffer);
}

/// Read one character from the stream
/** @return  The character read, or `EOF` on end-of-file or error.
*/
pragma(inline, true)
int hgetc(hFILE* fp)
{
    
    return (fp.end > fp.begin)? cast(ubyte) *(fp.begin++) : hgetc2(fp);
}

int hgetc2(hFILE* );

/// Read from the stream until the delimiter, up to a maximum length
/** @param buffer  The buffer into which bytes will be written
    @param size    The size of the buffer
    @param delim   The delimiter (interpreted as an `unsigned char`)
    @param fp      The file stream
    @return  The number of bytes read, or negative on error.
    @since   1.4

Bytes will be read into the buffer up to and including a delimiter, until
EOF is reached, or _size-1_ bytes have been written, whichever comes first.
The string will then be terminated with a NUL byte (`\0`).
*/
ssize_t hgetdelim(char* buffer, size_t size, int delim, hFILE* fp);

/// Read a line from the stream, up to a maximum length
/** @param buffer  The buffer into which bytes will be written
    @param size    The size of the buffer
    @param fp      The file stream
    @return  The number of bytes read, or negative on error.
    @since   1.4

Specialization of hgetdelim() for a `\n` delimiter.
*/
pragma(inline, true)
ssize_t hgetln(char* buffer, size_t size, hFILE* fp)
{
    return hgetdelim(buffer, size, '\n', fp);
}

/// Read a line from the stream, up to a maximum length
/** @param buffer  The buffer into which bytes will be written
    @param size    The size of the buffer (must be > 1 to be useful)
    @param fp      The file stream
    @return  _buffer_ on success, or `NULL` if an error occurred.
    @since   1.4

This function can be used as a replacement for `fgets(3)`, or together with
kstring's `kgetline()` to read arbitrarily-long lines into a _kstring_t_.
*/
char* hgets(char* buffer, int size, hFILE* fp);

/// Peek at characters to be read without removing them from buffers
/** @param fp      The file stream
    @param buffer  The buffer to which the peeked bytes will be written
    @param nbytes  The number of bytes to peek at; limited by the size of the
                   internal buffer, which could be as small as 4K.
    @return  The number of bytes peeked, which may be less than _nbytes_
             if EOF is encountered; or negative, if there was an I/O error.

The characters peeked at remain in the stream's internal buffer, and will be
returned by later hread() etc calls.
*/
ssize_t hpeek(hFILE* fp, void* buffer, size_t nbytes);

/// Read a block of characters from the file
/** @return  The number of bytes read, or negative if an error occurred.

The full _nbytes_ requested will be returned, except as limited by EOF
or I/O errors.
*/
pragma(inline, true)
ssize_t hread(hFILE* fp, void* buffer, size_t nbytes)
{
    size_t n = fp.end - fp.begin;
    if (n > nbytes) n = nbytes;
    memcpy(buffer, fp.begin, n);
    fp.begin += n;
    return (n == nbytes || !fp.mobile)? cast(ssize_t) n : hread2(fp, buffer, nbytes, n);
}
/// ditto
ssize_t hread2(hFILE* , void* , size_t, size_t);

/// Write a character to the stream
/** @return  The character written, or `EOF` if an error occurred.
*/
pragma(inline, true)
int hputc(int c, hFILE* fp)
{
    if (fp.begin < fp.limit) *(fp.begin++) = cast(char) c;
    else c = hputc2(c, fp);
    return c;
}
/// ditto
int hputc2(int, hFILE* );

/// Write a string to the stream
/** @return  0 if successful, or `EOF` if an error occurred.
*/
pragma(inline, true)
int hputs(const(char)* text, hFILE* fp)
{

    size_t nbytes = strlen(text), n = fp.limit - fp.begin;
    if (n > nbytes) n = nbytes;
    memcpy(fp.begin, text, n);
    fp.begin += n;
    return (n == nbytes)? 0 : hputs2(text, nbytes, n, fp);
}
/// ditto
int hputs2(const (char)*, size_t, size_t, hFILE* );

/// Write a block of characters to the file
/** @return  Either _nbytes_, or negative if an error occurred.

In the absence of I/O errors, the full _nbytes_ will be written.
*/

// Go straight to hwrite2 if the buffer is empty and the request
// won't fit.
pragma(inline, true)
ssize_t hwrite(hFILE* fp, const(void)* buffer, size_t nbytes)
{

    if (!fp.mobile) {
        size_t n = fp.limit - fp.begin;
        if (n < nbytes) {
            hfile_set_blksize(fp, fp.limit - fp.buffer + nbytes);
            fp.end = fp.limit;
        }
    }

    size_t n = fp.limit - fp.begin;
    if (nbytes >= n && fp.begin == fp.buffer) {
        // Go straight to hwrite2 if the buffer is empty and the request
        // won't fit.
        return hwrite2(fp, buffer, nbytes, 0);
    }

    if (n > nbytes) n = nbytes;
    memcpy(fp.begin, buffer, n);
    fp.begin += n;
    return (n==nbytes)? cast(ssize_t) n : hwrite2(fp, buffer, nbytes, n);
}
/// ditto
ssize_t hwrite2(hFILE* , const(void)* , size_t, size_t);
/// set hfile blocksize
int hfile_set_blksize(hFILE* fp, size_t bufsiz);

/// For writing streams, flush buffered output to the underlying stream
/** @return  0 if successful, or `EOF` if an error occurred.

This includes low-level flushing such as via `fdatasync(2)`.
*/
int hflush(hFILE* fp);

/// For hfile_mem: get the internal buffer and it's size from a hfile
/** @return  buffer if successful, or NULL if an error occurred

The buffer returned should not be freed as this will happen when the
hFILE is closed.
*/
char* hfile_mem_get_buffer(hFILE* file, size_t* length);

/// For hfile_mem: get the internal buffer and it's size from a hfile.
/** @return  buffer if successful, or NULL if an error occurred

This is similar to hfile_mem_get_buffer except that ownership of the
buffer is granted to the caller, who now has responsibility for freeing
it.  From this point onwards, the hFILE should not be used for any
purpose other than closing.
*/
char* hfile_mem_steal_buffer(hFILE* file, size_t* length);

/// Fills out sc_list[] with the list of known URL schemes.
/**
 * @param plugin   [in]     Restricts schemes to only those from 'plugin.
 * @param sc_list  [out]    Filled out with the scheme names
 * @param nschemes [in/out] Size of sc_list (in) and number returned (out)
 *
 * Plugin may be passed in as NULL in which case all schemes are returned.
 * Use plugin "built-in" to list the built in schemes.
 * The size of sc_list is determined by the input value of *nschemes.
 * This is updated to return the output size.  It is up to the caller to
 * determine whether to call again with a larger number if this is too small.
 *
 * The return value represents the total number found matching plugin, which
 * may be larger than *nschemes if too small a value was specified.
 *
 * @return the number of schemes found on success.
 *         -1 on failure
 */
int hfile_list_schemes(
    const(char)* plugin,
    const(char)** sc_list,
    int* nschemes);

/// Fills out plist[] with the list of known hFILE plugins.
/*
 * @param plist    [out]    Filled out with the plugin names
 * @param nplugins [in/out] Size of plist (in) and number returned (out)
 *
 * The size of plist is determined by the input value of *nplugins.
 * This is updated to return the output size.  It is up to the caller to
 * determine whether to call again with a larger number if this is too small.
 *
 * The return value represents the total number found, which may be
 * larger than *nplugins if too small a value was specified.
 *
 * @return the number of plugins found on success.
 *         -1 on failure
 */
int hfile_list_plugins(const(char)** plist, int* nplugins);

/// Tests for the presence of a specific hFILE plugin.
/*
 * @param name     The name of the plugin to query.
 *
 * @return 1 if found, 0 otherwise.
 */
int hfile_has_plugin(const(char)* name);

