/* The MIT License

   Copyright (c) 2008, 2012, 2014, 2021 Genome Research Ltd (GRL).
                 2010 by Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/
module htslib.knetfile;

import core.sys.posix.fcntl;
import core.sys.posix.sys.types;

extern (C):

// alias netread = read;
// alias netwrite = write;
// alias netclose = close;

// FIXME: currently I/O is unbuffered

enum KNF_TYPE_LOCAL = 1;
enum KNF_TYPE_FTP = 2;
enum KNF_TYPE_HTTP = 3;

// Kept for API/ABI compatability only.  Do not use directly!
struct knetFile_s
{
    int type;
    int fd;
    long offset;
    char* host;
    char* port;

    // the following are for FTP only
    int ctrl_fd;
    int[4] pasv_ip;
    int pasv_port;
    int max_response;
    int no_reconnect;
    int is_ready;
    char* response;
    char* retr;
    char* size_cmd;
    long seek_offset; // for lazy seek
    long file_size;

    // the following are for HTTP only
    char* path;
    char* http_host;
}

alias knetFile = knetFile_s;

extern (D) auto knet_tell(T)(auto ref T fp)
{
    return fp.offset;
}

extern (D) auto knet_fileno(T)(auto ref T fp)
{
    return fp.fd;
}

knetFile* knet_open (const(char)* fn, const(char)* mode);

/*
   This only works with local files.
 */
knetFile* knet_dopen (int fd, const(char)* mode);

/*
  If ->is_ready==0, this routine updates ->fd; otherwise, it simply
  reads from ->fd.
 */
ssize_t knet_read (knetFile* fp, void* buf, size_t len);

/*
  This routine only sets ->offset and ->is_ready=0. It does not
  communicate with the FTP server.
 */
off_t knet_seek (knetFile* fp, off_t off, int whence);
int knet_close (knetFile* fp);

