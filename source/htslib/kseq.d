/* The MIT License

   Copyright (c) 2008, 2009, 2011 Attractive Chaos <attractor@live.co.uk>
   Copyright (C) 2013, 2018, 2020 Genome Research Ltd.

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

/* Last Modified: 05MAR2012 */
module htslib.kseq;

extern (C):

/* klib_unused */

enum KS_SEP_SPACE = 0; // isspace(): \t, \n, \v, \f, \r
enum KS_SEP_TAB = 1; // isspace() && !' '
enum KS_SEP_LINE = 2; // line separator: "\n" (Unix) or "\r\n" (Windows)
enum KS_SEP_MAX = 2;

extern (D) auto ks_eof(T)(auto ref T ks)
{
    return ks.is_eof && ks.begin >= ks.end;
}

/* never come to here! */

/******************
 * FASTA/Q parser *
 ******************/

/* Return value:
   >=0  length of the sequence (normal)
   -1   end-of-file
   -2   truncated quality string
 */

/* then jump to the next header line */

/* end of file */

/* else: the first header char has been read in the previous call */
/* reset all members */
/* normal exit: EOF */
/* read FASTA/Q comment */
/* we can do this in the loop below, but that is slower */

/* skip empty lines */
/* this is safe: we always have enough space for 1 char */
/* read the rest of the line */

/* the first header char has been read */
/* seq->seq.s[seq->seq.l] below may be out of boundary */

/* rounded to the next closest 2^k */

/* null terminated string */
/* FASTA */
/* allocate memory for qual in case insufficient */

/* skip the rest of '+' line */
/* error: no quality string */

/* we have not come to the next header line */
/* error: qual string is of a different length */

