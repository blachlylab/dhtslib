/* The MIT License

   Copyright (C) 2011 by Attractive Chaos <attractor@live.co.uk>
   Copyright (C) 2013-2014, 2016, 2018-2020 Genome Research Ltd.

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

import core.stdc.config;
import core.stdc.stdio;
import core.stdc.stdlib;

extern (C):

// __MINGW_PRINTF_FORMAT

enum HAVE___BUILTIN_CLZ = 1;

/* kstring_t is a simple non-opaque type whose fields are likely to be
 * used directly by user code (but see also ks_str() and ks_len() below).
 * A kstring_t object is initialised by either of
 *       kstring_t str = KS_INITIALIZE;
 *       kstring_t str; ...; ks_initialize(&str);
 * and either ownership of the underlying buffer should be given away before
 * the object disappears (see ks_release() below) or the kstring_t should be
 * destroyed with  ks_free(&str) or free(str.s) */

alias KSTRING_T = kstring_t;

struct kstring_t
{
    size_t l;
    size_t m;
    char* s;
}

struct ks_tokaux_t
{
    ulong[4] tab;
    int sep;
    int finished;
    const(char)* p; // end of the current token
}

int kvsprintf (kstring_t* s, const(char)* fmt, va_list ap);

int ksprintf (kstring_t* s, const(char)* fmt, ...);

int kputd (double d, kstring_t* s); // custom %g only handler

int ksplit_core (char* s, int delimiter, int* _max, int** _offsets);

char* kstrstr (const(char)* str, const(char)* pat, int** _prep);

char* kstrnstr (const(char)* str, const(char)* pat, int n, int** _prep);

void* kmemmem (const(void)* _str, int n, const(void)* _pat, int m, int** _prep);

/* kstrtok() is similar to strtok_r() except that str is not
	 * modified and both str and sep can be NULL. For efficiency, it is
	 * actually recommended to set both to NULL in the subsequent calls
	 * if sep is not changed. */
char* kstrtok (const(char)* str, const(char)* sep, ks_tokaux_t* aux);

/* kgetline() uses the supplied fgets()-like function to read a "\n"-
 * or "\r\n"-terminated line from fp.  The line read is appended to the
 * kstring without its terminator and 0 is returned; EOF is returned at
 * EOF or on error (determined by querying fp, as per fgets()). */
alias kgets_func = char* function (char*, int, void*);
int kgetline (kstring_t* s, char* function () fgets_fn, void* fp);

/* kgetline2() uses the supplied hgetln()-like function to read a "\n"-
 * or "\r\n"-terminated line from fp.  The line read is appended to the
 * ksring without its terminator and 0 is returned; EOF is returned at
 * EOF or on error (determined by querying fp, as per fgets()). */
alias kgets_func2 = c_long function (char*, size_t, void*);
int kgetline2 (kstring_t* s, ssize_t function () fgets_fn, void* fp);

/// kstring initializer for structure assignment

/// kstring initializer for pointers
/**
   @note Not to be used if the buffer has been allocated.  Use ks_release()
   or ks_clear() instead.
*/

void ks_initialize (kstring_t* s);

/// Resize a kstring to a given capacity
int ks_resize (kstring_t* s, size_t size);

/// Increase kstring capacity by a given number of bytes

// Overflow check
int ks_expand (kstring_t* s, size_t expansion);

/// Returns the kstring buffer
char* ks_str (kstring_t* s);

/// Returns the kstring buffer, or an empty string if l == 0
/**
 * Unlike ks_str(), this function will never return NULL.  If the kstring is
 * empty it will return a read-only empty string.  As the returned value
 * may be read-only, the caller should not attempt to modify it.
 */
const(char)* ks_c_str (kstring_t* s);

size_t ks_len (kstring_t* s);

/// Reset kstring length to zero
/**
   @return The kstring itself

   Example use: kputsn(string, len, ks_clear(s))
*/
kstring_t* ks_clear (kstring_t* s);

// Give ownership of the underlying buffer away to something else (making
// that something else responsible for freeing it), leaving the kstring_t
// empty and ready to be used again, or ready to go out of scope without
// needing  free(str.s)  to prevent a memory leak.
char* ks_release (kstring_t* s);

/// Safely free the underlying buffer in a kstring.
void ks_free (kstring_t* s);

int kputsn (const(char)* p, size_t l, kstring_t* s);

int kputs (const(char)* p, kstring_t* s);

int kputc (int c, kstring_t* s);

int kputc_ (int c, kstring_t* s);

int kputsn_ (const(void)* p, size_t l, kstring_t* s);

// Trivial case - also prevents __builtin_clz(0), which is undefined

// Find out how many digits are to be printed.

/*
 * Table method - should be quick if clz can be done in hardware.
 * Find the most significant bit of the value to print and look
 * up in a table to find out how many decimal digits are needed.
 * This number needs to be adjusted by 1 for cases where the decimal
 * length could vary for a given number of bits (for example,
 * a four bit number could be between 8 and 15).
 */

// Fallback for when clz is not available

// Add digits two at a time

// Last one (if necessary).  We know that x < 10 by now.
int kputuw (uint x, kstring_t* s);

int kputw (int c, kstring_t* s);

int kputll (long c, kstring_t* s);

int kputl (c_long c, kstring_t* s);

/*
 * Returns 's' split by delimiter, with *n being the number of components;
 *         NULL on failure.
 */
int* ksplit (kstring_t* s, int delimiter, int* n);

