/* The MIT License

   Copyright (c) 2008, 2012-2013, 2017-2019 Genome Research Ltd (GRL).

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

/* Contact: Heng Li <lh3@sanger.ac.uk> */

/*
  2012-12-11 (0.1.4):

    * Defined __ks_insertsort_##name as static to compile with C99.

  2008-11-16 (0.1.4):

    * Fixed a bug in introsort() that happens in rare cases.

  2008-11-05 (0.1.3):

    * Fixed a bug in introsort() for complex comparisons.

	* Fixed a bug in mergesort(). The previous version is not stable.

  2008-09-15 (0.1.2):

	* Accelerated introsort. On my Mac (not on another Linux machine),
	  my implementation is as fast as the C++ standard library's sort()
	  on random input.

	* Added combsort and in introsort, switch to combsort if the
	  recursion is too deep.

  2008-09-13 (0.1.1):

	* Added k-small algorithm

  2008-09-05 (0.1.0):

	* Initial version

*/
module htslib.ksort;
import core.stdc.string;

extern (C):

/* klib_unused */

// Use our own drand48() symbol (used by ks_shuffle) to avoid portability
// problems on Windows.  Don't include htslib/hts_os.h for this as it
// may not get on with older attempts to fix this in code that includes
// this file.
double hts_drand48 ();

struct ks_isort_stack_t
{
    void* left;
    void* right;
    int depth;
}

/* This function is adapted from: http://ndevilla.free.fr/median/ */
/* 0 <= kk < n */

extern (D) auto ks_lt_generic(T0, T1)(auto ref T0 a, auto ref T1 b)
{
    return a < b;
}

extern (D) auto ks_lt_str(T0, T1)(auto ref T0 a, auto ref T1 b)
{
    return strcmp(a, b) < 0;
}

alias ksstr_t = const(char)*;

// enum KSORT_INIT_STR = KSORT_INIT(str, ksstr_t, ks_lt_str);

// enum KSORT_INIT_STATIC_STR = KSORT_INIT_STATIC(str, ksstr_t, ks_lt_str);

// enum KSORT_INIT2_STR = KSORT_INIT2(str, SCOPE, ksstr_t, ks_lt_str);

