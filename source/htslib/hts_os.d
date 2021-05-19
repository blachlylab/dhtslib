/// @file hts_os.h
/// Operating System specific tweaks, for compatibility with POSIX.
/*
   Copyright (C) 2017, 2019-2020 Genome Research Ltd.

    Author: James Bonfield <jkb@sanger.ac.uk>

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
module htslib.hts_os;

import core.stdc.config;

extern (C):

/* This is srand48_deterministic() on platforms that provide it, or srand48()
   otherwise (or our own POSIX srand48() on platforms that provide neither).
   Hence calling hts_srand48() will always set up the same POSIX-determined
   sequence of pseudo-random numbers on any platform, while calling srand48()
   may (e.g., on OpenBSD) set up a different non-deterministic sequence. */
void hts_srand48 (c_long seed);

double hts_erand48 (ref ushort[3] xseed);

double hts_drand48 ();

c_long hts_lrand48 ();

// Windows usually lacks *rand48(), but cygwin provides them.

/* def _WIN32 - disabled for now, not currently used */
/* Check if the fd is a cygwin/msys's pty. */

// HTSLIB_HTS_OS_H
