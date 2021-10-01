/*  test_kstring.c -- kstring unit tests

    Copyright (C) 2018, 2020 Genome Research Ltd.

    Author: Rob Davies <rmd@sanger.ac.uk>

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
module htslib.test.kstring;

/**
* Adapted from htslib/test/test_kstring.c. License is included above.
*/

// #include <config.h>

import core.stdc.stdio;
import core.stdc.stdlib;
import core.stdc.string;
// #include <stdio.h>
// #include <limits.h>
// #include <stdint.h>
// #include <inttypes.h>
// #include <unistd.h>

// #include "../htslib/kstring.h"
import htslib.kstring;
import htslib.kroundup;

nothrow @nogc @system:

pragma(inline, true)
void clamp(long *val, long min, long max) {
    if (*val < min) *val = min;
    if (*val > max) *val = max;
}

int test_kroundup_size_t(int verbose) {
    size_t val, exp;
    int ret = 0;

    val = 0;
    kroundup_size_t(val);
    if (verbose) {
        printf("kroundup_size_t(0) = 0x%zx\n", val);
    }
    if (val != 0) {
        fprintf(stderr, "kroundup_size_t(0) produced 0x%zx, expected 0\n", val);
        ret = -1;
    }

    for (exp = 0; exp < val.sizeof * 8; exp++) {
        size_t expected = (cast(size_t) 1) << exp;
        ssize_t delta;
        for (delta = exp > 1 ? -1 : 0; delta <= (exp < 2 ? 0 : 1); delta++) {
            size_t val_in = expected + delta;
            val = val_in;
            kroundup_size_t(val);
            if (verbose) {
                printf("kroundup_size_t(0x%zx) = 0x%zx\n", val_in, val);
            }
            if (delta <= 0) {
                if (val != expected) {
                    fprintf(stderr, "kroundup_size_t(0x%zx) produced 0x%zx, expected 0x%zx\n",
                            val_in, val, expected);
                    ret = -1;
                }
            } else {
                expected *= 2;
                if (!expected) --expected;
                if (val != expected) {
                    fprintf(stderr, "kroundup_size_t(0x%zx) produced 0x%zx, expected 0x%zx\n",
                            val_in, val, expected);
                    ret = -1;
                }
            }
        }
    }
    return ret;
}

int test_kroundup_signed(int verbose) {
    int val, ret = 0;
    size_t exp;
    for (exp = 0; exp < val.sizeof * 8 - 1; exp++) {
        uint expected = (cast(uint) 1) << exp;
        ssize_t delta;
        for (delta = exp > 1 ? -1 : 0; delta <= (exp < 2 ? 0 : 1); delta++) {
            int val_in = cast(int)(expected + delta);
            val = val_in;
            kroundup32(val);
            if (verbose) {
                printf("kroundup32(%d) = %d\n", val_in, val);
            }
            if (delta <= 0) {
                if (cast(uint) val != expected) {
                    fprintf(stderr, "kroundup32(%d) produced %d, expected %u\n",
                            val_in, val, expected);
                    ret = -1;
                }
            } else {
                if (exp < val.sizeof * 8 - 2) {
                    expected *= 2;
                } else {
                    expected = ((expected - 1) << 1 | 1);
                }
                if (cast(uint) val != expected) {
                    fprintf(stderr, "kroundup32(%d) produced %d, expected %u\n",
                            val_in, val, expected);
                    ret = -1;
                }
            }
        }
    }
    return ret;
}

int test_kputuw_from_to(T)(kstring_t *str, T s, T e) {
    uint i = cast(uint)s;

    for (;;) {
        str.l = 0;
        memset(str.s, 0xff, str.m);
        if (kputuw(i, str) < 0 || !str.s) {
            perror("kputuw");
            return -1;
        }
        if (str.l >= str.m || str.s[str.l] != '\0') {
            fprintf(stderr, "No NUL termination on string from kputuw\n");
            return -1;
        }
        if (i != strtoul(str.s, null, 10)) {
            fprintf(stderr,
                    "kputuw wrote the wrong value, expected %u, got %s\n",
                    i, str.s);
            return -1;
        }
        if (i >= e) break;
        i++;
    }
    return 0;
}

int test_kputuw(T)(T start, T end) {
    kstring_t str;
    ks_initialize(&str);
    long val;

    str.s = cast(char*)malloc(2);
    if (!str.s) {
        perror("malloc");
        return -1;
    }
    str.m = 2;

    for (val = 0; val < uint.max; val = val == 0 ? 1 : val * 10) {
        auto s = val == 0 ? 0 : val - 5;
        auto e = val + 5;

        if (test_kputuw_from_to(&str, s, e) < 0) {
            free(ks_release(&str));
            return -1;
        }
    }

    if (test_kputuw_from_to(&str, uint.max - 5, uint.max) < 0) {
        free(ks_release(&str));
        return -1;
    }

    str.m = 1; // Force a resize
    clamp(&start, 0, uint.max);
    clamp(&end,   0, uint.max);

    if (test_kputuw_from_to(&str, start, end) < 0) {
        free(ks_release(&str));
        return -1;
    }

    free(ks_release(&str));

    return 0;
}

int test_kputw_from_to(T)(kstring_t *str, T s, T e) {
    int i = cast(int)s;

    for (;;) {
        str.l = 0;
        memset(str.s, 0xff, str.m);
        if (kputw(i, str) < 0 || !str.s) {
            perror("kputw");
            return -1;
        }
        if (str.l >= str.m || str.s[str.l] != '\0') {
            fprintf(stderr, "No NUL termination on string from kputw\n");
            return -1;
        }
        if (i != strtol(str.s, null, 10)) {
            fprintf(stderr,
                    "kputw wrote the wrong value, expected %d, got %s\n",
                    i, str.s);
            return -1;
        }
        if (i >= e) break;
        i++;
    }
    return 0;
}

int test_kputw(long start, long end) {
    kstring_t str = { 0, 0, null };
    long val;

    str.s = cast(char*)malloc(2);
    if (!str.s) {
        perror("malloc");
        return -1;
    }
    str.m = 2;

    for (val = 1; val < int.max; val *= 10) {
        if (test_kputw_from_to(&str, val > 5 ? val - 5 : 0, val + 5) < 0) {
            free(ks_release(&str));
            return -1;
        }
    }

    for (val = -1; val > int.min; val *= 10) {
        if (test_kputw_from_to(&str, val - 5, val < -5 ? val + 5 : 0) < 0) {
            free(ks_release(&str));
            return -1;
        }
    }

    if (test_kputw_from_to(&str, int.max - 5, int.max) < 0) {
        free(ks_release(&str));
        return -1;
    }

    if (test_kputw_from_to(&str, int.min, int.min + 5) < 0) {
        free(ks_release(&str));
        return -1;
    }

    str.m = 1; // Force a resize
    clamp(&start, int.min, int.max);
    clamp(&end,   int.min, int.max);

    if (test_kputw_from_to(&str, start, end) < 0) {
        free(ks_release(&str));
        return -1;
    }

    free(ks_release(&str));

    return 0;
}

// callback used by test_kgetline
extern(C) char *mock_fgets(char *str, int num, void *p) {
    int *mock_state = cast(int*)p;
    (*mock_state)++;
    switch (*mock_state) {
    case 1:
    case 4:
    case 7:
        // a few characters, no endline
        strcpy(str, "ABCD");
        break;
    case 2:
    case 3:
        // \n endline
        strcpy(str, "\n");
        break;
    case 5:
    case 6:
        // \r\n endline
        strcpy(str, "\r\n");
        break;
    default:
        // eof
        return null;
    }

    return str;
}

unittest
{
    kstring_t s;
    int mock_state = 0;

    // normal line, \n terminated, called with non-empty s
    kputs("_", &s);
    assert(0 == kgetline(&s, cast(kgets_func)&mock_fgets, cast(void *)&mock_state) || 0 != strcmp("_ABCD", s.s) || 5 != s.l);
    s.l = 0;
    // empty line, \n terminated
    assert(0 == kgetline(&s, cast(kgets_func)&mock_fgets, cast(void *)&mock_state) || 0 != strcmp("", s.s) || 0 != s.l);
    s.l = 0;
    // normal line, \r\n terminated
    assert(0 == kgetline(&s, cast(kgets_func)&mock_fgets, cast(void *)&mock_state) || 0 != strcmp("ABCD", s.s) || 4 != s.l);
    s.l = 0;
    // empty line, \r\n terminated
    assert(0 == kgetline(&s, cast(kgets_func)&mock_fgets, cast(void *)&mock_state) || 0 != strcmp("", s.s) || 0 != s.l);
    s.l = 0;
    // line terminated by EOF
    assert(0 == kgetline(&s, cast(kgets_func)&mock_fgets, cast(void *)&mock_state) || 0 != strcmp("ABCD", s.s) || 4 != s.l);
    s.l = 0;
    // EOF
    assert(EOF == kgetline(&s, cast(kgets_func)&mock_fgets, cast(void *)&mock_state) || 0 != s.l);

    ks_free(&s);
}

// callback used by test_kgetline2
extern(C) ssize_t mock_fgets2(char *str, size_t num, void *p) {
    int *mock_state = cast(int*)p;
    (*mock_state)++;
    switch (*mock_state) {
    case 1:
    case 4:
    case 7:
        // a few characters, no endline
        strcpy(str, "ABCD");
        break;
    case 2:
    case 3:
        // \n endline
        strcpy(str, "\n");
        break;
    case 5:
    case 6:
        // \r\n endline
        strcpy(str, "\r\n");
        break;
    default:
        // eof
        return 0;
    }

    return strlen(str);
}

unittest 
{
    kstring_t s;
    int mock_state = 0;

    // normal line, \n terminated, called with non-empty s
    kputs("_", &s);
    assert(0 == kgetline2(&s, cast(kgets_func2)&mock_fgets2, cast(void *)&mock_state) || 0 != strcmp("_ABCD", s.s) || 5 != s.l);
    s.l = 0;
    // empty line, \n terminated
    assert(0 == kgetline2(&s, cast(kgets_func2)&mock_fgets2, cast(void *)&mock_state) || 0 != strcmp("", s.s) || 0 != s.l);
    s.l = 0;
    // normal line, \r\n terminated
    assert(0 == kgetline2(&s, cast(kgets_func2)&mock_fgets2, cast(void *)&mock_state) || 0 != strcmp("ABCD", s.s) || 4 != s.l);
    s.l = 0;
    // empty line, \r\n terminated
    assert(0 == kgetline2(&s, cast(kgets_func2)&mock_fgets2, cast(void *)&mock_state) || 0 != strcmp("", s.s) || 0 != s.l);
    s.l = 0;
    // line terminated by EOF
    assert(0 == kgetline2(&s, cast(kgets_func2)&mock_fgets2, cast(void *)&mock_state) || 0 != strcmp("ABCD", s.s) || 4 != s.l);
    s.l = 0;
    // EOF
    assert(EOF == kgetline2(&s, cast(kgets_func2)&mock_fgets2, cast(void *)&mock_state) || 0 != s.l);

    ks_free(&s);
}

unittest
{
    int verbose = 0;
    long start = 0;
    long end = 0;
    assert(test_kroundup_size_t(verbose) == 0);
    assert(test_kroundup_signed(verbose) == 0);
    assert(test_kputuw(start, end) == 0);
    assert(test_kputw(start, end) == 0);
}