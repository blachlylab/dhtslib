/*  test/hts_endian.c -- hts_endian.h unit tests

    Copyright (C) 2017 Genome Research Ltd.

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

module htslib.test.hts_endian;

/**
* Adapted from htslib/test/hts_endian.c. License is included above.
*/

import core.stdc.stdio;
import core.stdc.stdlib;
import core.stdc.string;
import core.stdc.stdint;
import core.stdc.inttypes;

import std.conv : to;

import htslib.hts_endian;
import htslib.hts_log;

struct Test16
{
    uint8_t[2] u8;
    uint8_t[3] u8_unaligned;
    int16_t  i16;
    uint16_t u16;
}

struct Test32
{
    uint8_t[4] u8;
    uint8_t[5] u8_unaligned;
    int32_t  i32;
    uint32_t u32;
}

struct Test64
{
    uint8_t[8] u8;
    uint8_t[9] u8_unaligned;
    int64_t  i64;
    uint64_t u64;
}

struct Test_float{
    uint8_t[4] u8;
    uint8_t[5] u8_unaligned;
    float f;
}

struct Test_double
{
    uint8_t[8] u8;
    uint8_t[9] u8_unaligned;
    double d;
}

auto T16(T1, T2, T3, T4)(T1 b0, T2 b1, T3 sgn, T4 unsgn)
{ 
    return Test16([ b0, b1 ].to!(ubyte[2]), [ 0x00, b0, b1 ].to!(ubyte[3]), sgn.to!short, unsgn.to!ushort);
}

Test16[] tests_16_bit = [
    T16(0x00, 0x00,      0,     0),
    T16(0x01, 0x00,      1,     1),
    T16(0x00, 0x01,    256,   256),
    T16(0xff, 0x7f,  32767, 32767),
    T16(0x00, 0x80, -32768, 32768),
    T16(0xff, 0xff,     -1, 65535),
];

auto T32(T1, T2, T3, T4, T5, T6)(T1 b0, T2 b1, T3 b2, T4 b3, T5 sgn, T6 unsgn)
{ 
    return Test32([ b0, b1, b2, b3 ].to!(ubyte[4]), [ 0x00, b0, b1, b2, b3 ].to!(ubyte[5]), sgn.to!int, unsgn.to!uint);
}

Test32[] tests_32_bit = [
    T32(0x00, 0x00, 0x00, 0x00,           0,              0),
    T32(0x01, 0x00, 0x00, 0x00,           1,              1),
    T32(0x00, 0x01, 0x00, 0x00,         256,            256),
    T32(0x00, 0x00, 0x01, 0x00,       65536,          65536),
    T32(0xff, 0xff, 0xff, 0x7f,  2147483647,     2147483647),
    // Odd coding of signed result below avoids a compiler warning
    // as 2147483648 can't fit in a signed 32-bit number
    T32(0x00, 0x00, 0x00, 0x80, -2147483647 - 1, 2147483648U),
    T32(0xff, 0xff, 0xff, 0xff,          -1,     4294967295U),
];

auto T64(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10)(T1 b0, T2 b1, T3 b2, T4 b3, T5 b4, T6 b5, T7 b6, T8 b7, T9 sgn, T10 unsgn) 
{
    return Test64([ b0, b1, b2, b3, b4, b5, b6, b7].to!(ubyte[8]), [ 0x00, b0, b1, b2, b3, b4, b5, b6, b7].to!(ubyte[9]), sgn.to!long, unsgn.to!ulong);
}


Test64[] tests_64_bit = [
    T64(0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0, 0),
    T64(0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 1, 1),
    T64(0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 256, 256),
    T64(0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 65536, 65536),
    T64(0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 4294967296L, 4294967296UL),
    T64(0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f,
        9223372036854775807L, 9223372036854775807UL),
    // Odd coding of signed result below avoids a compiler warning
    // as 9223372036854775808 can't fit in a signed 64-bit number
    T64(0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80,
        -9223372036854775807L - 1L, 9223372036854775808UL),
    T64(0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        -1, 18446744073709551615UL),
];

auto TF(T1, T2, T3, T4, T5)(T1 b0, T2 b1, T3 b2, T4 b3, T5 f)
{
    return Test_float([b0, b1, b2, b3 ].to!(ubyte[4]),[ 0x00, b0, b1, b2, b3].to!(ubyte[5]), f);
}

Test_float[] tests_float = [
    TF(0x00, 0x00, 0x00, 0x00,  0.0f),
    TF(0x00, 0x00, 0x80, 0x3f,  1.0f),
    TF(0x00, 0x00, 0x80, 0xbf, -1.0f),
    TF(0x00, 0x00, 0x20, 0x41, 10.0f),
    TF(0xd0, 0x0f, 0x49, 0x40,  3.14159f),
    TF(0xa8, 0x0a, 0xff, 0x66,  6.022e23f),
    TF(0xcd, 0x84, 0x03, 0x13,  1.66e-27f),
];

auto TD(T1, T2, T3, T4, T5, T6, T7, T8, T9)(T1 b0, T2 b1, T3 b2, T4 b3, T5 b4, T6 b5, T7 b6, T8 b7, T9 d)
{
    return Test_double([b0, b1, b2, b3, b4, b5, b6, b7].to!(ubyte[8]), [0x00, b0, b1, b2, b3, b4, b5, b6, b7 ].to!(ubyte[9]), d);
}

Test_double[] tests_double = [
    TD(0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,  0.0),
    TD(0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f,  1.0),
    TD(0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0xbf, -1.0),
    TD(0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x24, 0x40, 10.0),
    TD(0x18, 0x2d, 0x44, 0x54, 0xfb, 0x21, 0x09, 0x40,  3.141592653589793),
    TD(0x2b, 0x08, 0x0c, 0xd3, 0x85, 0xe1, 0xdf, 0x44,  6.022140858e23),
    TD(0x55, 0xfa, 0x81, 0x74, 0xf7, 0x71, 0x60, 0x3a,  1.66053904e-27),
];

auto NELE(T)(T x){return x.length;}

// char * to_hex(uint8_t *buf.ptr, int len) {
//     char[64] str;
//     int i, o;

//     for (i = 0, o = 0; i < len; i++, o += 3) {
//         snprintf(str + o, sizeof(str) - o, "%02x ", buf.ptr[i]);
//     }
//     return str;
// }

debug(dhtslib_unittest) unittest
{
    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
    hts_log_info(__FUNCTION__,"hts_endian tests: 16 bit");
    uint8_t[9] buf;
    size_t i;
    int errors = 0;

    for (i = 0; i < NELE(tests_16_bit); i++) {
        uint16_t u16;
        int16_t  i16;

        // if (verbose) {
        //     fprintf(stderr, "%s %6"~PRId16~" %6"~PRId16~"\n",
        //             to_hex(tests_16_bit[i].u8.ptr, 2),
        //             tests_16_bit[i].i16, tests_16_bit[i].u16);
        // }

        u16 = le_to_u16(tests_16_bit[i].u8.ptr);
        assert(u16 == tests_16_bit[i].u16);

        i16 = le_to_i16(tests_16_bit[i].u8.ptr);
        assert(i16 == tests_16_bit[i].i16);

        u16 = le_to_u16(tests_16_bit[i].u8_unaligned.ptr + 1);
        assert(u16 == tests_16_bit[i].u16);

        i16 = le_to_i16(tests_16_bit[i].u8_unaligned.ptr + 1);
        assert(i16 == tests_16_bit[i].i16);

        u16_to_le(tests_16_bit[i].u16, buf.ptr);
        assert(memcmp(buf.ptr, tests_16_bit[i].u8.ptr, 2) == 0);

        i16_to_le(tests_16_bit[i].i16, buf.ptr);
        assert(memcmp(buf.ptr, tests_16_bit[i].u8.ptr, 2) == 0);

        u16_to_le(tests_16_bit[i].u16, buf.ptr + 1);
        assert(memcmp(buf.ptr + 1, tests_16_bit[i].u8.ptr, 2) == 0);

        i16_to_le(tests_16_bit[i].i16, buf.ptr + 1);
        assert(memcmp(buf.ptr + 1, tests_16_bit[i].u8.ptr, 2) == 0);
    }
}

debug(dhtslib_unittest) unittest
{
    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
    hts_log_info(__FUNCTION__,"hts_endian tests: 32 bit");
    uint8_t[9] buf;
    size_t i;
    int errors = 0;

    for (i = 0; i < NELE(tests_32_bit); i++) {
        uint32_t u32;
        int32_t  i32;

        // if (verbose) {
        //     fprintf(stderr, "%s %11"~PRId32~" %11"~PRIu32~"\n",
        //             to_hex(tests_32_bit[i].u8.ptr, 4),
        //             tests_32_bit[i].i32, tests_32_bit[i].u32);
        // }

        u32 = le_to_u32(tests_32_bit[i].u8.ptr);
        assert(u32 == tests_32_bit[i].u32);
        i32 = le_to_i32(tests_32_bit[i].u8.ptr);
        assert(i32 == tests_32_bit[i].i32);

        u32 = le_to_u32(tests_32_bit[i].u8_unaligned.ptr + 1);
        assert(u32 == tests_32_bit[i].u32);

        i32 = le_to_i32(tests_32_bit[i].u8_unaligned.ptr + 1);
        assert(i32 == tests_32_bit[i].i32);

        u32_to_le(tests_32_bit[i].u32, buf.ptr);
        assert(memcmp(buf.ptr, tests_32_bit[i].u8.ptr, 4) == 0);

        i32_to_le(tests_32_bit[i].i32, buf.ptr);
        assert(memcmp(buf.ptr, tests_32_bit[i].u8.ptr, 4) == 0);

        u32_to_le(tests_32_bit[i].u32, buf.ptr + 1);
        assert(memcmp(buf.ptr + 1, tests_32_bit[i].u8.ptr, 4) == 0);

        i32_to_le(tests_32_bit[i].i32, buf.ptr + 1);
        assert(memcmp(buf.ptr + 1, tests_32_bit[i].u8.ptr, 4) == 0);
    }
}

debug(dhtslib_unittest) unittest
{
    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
    hts_log_info(__FUNCTION__,"hts_endian tests: 64 bit");
    uint8_t[9] buf;
    size_t i;
    int errors = 0;

    for (i = 0; i < NELE(tests_64_bit); i++) {
        uint64_t u64;
        int64_t  i64;

        // if (verbose) {
        //     fprintf(stderr, "%s %20"~PRId64~" %20"~PRIu64~"\n",
        //             to_hex(tests_64_bit[i].u8.ptr, 8),
        //             tests_64_bit[i].i64, tests_64_bit[i].u64);
        // }

        u64 = le_to_u64(tests_64_bit[i].u8.ptr);
        assert(u64 == tests_64_bit[i].u64);

        i64 = le_to_i64(tests_64_bit[i].u8.ptr);
        assert(i64 == tests_64_bit[i].i64);

        u64 = le_to_u64(tests_64_bit[i].u8_unaligned.ptr + 1);
        assert(u64 == tests_64_bit[i].u64);

        i64 = le_to_i64(tests_64_bit[i].u8_unaligned.ptr + 1);
        assert(i64 == tests_64_bit[i].i64);

        u64_to_le(tests_64_bit[i].u64, buf.ptr);
        assert(memcmp(buf.ptr, tests_64_bit[i].u8.ptr, 8) == 0);

        i64_to_le(tests_64_bit[i].i64, buf.ptr);
        assert(memcmp(buf.ptr, tests_64_bit[i].u8.ptr, 8) == 0);

        u64_to_le(tests_64_bit[i].u64, buf.ptr + 1);
        assert(memcmp(buf.ptr + 1, tests_64_bit[i].u8.ptr, 8) == 0);

        i64_to_le(tests_64_bit[i].i64, buf.ptr + 1);
        assert(memcmp(buf.ptr + 1, tests_64_bit[i].u8.ptr, 8) == 0);
    }
}

debug(dhtslib_unittest) unittest
{
    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
    hts_log_info(__FUNCTION__,"hts_endian tests: float");
    uint8_t[9] buf;
    size_t i;
    int errors = 0;

    for (i = 0; i < NELE(tests_float); i++) {
        float f;

        // if (verbose) {
        //     fprintf(stderr, "%s %g\n",
        //             to_hex(tests_float[i].u8.ptr, 4), tests_float[i].f);
        // }

        f = le_to_float(tests_float[i].u8.ptr);
        assert(f == tests_float[i].f);

        f = le_to_float(tests_float[i].u8_unaligned.ptr + 1);
        assert(f == tests_float[i].f);

        float_to_le(tests_float[i].f, buf.ptr);
        assert(memcmp(tests_float[i].u8.ptr, buf.ptr, 4) == 0);

        float_to_le(tests_float[i].f, buf.ptr + 1);
        assert(memcmp(tests_float[i].u8.ptr, buf.ptr + 1, 4) == 0);
    }
}

debug(dhtslib_unittest) unittest
{
    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
    hts_log_info(__FUNCTION__,"hts_endian tests: double bit");
    uint8_t[9] buf;
    size_t i;
    int errors = 0;

    for (i = 0; i < NELE(tests_double); i++) {
        double f;

        // if (verbose) {
        //     fprintf(stderr, "%s %.15g\n",
        //             to_hex(tests_double[i].u8.ptr, 8), tests_double[i].d);
        // }

        f = le_to_double(tests_double[i].u8.ptr);
        assert(f == tests_double[i].d);

        f = le_to_double(tests_double[i].u8_unaligned.ptr + 1);
        assert(f == tests_double[i].d);

        double_to_le(tests_double[i].d, buf.ptr);
        assert(memcmp(tests_double[i].u8.ptr, buf.ptr, 8) == 0);

        double_to_le(tests_double[i].d, buf.ptr + 1);
        assert(memcmp(tests_double[i].u8.ptr, buf.ptr + 1, 8) == 0);
    }
}

// int main(int argc, char **argv) {
//     int verbose = 0;
//     int errors = 0;

//     if (argc > 1 && strcmp(argv[1], "-v") == 0) verbose = 1;

//     errors += t16_bit(verbose);
//     errors += t32_bit(verbose);
//     errors += t64_bit(verbose);
//     errors += t_float(verbose);
//     errors += t_double(verbose);
//     if (errors) {
//         fprintf(stderr, "%d errors\n", errors);
//         return EXIT_FAILURE;
//     }

//     return EXIT_SUCCESS;
// }
