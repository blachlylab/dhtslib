/// @file hts_endian.h
/// Byte swapping and unaligned access functions.
/*
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

import core.stdc.config;

extern (C):

/*
 * Compile-time endianness tests.
 *
 * Note that these tests may fail.  They should only be used to enable
 * faster versions of endian-neutral implementations.  The endian-neutral
 * version should always be available as a fall-back.
 *
 * See https://sourceforge.net/p/predef/wiki/Endianness/
 */

/* Save typing as both endian and unaligned tests want to know about x86 */ /* x86 and x86_64 platform */

/** @def HTS_LITTLE_ENDIAN
 *  @brief Defined if platform is known to be little-endian
 */

/** @def HTS_BIG_ENDIAN
 *  @brief Defined if platform is known to be big-endian
 */

/** @def HTS_ENDIAN_NEUTRAL
 *  @brief Define this to disable any endian-specific optimizations
 */

/* Disable all endian-specific code. */

/** @def HTS_ALLOW_UNALIGNED
 *  @brief Control use of unaligned memory access.
 *
 * Defining HTS_ALLOW_UNALIGNED=1 converts shift-and-or to simple casts on
 * little-endian platforms that can tolerate unaligned access (notably Intel
 * x86).
 *
 * Defining HTS_ALLOW_UNALIGNED=0 forces shift-and-or.
 */

// Consider using AX_CHECK_ALIGNED_ACCESS_REQUIRED in autoconf.

enum HTS_ALLOW_UNALIGNED = 1;

// This prevents problems with gcc's vectoriser generating the wrong
// instructions for unaligned data.

alias uint16_u = ushort;
alias uint32_u = uint;
alias uint64_u = c_ulong;

/// Get a uint16_t value from an unsigned byte array
/** @param buf Pointer to source byte, may be unaligned
 *  @return A 16 bit unsigned integer
 *  The input is read in little-endian byte order.
 */
ushort le_to_u16 (const(ubyte)* buf);

/// Get a uint32_t value from an unsigned byte array
/** @param buf Pointer to source byte array, may be unaligned
 *  @return A 32 bit unsigned integer
 *  The input is read in little-endian byte order.
 */
uint le_to_u32 (const(ubyte)* buf);

/// Get a uint64_t value from an unsigned byte array
/** @param buf Pointer to source byte array, may be unaligned
 *  @return A 64 bit unsigned integer
 *  The input is read in little-endian byte order.
 */
ulong le_to_u64 (const(ubyte)* buf);

/// Store a uint16_t value in little-endian byte order
/** @param val The value to store
 *  @param buf Where to store it (may be unaligned)
 */
void u16_to_le (ushort val, ubyte* buf);

/// Store a uint32_t value in little-endian byte order
/** @param val The value to store
 *  @param buf Where to store it (may be unaligned)
 */
void u32_to_le (uint val, ubyte* buf);

/// Store a uint64_t value in little-endian byte order
/** @param val The value to store
 *  @param buf Where to store it (may be unaligned)
 */
void u64_to_le (ulong val, ubyte* buf);

/* Signed values.  Grab the data as unsigned, then convert to signed without
 * triggering undefined behaviour.  On any sensible platform, the conversion
 * should optimise away to nothing.
 */

/// Get an int8_t value from an unsigned byte array
/** @param buf Pointer to source byte array, may be unaligned
 *  @return A 8 bit signed integer
 *  The input data is interpreted as 2's complement representation.
 */
byte le_to_i8 (const(ubyte)* buf);

/// Get an int16_t value from an unsigned byte array
/** @param buf Pointer to source byte array, may be unaligned
 *  @return A 16 bit signed integer
 *  The input data is interpreted as 2's complement representation in
 *  little-endian byte order.
 */
short le_to_i16 (const(ubyte)* buf);

/// Get an int32_t value from an unsigned byte array
/** @param buf Pointer to source byte array, may be unaligned
 *  @return A 32 bit signed integer
 *  The input data is interpreted as 2's complement representation in
 *  little-endian byte order.
 */
int le_to_i32 (const(ubyte)* buf);

/// Get an int64_t value from an unsigned byte array
/** @param buf Pointer to source byte array, may be unaligned
 *  @return A 64 bit signed integer
 *  The input data is interpreted as 2's complement representation in
 *  little-endian byte order.
 */
long le_to_i64 (const(ubyte)* buf);

// Converting the other way is easier as signed -> unsigned is well defined.

/// Store a uint16_t value in little-endian byte order
/** @param val The value to store
 *  @param buf Where to store it (may be unaligned)
 */
void i16_to_le (short val, ubyte* buf);

/// Store a uint32_t value in little-endian byte order
/** @param val The value to store
 *  @param buf Where to store it (may be unaligned)
 */
void i32_to_le (int val, ubyte* buf);

/// Store a uint64_t value in little-endian byte order
/** @param val The value to store
 *  @param buf Where to store it (may be unaligned)
 */
void i64_to_le (long val, ubyte* buf);

/* Floating point.  Assumptions:
 *  Platform uses IEEE 754 format
 *  sizeof(float) == sizeof(uint32_t)
 *  sizeof(double) == sizeof(uint64_t)
 *  Endian-ness is the same for both floating point and integer
 *  Type-punning via a union is allowed
 */

/// Get a float value from an unsigned byte array
/** @param buf Pointer to source byte array, may be unaligned
 *  @return A 32 bit floating point value
 *  The input is interpreted as an IEEE 754 format float in little-endian
 *  byte order.
 */
float le_to_float (const(ubyte)* buf);

/// Get a double value from an unsigned byte array
/** @param buf Pointer to source byte array, may be unaligned
 *  @return A 64 bit floating point value
 *  The input is interpreted as an IEEE 754 format double in little-endian
 *  byte order.
 */
double le_to_double (const(ubyte)* buf);

/// Store a float value in little-endian byte order
/** @param val The value to store
 *  @param buf Where to store it (may be unaligned)
 */
void float_to_le (float val, ubyte* buf);

/// Store a double value in little-endian byte order
/** @param val The value to store
 *  @param buf Where to store it (may be unaligned)
 */
void double_to_le (double val, ubyte* buf);

/* HTS_ENDIAN_H */
