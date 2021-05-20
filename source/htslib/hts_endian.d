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
module htslib.hts_endian;

import core.stdc.config;

/*
 * Compile-time endianness tests.
 *
 * Note that these tests may fail.  They should only be used to enable
 * faster versions of endian-neutral implementations.  The endian-neutral
 * version should always be available as a fall-back.
 *
 * See https://sourceforge.net/p/predef/wiki/Endianness/
 */

version(X86) enum HTS_x86 = true;
version(X86_64) enum HTS_x86 = true;
else enum HTS_x86 = false;

/* Save typing as both endian and unaligned tests want to know about x86 */ /* x86 and x86_64 platform */

/** @def HTS_LITTLE_ENDIAN
 *  @brief Defined if platform is known to be little-endian
 */

version(LittleEndian) enum HTS_LITTLE_ENDIAN = true;
else enum HTS_LITTLE_ENDIAN = false;

/** @def HTS_BIG_ENDIAN
 *  @brief Defined if platform is known to be big-endian
 */

version(BigEndian) enum HTS_BIG_ENDIAN = true;
else enum HTS_BIG_ENDIAN = false;

/** @def HTS_ENDIAN_NEUTRAL
 *  @brief Define this to disable any endian-specific optimizations
 */

/* Disable all endian-specific code. */
version(HTS_ENDIAN_NEUTRAL)
{
    enum HTS_LITTLE_ENDIAN = false;
    enum HTS_BIG_ENDIAN = false;
}

/** @def HTS_ALLOW_UNALIGNED
 *  @brief Control use of unaligned memory access.
 *
 * Defining HTS_ALLOW_UNALIGNED=1 converts shift-and-or to simple casts on
 * little-endian platforms that can tolerate unaligned access (notably Intel
 * x86).
 *
 * Defining HTS_ALLOW_UNALIGNED=0 forces shift-and-or.
 */

static if(HTS_x86) enum HTS_ALLOW_UNALIGNED = true;
else enum HTS_ALLOW_UNALIGNED = false;

// Consider using AX_CHECK_ALIGNED_ACCESS_REQUIRED in autoconf.

// This prevents problems with gcc's vectoriser generating the wrong
// instructions for unaligned data.

static if(HTS_ALLOW_UNALIGNED){
    alias uint16_u = align(1) ushort;
    alias uint32_u = align(1) uint;
    alias uint64_u = align(1) c_ulong;
}else{
    alias uint16_u = ushort;
    alias uint32_u = uint;
    alias uint64_u = c_ulong;    
}

/// Basically just a byte
/// This workaround to avoid int promotion issues
/// implmentation direvied from this post on the dlang forums
/// https://forum.dlang.org/post/r0imk7$14b7$3@digitalmars.com
struct int8_t
{
    byte _val;
    alias _val this;

    pragma(inline, true)
    this(byte x) @nogc
    {
        _val = x;
    }

    pragma(inline, true)
    int8_t opBinary(string op, T)(T other)
        if ((is(T == int8_t) || is(T == byte) || is(T == ubyte)))
    {
        return mixin("int8_t(cast(byte)(_val " ~ op ~ " cast(byte)other))");
    }

    pragma(inline, true)
    int8_t opBinaryRight(string op, T)(T other)
        if ((is(T == int8_t) || is(T == byte) || is(T == ubyte)))
    {
        return mixin("int8_t(cast(byte)(_val " ~ op ~ " cast(byte)other))");
    }

    pragma(inline, true)
    int8_t opUnary(string op: "-")()
    {
        return (cast(int8_t) _val ^ cast(ubyte) 0xFF) + cast(ubyte) 1;
    }
}

unittest
{
    assert(-int8_t(cast(byte) 8) == int8_t(cast(byte) -8));
}

/// Basically just a short
/// This workaround to avoid int promotion issues
/// implmentation direvied from this post on the dlang forums
/// https://forum.dlang.org/post/r0imk7$14b7$3@digitalmars.com
struct int16_t
{
    short _val;
    alias _val this;

    pragma(inline, true)
    this(short x) @nogc
    {
        _val = x;
    }

    pragma(inline, true)
    int16_t opBinary(string op, T)(T other)
        if ((is(T == int16_t) || is(T == short) || is(T == ushort)))
    {
        return mixin("int16_t(cast(short)(_val " ~ op ~ " cast(short)other))");
    }

    pragma(inline, true)
    int16_t opBinaryRight(string op, T)(T other)
        if ((is(T == int16_t) || is(T == short) || is(T == ushort)))
    {
        return mixin("int16_t(cast(short)(_val " ~ op ~ " cast(short)other))");
    }

    pragma(inline, true)
    int16_t opUnary(string op: "-")()
    {
        return (cast(int16_t) _val ^ cast(ushort) 0xFFFF) + cast(ushort) 1;
    }
}

unittest
{
    assert(-int16_t(cast(short) 8) == int16_t(cast(short) -8));
}

pragma(inline, true):
@nogc:
@system:
/// Get a ushort value from an unsigned byte array
/** @param buf Pointer to source byte, may be unaligned
 *  @return A 16 bit unsigned integer
 *  The input is read in little-endian byte order.
 */
ushort le_to_u16(const(ubyte)* buf)
{
    static if(HTS_LITTLE_ENDIAN && HTS_ALLOW_UNALIGNED)
        return *(cast(uint16_u *) buf);
    else
        return cast(ushort) buf[0] | (cast(ushort) buf[1] << 8);
}

/// Get a uint value from an unsigned byte array
/** @param buf Pointer to source byte array, may be unaligned
 *  @return A 32 bit unsigned integer
 *  The input is read in little-endian byte order.
 */
uint le_to_u32(const(ubyte)* buf)
{
    static if(HTS_LITTLE_ENDIAN && HTS_ALLOW_UNALIGNED)
        return *(cast(uint32_u *) buf);
    else
        return (cast(uint) buf[0] |
            (cast(uint) buf[1] << 8) |
            (cast(uint) buf[2] << 16) |
            (cast(uint) buf[3] << 24));
}

/// Get a ulong value from an unsigned byte array
/** @param buf Pointer to source byte array, may be unaligned
 *  @return A 64 bit unsigned integer
 *  The input is read in little-endian byte order.
 */
ulong le_to_u64(const(ubyte)* buf)
{
    static if(HTS_LITTLE_ENDIAN && HTS_ALLOW_UNALIGNED)
        return *(cast(uint64_u *) buf);
    else
        return (cast(ulong) buf[0] |
            (cast(ulong) buf[1] << 8) |
            (cast(ulong) buf[2] << 16) |
            (cast(ulong) buf[3] << 24) |
            (cast(ulong) buf[4] << 32) |
            (cast(ulong) buf[5] << 40) |
            (cast(ulong) buf[6] << 48) |
            (cast(ulong) buf[7] << 56));
}

/// Store a ushort value in little-endian byte order
/** @param val The value to store
 *  @param buf Where to store it (may be unaligned)
 */
void u16_to_le(ushort val, ubyte* buf)
{
    static if(HTS_LITTLE_ENDIAN && HTS_ALLOW_UNALIGNED)
        *(cast(uint16_u *) buf) = val;
    else{
        buf[0] = val & 0xff;
        buf[1] = (val >> 8) & 0xff;
    }
}

/// Store a uint value in little-endian byte order
/** @param val The value to store
 *  @param buf Where to store it (may be unaligned)
 */
void u32_to_le(uint val, ubyte* buf)
{
    static if(HTS_LITTLE_ENDIAN && HTS_ALLOW_UNALIGNED)
        *(cast(uint32_u *) buf) = val;
    else{
        buf[0] = val & 0xff;
        buf[1] = (val >> 8) & 0xff;
        buf[2] = (val >> 16) & 0xff;
        buf[3] = (val >> 24) & 0xff;
    }
}

/// Store a ulong value in little-endian byte order
/** @param val The value to store
 *  @param buf Where to store it (may be unaligned)
 */
void u64_to_le(ulong val, ubyte* buf)
{
    static if(HTS_LITTLE_ENDIAN && HTS_ALLOW_UNALIGNED)
        *(cast(uint64_u *) buf) = val;
    else{
        buf[0] = val & 0xff;
        buf[1] = (val >> 8) & 0xff;
        buf[2] = (val >> 16) & 0xff;
        buf[3] = (val >> 24) & 0xff;
        buf[4] = (val >> 32) & 0xff;
        buf[5] = (val >> 40) & 0xff;
        buf[6] = (val >> 48) & 0xff;
        buf[7] = (val >> 56) & 0xff;
    }
}

/* Signed values.  Grab the data as unsigned, then convert to signed without
 * triggering undefined behaviour.  On any sensible platform, the conversion
 * should optimise away to nothing.
 */

/// Get an int8_t value from an unsigned byte array
/** @param buf Pointer to source byte array, may be unaligned
 *  @return A 8 bit signed integer
 *  The input data is interpreted as 2's complement representation.
 */
byte le_to_i8(const(ubyte)* buf)
{
    return *buf < 0x80 ? cast(int8_t) *buf : -((int8_t(cast(byte)0xff) - cast(int8_t)*buf)) - int8_t(1);
}

/// Get an short value from an unsigned byte array
/** @param buf Pointer to source byte array, may be unaligned
 *  @return A 16 bit signed integer
 *  The input data is interpreted as 2's complement representation in
 *  little-endian byte order.
 */
short le_to_i16(const(ubyte)* buf)
{
    ushort v = le_to_u16(buf);
    return v < 0x8000 ? cast(int16_t) v : -((int16_t(cast(short)0xffff) - v)) - cast(int16_t)1;
}

/// Get an int value from an unsigned byte array
/** @param buf Pointer to source byte array, may be unaligned
 *  @return A 32 bit signed integer
 *  The input data is interpreted as 2's complement representation in
 *  little-endian byte order.
 */
int le_to_i32(const(ubyte)* buf)
{
    uint v = le_to_u32(buf);
    return v < 0x80000000U ? cast(int) v : -(cast(int) (0xffffffffU - v)) - 1;
}

/// Get an long value from an unsigned byte array
/** @param buf Pointer to source byte array, may be unaligned
 *  @return A 64 bit signed integer
 *  The input data is interpreted as 2's complement representation in
 *  little-endian byte order.
 */
long le_to_i64(const(ubyte)* buf)
{
    ulong v = le_to_u64(buf);
    return (v < 0x8000000000000000UL
            ? cast(long) v : -(cast(long) (0xffffffffffffffffUL - v)) - 1);
}

// Converting the other way is easier as signed -> unsigned is well defined.

/// Store a ushort value in little-endian byte order
/** @param val The value to store
 *  @param buf Where to store it (may be unaligned)
 */
void i16_to_le(short val, ubyte* buf)
{
    u16_to_le(val, buf);
}

/// Store a uint value in little-endian byte order
/** @param val The value to store
 *  @param buf Where to store it (may be unaligned)
 */
void i32_to_le(int val, ubyte* buf)
{
    u32_to_le(val, buf);
}

/// Store a ulong value in little-endian byte order
/** @param val The value to store
 *  @param buf Where to store it (may be unaligned)
 */
void i64_to_le(long val, ubyte* buf)
{
    u64_to_le(val, buf);
}

/* Floating point.  Assumptions:
 *  Platform uses IEEE 754 format
 *  sizeof(float) == sizeof(uint)
 *  sizeof(double) == sizeof(ulong)
 *  Endian-ness is the same for both floating point and integer
 *  Type-punning via a union is allowed
 */

/// Get a float value from an unsigned byte array
/** @param buf Pointer to source byte array, may be unaligned
 *  @return A 32 bit floating point value
 *  The input is interpreted as an IEEE 754 format float in little-endian
 *  byte order.
 */
float le_to_float(const(ubyte)* buf)
{
    union CONVERT 
    {
        uint u;
        float   f;
    }

    CONVERT convert;
    convert.u = le_to_u32(buf);
    return convert.f;
}

/// Get a double value from an unsigned byte array
/** @param buf Pointer to source byte array, may be unaligned
 *  @return A 64 bit floating point value
 *  The input is interpreted as an IEEE 754 format double in little-endian
 *  byte order.
 */
double le_to_double(const(ubyte)* buf)
{
    union CONVERT 
    {
        ulong u;
        double   f;
    }
    CONVERT convert;
    convert.u = le_to_u64(buf);
    return convert.f;
}

/// Store a float value in little-endian byte order
/** @param val The value to store
 *  @param buf Where to store it (may be unaligned)
 */
void float_to_le(float val, ubyte* buf)
{
    union CONVERT 
    {
        uint u;
        float f;
    }
    CONVERT convert;
    convert.f = val;
    u32_to_le(convert.u, buf);
}

/// Store a double value in little-endian byte order
/** @param val The value to store
 *  @param buf Where to store it (may be unaligned)
 */
void double_to_le(double val, ubyte* buf)
{
    union CONVERT 
    {
        ulong u;
        double f;
    }
    CONVERT convert;
    convert.f = val;
    u64_to_le(convert.u, buf);
}

/* HTS_ENDIAN_H */
