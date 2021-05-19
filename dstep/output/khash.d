/* The MIT License

   Copyright (c) 2008, 2009, 2011 by Attractive Chaos <attractor@live.co.uk>
   Copyright (C) 2014-2015, 2018 Genome Research Ltd.

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

/*
  An example:

#include "khash.h"
KHASH_MAP_INIT_INT(32, char)
int main() {
	int ret, is_missing;
	khiter_t k;
	khash_t(32) *h = kh_init(32);
	k = kh_put(32, h, 5, &ret);
	kh_value(h, k) = 10;
	k = kh_get(32, h, 10);
	is_missing = (k == kh_end(h));
	k = kh_get(32, h, 5);
	kh_del(32, h, k);
	for (k = kh_begin(h); k != kh_end(h); ++k)
		if (kh_exist(h, k)) kh_value(h, k) = 1;
	kh_destroy(32, h);
	return 0;
}
*/

/*
  2013-05-02 (0.2.8):

	* Use quadratic probing. When the capacity is power of 2, stepping function
	  i*(i+1)/2 guarantees to traverse each bucket. It is better than double
	  hashing on cache performance and is more robust than linear probing.

	  In theory, double hashing should be more robust than quadratic probing.
	  However, my implementation is probably not for large hash tables, because
	  the second hash function is closely tied to the first hash function,
	  which reduce the effectiveness of double hashing.

	Reference: http://research.cs.vt.edu/AVresearch/hashing/quadratic.php

  2011-12-29 (0.2.7):

    * Minor code clean up; no actual effect.

  2011-09-16 (0.2.6):

	* The capacity is a power of 2. This seems to dramatically improve the
	  speed for simple keys. Thank Zilong Tan for the suggestion. Reference:

	   - http://code.google.com/p/ulib/
	   - http://nothings.org/computer/judy/

	* Allow to optionally use linear probing which usually has better
	  performance for random input. Double hashing is still the default as it
	  is more robust to certain non-random input.

	* Added Wang's integer hash function (not used by default). This hash
	  function is more robust to certain non-random input.

  2011-02-14 (0.2.5):

    * Allow to declare global functions.

  2009-09-26 (0.2.4):

    * Improve portability

  2008-09-19 (0.2.3):

	* Corrected the example
	* Improved interfaces

  2008-09-11 (0.2.2):

	* Improved speed a little in kh_put()

  2008-09-10 (0.2.1):

	* Added kh_clear()
	* Fixed a compiling error

  2008-09-02 (0.2.0):

	* Changed to token concatenation which increases flexibility.

  2008-08-31 (0.1.2):

	* Fixed a bug in kh_get(), which has not been tested previously.

  2008-08-31 (0.1.1):

	* Added destructor
*/

import core.stdc.config;
import core.stdc.string;

extern (C):

/*!
  @header

  Generic hash table library.
 */

enum AC_VERSION_KHASH_H = "0.2.8";

/* compiler specific configuration */

alias khint32_t = uint;

alias khint64_t = c_ulong;

/* kh_inline */

/* klib_unused */

alias khint_t = uint;
alias khiter_t = uint;

extern (D) auto __ac_isempty(T0, T1)(auto ref T0 flag, auto ref T1 i)
{
    return (flag[i >> 4] >> ((i & 0xfU) << 1)) & 2;
}

extern (D) auto __ac_isdel(T0, T1)(auto ref T0 flag, auto ref T1 i)
{
    return (flag[i >> 4] >> ((i & 0xfU) << 1)) & 1;
}

extern (D) auto __ac_iseither(T0, T1)(auto ref T0 flag, auto ref T1 i)
{
    return (flag[i >> 4] >> ((i & 0xfU) << 1)) & 3;
}

extern (D) auto __ac_fsize(T)(auto ref T m)
{
    return m < 16 ? 1 : m >> 4;
}

alias kcalloc = calloc;

alias kmalloc = malloc;

alias krealloc = realloc;

alias kfree = free;

extern __gshared const double __ac_HASH_UPPER;

/* This function uses 0.25*n_buckets bytes of working space instead of [sizeof(key_t+val_t)+.25]*n_buckets. */

/* requested size is too small */
/* hash table size to be changed (shrink or expand); rehash */

/* expand */

/* otherwise shrink */

/* rehashing is needed */

/* kick-out process; sort of like in Cuckoo hashing */

/* kick out the existing element */

/* mark it as deleted in the old hash table */
/* write the element and jump out of the loop */

/* shrink the hash table */

/* free the working space */

/* update the hash table */

/* clear "deleted" elements */

/* expand the hash table */

/* TODO: to implement automatically shrinking; resize() already support shrinking */

/* for speed up */

/* not present at all */

/* deleted */

/* Don't touch h->keys[x] if present and not deleted */

/* --- BEGIN OF HASH FUNCTIONS --- */

/*! @function
  @abstract     Integer hash function
  @param  key   The integer [khint32_t]
  @return       The hash value [khint_t]
 */
extern (D) auto kh_int_hash_func(T)(auto ref T key)
{
    return cast(khint32_t) key;
}

/*! @function
  @abstract     Integer comparison function
 */
extern (D) auto kh_int_hash_equal(T0, T1)(auto ref T0 a, auto ref T1 b)
{
    return a == b;
}

/*! @function
  @abstract     64-bit integer hash function
  @param  key   The integer [khint64_t]
  @return       The hash value [khint_t]
 */
extern (D) auto kh_int64_hash_func(T)(auto ref T key)
{
    return cast(khint32_t) key >> 33 ^ key ^ key << 11;
}

/*! @function
  @abstract     64-bit integer comparison function
 */
extern (D) auto kh_int64_hash_equal(T0, T1)(auto ref T0 a, auto ref T1 b)
{
    return a == b;
}

/*! @function
  @abstract     const char* hash function
  @param  s     Pointer to a null terminated string
  @return       The hash value
 */
khint_t __ac_X31_hash_string (const(char)* s);
/*! @function
  @abstract     Another interface to const char* hash function
  @param  key   Pointer to a nul terminated string [const char*]
  @return       The hash value [khint_t]
 */
alias kh_str_hash_func = __ac_X31_hash_string;
/*! @function
  @abstract     Const char* comparison function
 */
extern (D) auto kh_str_hash_equal(T0, T1)(auto ref T0 a, auto ref T1 b)
{
    return strcmp(a, b) == 0;
}

/*! @function
  @abstract     Kstring hash function
  @param  s     Pointer to a kstring
  @return       The hash value
 */
khint_t __ac_X31_hash_kstring (const kstring_t ks);
/*! @function
  @abstract     Interface to kstring hash function.
  @param  key   Pointer to a khash; permits hashing on non-nul terminated strings.
  @return       The hash value [khint_t]
 */
alias kh_kstr_hash_func = __ac_X31_hash_kstring;
/*! @function
  @abstract     kstring comparison function
 */
extern (D) auto kh_kstr_hash_equal(T0, T1)(auto ref T0 a, auto ref T1 b)
{
    return a.l == b.l && strncmp(a.s, b.s, a.l) == 0;
}

khint_t __ac_Wang_hash (khint_t key);

extern (D) auto kh_int_hash_func2(T)(auto ref T k)
{
    return __ac_Wang_hash(cast(khint_t) key);
}

/* --- END OF HASH FUNCTIONS --- */

/* Other convenient macros... */

/*!
  @abstract Type of the hash table.
  @param  name  Name of the hash table [symbol]
 */
extern (D) string khash_t(T)(auto ref T name)
{
    import std.conv : to;

    return "kh_" ~ to!string(name) ~ "_t";
}

/*! @function
  @abstract     Initiate a hash table.
  @param  name  Name of the hash table [symbol]
  @return       Pointer to the hash table [khash_t(name)*]
 */

/*! @function
  @abstract     Destroy a hash table.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [khash_t(name)*]
 */

/*! @function
  @abstract     Reset a hash table without deallocating memory.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [khash_t(name)*]
 */

/*! @function
  @abstract     Resize a hash table.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  s     New size [khint_t]
 */

/*! @function
  @abstract     Insert a key to the hash table.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  k     Key [type of keys]
  @param  r     Extra return code: -1 if the operation failed;
                0 if the key is present in the hash table;
                1 if the bucket is empty (never used); 2 if the element in
				the bucket has been deleted [int*]
  @return       Iterator to the inserted element [khint_t]
 */

/*! @function
  @abstract     Retrieve a key from the hash table.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  k     Key [type of keys]
  @return       Iterator to the found element, or kh_end(h) if the element is absent [khint_t]
 */

/*! @function
  @abstract     Remove a key from the hash table.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  k     Iterator to the element to be deleted [khint_t]
 */

/*! @function
  @abstract     Test whether a bucket contains data.
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  x     Iterator to the bucket [khint_t]
  @return       1 if containing data; 0 otherwise [int]
 */
extern (D) auto kh_exist(T0, T1)(auto ref T0 h, auto ref T1 x)
{
    return !__ac_iseither(h.flags, x);
}

/*! @function
  @abstract     Get key given an iterator
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  x     Iterator to the bucket [khint_t]
  @return       Key [type of keys]
 */
extern (D) auto kh_key(T0, T1)(auto ref T0 h, auto ref T1 x)
{
    return h.keys[x];
}

/*! @function
  @abstract     Get value given an iterator
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  x     Iterator to the bucket [khint_t]
  @return       Value [type of values]
  @discussion   For hash sets, calling this results in segfault.
 */
extern (D) auto kh_val(T0, T1)(auto ref T0 h, auto ref T1 x)
{
    return h.vals[x];
}

/*! @function
  @abstract     Alias of kh_val()
 */
extern (D) auto kh_value(T0, T1)(auto ref T0 h, auto ref T1 x)
{
    return h.vals[x];
}

/*! @function
  @abstract     Get the start iterator
  @param  h     Pointer to the hash table [khash_t(name)*]
  @return       The start iterator [khint_t]
 */
extern (D) auto kh_begin(T)(auto ref T h)
{
    return cast(khint_t) 0;
}

/*! @function
  @abstract     Get the end iterator
  @param  h     Pointer to the hash table [khash_t(name)*]
  @return       The end iterator [khint_t]
 */
extern (D) auto kh_end(T)(auto ref T h)
{
    return h.n_buckets;
}

/*! @function
  @abstract     Get the number of elements in the hash table
  @param  h     Pointer to the hash table [khash_t(name)*]
  @return       Number of elements in the hash table [khint_t]
 */
extern (D) auto kh_size(T)(auto ref T h)
{
    return h.size;
}

/*! @function
  @abstract     Get the number of buckets in the hash table
  @param  h     Pointer to the hash table [khash_t(name)*]
  @return       Number of buckets in the hash table [khint_t]
 */
extern (D) auto kh_n_buckets(T)(auto ref T h)
{
    return h.n_buckets;
}

/*! @function
  @abstract     Iterate over the entries in the hash table
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  kvar  Variable to which key will be assigned
  @param  vvar  Variable to which value will be assigned
  @param  code  Block of code to execute
 */

/*! @function
  @abstract     Iterate over the values in the hash table
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  vvar  Variable to which value will be assigned
  @param  code  Block of code to execute
 */

/* More convenient interfaces */

/*! @function
  @abstract     Instantiate a hash set containing integer keys
  @param  name  Name of the hash table [symbol]
 */
extern (D) auto KHASH_SET_INIT_INT(T)(auto ref T name)
{
    return KHASH_INIT(name, khint32_t, char, 0, kh_int_hash_func, kh_int_hash_equal);
}

/*! @function
  @abstract     Instantiate a hash map containing integer keys
  @param  name  Name of the hash table [symbol]
  @param  khval_t  Type of values [type]
 */
extern (D) auto KHASH_MAP_INIT_INT(T0, T1)(auto ref T0 name, auto ref T1 khval_t)
{
    return KHASH_INIT(name, khint32_t, khval_t, 1, kh_int_hash_func, kh_int_hash_equal);
}

/*! @function
  @abstract     Instantiate a hash set containing 64-bit integer keys
  @param  name  Name of the hash table [symbol]
 */
extern (D) auto KHASH_SET_INIT_INT64(T)(auto ref T name)
{
    return KHASH_INIT(name, khint64_t, char, 0, kh_int64_hash_func, kh_int64_hash_equal);
}

/*! @function
  @abstract     Instantiate a hash map containing 64-bit integer keys
  @param  name  Name of the hash table [symbol]
  @param  khval_t  Type of values [type]
 */
extern (D) auto KHASH_MAP_INIT_INT64(T0, T1)(auto ref T0 name, auto ref T1 khval_t)
{
    return KHASH_INIT(name, khint64_t, khval_t, 1, kh_int64_hash_func, kh_int64_hash_equal);
}

alias kh_cstr_t = const(char)*;
/*! @function
  @abstract     Instantiate a hash set containing const char* keys
  @param  name  Name of the hash table [symbol]
 */
extern (D) auto KHASH_SET_INIT_STR(T)(auto ref T name)
{
    return KHASH_INIT(name, kh_cstr_t, char, 0, kh_str_hash_func, kh_str_hash_equal);
}

/*! @function
  @abstract     Instantiate a hash map containing const char* keys
  @param  name  Name of the hash table [symbol]
  @param  khval_t  Type of values [type]
 */
extern (D) auto KHASH_MAP_INIT_STR(T0, T1)(auto ref T0 name, auto ref T1 khval_t)
{
    return KHASH_INIT(name, kh_cstr_t, khval_t, 1, kh_str_hash_func, kh_str_hash_equal);
}

/*! @function
  @abstract     Instantiate a hash set containing kstring_t keys
  @param  name  Name of the hash table [symbol]
 */
extern (D) auto KHASH_SET_INIT_KSTR(T)(auto ref T name)
{
    return KHASH_INIT(name, kstring_t, char, 0, kh_kstr_hash_func, kh_kstr_hash_equal);
}

/*! @function
  @abstract     Instantiate a hash map containing kstring_t keys
  @param  name  Name of the hash table [symbol]
  @param  khval_t  Type of values [type]
 */
extern (D) auto KHASH_MAP_INIT_KSTR(T0, T1)(auto ref T0 name, auto ref T1 khval_t)
{
    return KHASH_INIT(name, kstring_t, khval_t, 1, kh_kstr_hash_func, kh_kstr_hash_equal);
}

/* __AC_KHASH_H */
