/*  khash_str2int.h -- C-string to integer hash table.

    Copyright (C) 2013-2014,2020 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

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
module htslib.khash_str2int;

extern (C):
alias khint32_t = uint;

alias khint64_t = ulong;

alias khint_t = uint;
alias khiter_t = uint;
alias kh_cstr_t = const(char) *;

struct kh_str2int_s
{
    khint_t n_buckets;
    khint_t size;
    khint_t n_occupied;
    khint_t upper_bound;
    khint32_t* flags;
    kh_cstr_t* keys;
    int* vals;
}

alias kh_str2int_t = kh_str2int_s;
kh_str2int_t* kh_init_str2int ();
void kh_destroy_str2int (kh_str2int_t* h);
void kh_clear_str2int (kh_str2int_t* h);
khint_t kh_get_str2int (const(kh_str2int_t)* h, kh_cstr_t key);
int kh_resize_str2int (kh_str2int_t* h, khint_t new_n_buckets);
khint_t kh_put_str2int (kh_str2int_t* h, kh_cstr_t key, int* ret);
void kh_del_str2int (kh_str2int_t* h, khint_t x);

/*
 *  Wrappers for khash dictionaries used by mpileup.
 */

void* khash_str2int_init ();

/*
 *  Destroy the hash structure, but not the keys
 */

// Note that strings are not freed.
void khash_str2int_destroy (void* _hash);

/*
 *  Destroys both the hash structure and the keys
 */
void khash_str2int_destroy_free (void* _hash);

/*
 *  Returns 1 if key exists or 0 if not
 */
int khash_str2int_has_key (void* _hash, const(char)* str);

/*
 *  Returns 0 on success and -1 when the key is not present. On success,
 *  *value is set, unless NULL is passed.
 */
int khash_str2int_get (void* _hash, const(char)* str, int* value);

/*
 *  Add a new string to the dictionary, auto-incrementing the value.
 *  On success returns the newly inserted integer id, on error -1
 *  is returned. Note that the key must continue to exist throughout
 *  the whole life of _hash.
 */
int khash_str2int_inc (void* _hash, const(char)* str);

/*
 *  Set a new key,value pair. On success returns the bin index, on
 *  error -1 is returned. Note that the key must continue to exist
 *  throughout the whole life of _hash.
 */
int khash_str2int_set (void* _hash, const(char)* str, int value);

/*
 *  Return the number of keys in the hash table.
 */
int khash_str2int_size (void* _hash);

