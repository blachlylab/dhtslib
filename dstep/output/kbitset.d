/* The MIT License

   Copyright (C) 2015, 2018 Genome Research Ltd.

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
import core.stdc.limits;

extern (C):

/* Example of using kbitset_t, which represents a subset of {0,..., N-1},
   where N is the size specified in kbs_init().

	kbitset_t *bset = kbs_init(100);
	kbs_insert(bset, 5);
	kbs_insert(bset, 68);
	kbs_delete(bset, 37);
	// ...

	if (kbs_exists(bset, 68)) printf("68 present\n");

	kbitset_iter_t itr;
	int i;
	kbs_start(&itr);
	while ((i = kbs_next(bset, &itr)) >= 0)
		printf("%d present\n", i);

	kbs_destroy(bset);

   Example of declaring a kbitset_t-using function in a header file, so that
   only source files that actually use process() need to include <kbitset.h>:

	struct kbitset_t;
	void process(struct kbitset_t *bset);
*/

enum KBS_ELTBITS = CHAR_BIT * c_ulong.sizeof;

extern (D) auto KBS_ELT(T)(auto ref T i)
{
    return i / KBS_ELTBITS;
}

extern (D) auto KBS_MASK(T)(auto ref T i)
{
    return 1UL << (i % KBS_ELTBITS);
}

struct kbitset_t
{
    size_t n;
    size_t n_max;
    c_ulong[1] b;
}

// (For internal use only.) Returns a mask (like 00011111) showing
// which bits are in use in the last slot (for the given ni) set.
c_ulong kbs_last_mask (size_t ni);

// Initialise a bit set capable of holding ni integers, 0 <= i < ni.
// The set returned is empty if fill == 0, or all of [0,ni) otherwise.

// b[n] is always non-zero (a fact used by kbs_next()).
kbitset_t* kbs_init2 (size_t ni, int fill);

// Initialise an empty bit set capable of holding ni integers, 0 <= i < ni.
kbitset_t* kbs_init (size_t ni);

// Resize an existing bit set to be capable of holding ni_new integers.
// Elements in [ni_old,ni_new) are added to the set if fill != 0.

// Need to clear excess bits when fill!=0 or n_new<n; always is simpler.
int kbs_resize2 (kbitset_t** bsp, size_t ni_new, int fill);

// Resize an existing bit set to be capable of holding ni_new integers.
// Returns negative on error.
int kbs_resize (kbitset_t** bsp, size_t ni_new);

// Destroy a bit set.
void kbs_destroy (kbitset_t* bs);

// Reset the bit set to empty.
void kbs_clear (kbitset_t* bs);

// Reset the bit set to all of [0,ni).
void kbs_insert_all (kbitset_t* bs);

// Insert an element into the bit set.
void kbs_insert (kbitset_t* bs, int i);

// Remove an element from the bit set.
void kbs_delete (kbitset_t* bs, int i);

// Test whether the bit set contains the element.
int kbs_exists (const(kbitset_t)* bs, int i);

struct kbitset_iter_t
{
    c_ulong mask;
    size_t elt;
    int i;
}

// Initialise or reset a bit set iterator.
void kbs_start (kbitset_iter_t* itr);

// Return the next element contained in the bit set, or -1 if there are no more.
int kbs_next (const(kbitset_t)* bs, kbitset_iter_t* itr);

