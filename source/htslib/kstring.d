module htslib.kstring;

extern (C):
/* The MIT License

   Copyright (C) 2011 by Attractive Chaos <attractor@live.co.uk>

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

import core.stdc.stdlib : realloc;
import core.stdc.string : memcpy, strlen;
// import core.stdc.stdarg;
import core.stdc.stdint;
import core.stdc.stdio : EOF;
// import core.stdc.limits;

/+
#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

+/
/// round 32 or 64 bit (u)int x to power of 2 that is equal or greater
pragma(inline, true)
void kroundup_size_t(ref size_t x) {
	x -= 1;
	x |= (x >> 1);
	x |= (x >> 2);
	x |= (x >> 4);
	x |= (x >> 8);
	x |= (x >> 16);

	static if (sizeof(size_t) == 8)
        x |= (x >> 32);

	return ++x;
}

/+
#if defined __GNUC__ && (__GNUC__ > 2 || (__GNUC__ == 2 && __GNUC_MINOR__ > 4))
#define KS_ATTR_PRINTF(fmt, arg) __attribute__((__format__ (__printf__, fmt, arg)))
#else
#define KS_ATTR_PRINTF(fmt, arg)
#endif
+/

/** kstring_t is a simple non-opaque type whose fields are likely to be
 * used directly by user code (but see also ks_str() and ks_len() below).
 * A kstring_t object is initialised by either of
 *       kstring_t str = { 0, 0, NULL };
 *       kstring_t str; ...; str.l = str.m = 0; str.s = NULL;
 * (NB: in Dlang it is automagically initialized to zeros unless declared = void)
 * and either ownership of the underlying buffer should be given away before
 * the object disappears (see ks_release() below) or the kstring_t should be
 * destroyed with  free(str.s);  */
struct __kstring_t { // @suppress(dscanner.style.phobos_naming_convention)
	/// length(l) of max(m)
	size_t l, m;
	char *s;	/// pointer to string -- must free(str.s) when done with it
} 
alias kstring_t = __kstring_t;

/// ???
struct ks_tokaux_t { // @suppress(dscanner.style.phobos_naming_convention)
	uint64_t[4] tab;	/// ???
	/// ???
	int sep, finished;
	const(char) *p;		/// end of the current token
}
/+
	int kvsprintf(kstring_t *s, const char *fmt, va_list ap) KS_ATTR_PRINTF(2,0);
	int ksprintf(kstring_t *s, const char *fmt, ...) KS_ATTR_PRINTF(2,3);
    int kputd(double d, kstring_t *s); // custom %g only handler
	int ksplit_core(char *s, int delimiter, int *_max, int **_offsets);
	char *kstrstr(const char *str, const char *pat, int **_prep);
	char *kstrnstr(const char *str, const char *pat, int n, int **_prep);
	void *kmemmem(const void *_str, int n, const void *_pat, int m, int **_prep);

	/* kstrtok() is similar to strtok_r() except that str is not
	 * modified and both str and sep can be NULL. For efficiency, it is
	 * actually recommended to set both to NULL in the subsequent calls
	 * if sep is not changed. */
	char *kstrtok(const char *str, const char *sep, ks_tokaux_t *aux);

	/* kgetline() uses the supplied fgets()-like function to read a "\n"-
	 * or "\r\n"-terminated line from fp.  The line read is appended to the
	 * kstring without its terminator and 0 is returned; EOF is returned at
	 * EOF or on error (determined by querying fp, as per fgets()). */
	typedef char *kgets_func(char *, int, void *);
	int kgetline(kstring_t *s, kgets_func *fgets, void *fp);
+/

pragma(inline, true)
int ks_resize(kstring_t *s, size_t size)
{
	if (s.m < size) {
		char *tmp;
		kroundup_size_t(size);
		tmp = cast(char*)realloc(s.s, size);
		if (!tmp)
			return -1;
		s.s = tmp;
		s.m = size;
	}
	return 0;
}

/+
static inline char *ks_str(kstring_t *s)
{
	return s->s;
}

static inline size_t ks_len(kstring_t *s)
{
	return s->l;
}

// Give ownership of the underlying buffer away to something else (making
// that something else responsible for freeing it), leaving the kstring_t
// empty and ready to be used again, or ready to go out of scope without
// needing  free(str.s)  to prevent a memory leak.
static inline char *ks_release(kstring_t *s)
{
	char *ss = s->s;
	s->l = s->m = 0;
	s->s = NULL;
	return ss;
}
+/

pragma(inline, true)
int kputsn(const(char) *p, size_t l, kstring_t *s)
{
    if (l > SIZE_MAX - 2 - s.l || ks_resize(s, s.l + l + 2) < 0)
		return EOF;
	memcpy(s.s + s.l, p, l);
	s.l += l;
	s.s[s.l] = 0;
	return l;
}

pragma(inline, true)
int kputs(const(char) *p, kstring_t *s)
{
	return kputsn(p, strlen(p), s);
}

pragma(inline, true)
int kputc(int c, kstring_t *s)
{
	if (ks_resize(s, s.l + 2) < 0)
		return EOF;
	s.s[s.l++] = c;
	s.s[s.l] = 0;
	return c;
}

/+
static inline int kputc_(int c, kstring_t *s)
{
	if (ks_resize(s, s->l + 1) < 0)
		return EOF;
	s->s[s->l++] = c;
	return 1;
}

static inline int kputsn_(const void *p, size_t l, kstring_t *s)
{
	if (l > SIZE_MAX - s->l || ks_resize(s, s->l + l) < 0)
		return EOF;
	memcpy(s->s + s->l, p, l);
	s->l += l;
	return l;
}
+/

pragma(inline, true)
int kputw(int c, kstring_t *s)
{
	char[16] buf;
	int i, l = 0;
	uint x = c;
	if (c < 0) x = -x;
	do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
	if (c < 0) buf[l++] = '-';
	if (ks_resize(s, s.l + l + 2) < 0)
		return EOF;
	for (i = l - 1; i >= 0; --i) s.s[s.l++] = buf[i];
	s.s[s.l] = 0;
	return 0;
}

/+
static inline int kputuw(unsigned c, kstring_t *s)
{
	char buf[16];
	int l, i;
	unsigned x;
	if (c == 0) return kputc('0', s);
	for (l = 0, x = c; x > 0; x /= 10) buf[l++] = x%10 + '0';
	if (ks_resize(s, s->l + l + 2) < 0)
		return EOF;
	for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
	s->s[s->l] = 0;
	return 0;
}

static inline int kputl(long c, kstring_t *s)
{
	char buf[32];
	int i, l = 0;
	unsigned long x = c;
	if (c < 0) x = -x;
	do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
	if (c < 0) buf[l++] = '-';
	if (ks_resize(s, s->l + l + 2) < 0)
		return EOF;
	for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
	s->s[s->l] = 0;
	return 0;
}

/*
 * Returns 's' split by delimiter, with *n being the number of components;
 *         NULL on failue.
 */
static inline int *ksplit(kstring_t *s, int delimiter, int *n)
{
	int max = 0, *offsets = 0;
	*n = ksplit_core(s->s, delimiter, &max, &offsets);
	return offsets;
}

+/