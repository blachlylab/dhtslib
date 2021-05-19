/* The MIT License

   Copyright (c) 2008-2009, by Attractive Chaos <attractor@live.co.uk>
   Copyright (C) 2013, 2015 Genome Research Ltd.

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

extern (C):

/* klib_unused */

extern (D) string kmempool_t(T)(auto ref T name)
{
    import std.conv : to;

    return "kmp_" ~ to!string(name) ~ "_t";
}

extern (D) string kliter_t(T)(auto ref T name)
{
    import std.conv : to;

    return "kl1_" ~ to!string(name);
}

extern (D) string klist_t(T)(auto ref T name)
{
    import std.conv : to;

    return "kl_" ~ to!string(name) ~ "_t";
}

extern (D) auto kl_val(T)(auto ref T iter)
{
    return iter.data;
}

extern (D) auto kl_next(T)(auto ref T iter)
{
    return iter.next;
}

extern (D) auto kl_begin(T)(auto ref T kl)
{
    return kl.head;
}

extern (D) auto kl_end(T)(auto ref T kl)
{
    return kl.tail;
}
