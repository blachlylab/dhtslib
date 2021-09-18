module dhtslib.memory;

import std.meta : AliasSeq, staticIndexOf;
import std.traits : isPointer, isSomeFunction, ReturnType;
import core.stdc.stdlib;
import htslib;

struct HtslibMemory(T, alias destroy)
if(!isPointer!T && isSomeFunction!destroy)
{
    /// data pointer
    T * ptr;

    /// reference count pointer
    private size_t * refct;

    /// ability to use this as the ptr directly
    alias ptr this;

    /// ptr ctor
    this(T * ptr)
    {
        this.ptr = ptr;
        this.refct = cast(size_t*)malloc(size_t.sizeof);
        *this.refct = 1;
    }

    /// struct ctor
    this(T data)
    {
        *this.ptr = data;
        this.refct = cast(size_t*)malloc(size_t.sizeof);
        *this.refct = 1;
    }

    /// postblit inc refct
    this(this) @safe pure nothrow @nogc
    {
        if (refct is null) return;
        ++(*this.refct);
    }

    /// dtor dec refct
    /// or destroy
    ~this()
    {
        if (this.refct is null) return;
        assert(*this.refct > 0);
        if (--(*this.refct))
            return;

        /// if destroy function return is void 
        /// just destroy
        /// else if int
        /// destroy then check return value 
        /// else don't compile
        static if(is(ReturnType!destroy == void))
            destroy(this.ptr);
        else static if(is(ReturnType!destroy == int))
        {
            auto success = destroy(this.ptr);
            if(!success) hts_log_error(__FUNCTION__,"Couldn't destroy "~T.stringof~" data using function "~destroy.stringof);
        }else{
            static assert(0, "HtslibMemory doesn't recognize destroy function return type");
        }
        
        /// free refct
        free(this.refct);
    }
}

alias Bam1tPtr = HtslibMemory!(bam1_t, bam_destroy1);

unittest
{
    // A pair of an `int` and a `size_t` - the latter being the
    // reference count - will be dynamically allocated
    auto rc1 = Bam1tPtr(bam_init1);
    assert(rc1.core.pos == 0);
    // No more allocation, add just one extra reference count
    auto rc2 = rc1;
    // Reference semantics
    rc2.core.pos = 42;
    assert(rc1.core.pos == 42);
    // the pair will be freed when rc1 and rc2 go out of scope
}