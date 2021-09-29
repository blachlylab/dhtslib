module dhtslib.memory;

import std.traits : isPointer, isSomeFunction, ReturnType, isSafe;
import core.stdc.stdlib : calloc, malloc, free;
import core.lifetime : move;
import std.typecons : RefCounted, RefCountedAutoInitialize;
import core.atomic : atomicOp;
import htslib;

/// can we use @live for scope checking? 
enum dip1000Enabled = isSafe!((int x) => *&x);

static if(dip1000Enabled)
    pragma(msg, "Using -dip1000 for scope checking and safety");

pragma(inline, true):
/**! Logs an event with severity HTS_LOG_ERROR and compile-time context. */
private void hts_log_errorNoGC(const(char)[] ctx)( string msg) @trusted @nogc nothrow
{
    string open_error_color = "\x1b[0;31m";
    string close_color      = "\x1b[0m";
    static newCtx = ctx ~ '\0';
    hts_log(htsLogLevel.HTS_LOG_ERROR, newCtx.ptr, "%.*s%.*s%.*s",
            cast(int)open_error_color.length, open_error_color.ptr,
            cast(int)msg.length, msg.ptr,
            cast(int)close_color.length, close_color.ptr);
}

/// Template struct that wraps an htslib
/// pointer and reference counts it and then
/// destroys with destroyFun when it goes 
/// truly out of scope
struct SafeHtslibPtr(T, alias destroyFun)
if(!isPointer!T && isSomeFunction!destroyFun)
{
    @safe @nogc nothrow:

    /// data pointer
    T * ptr;
    /// reference counting
    int* refct;
    /// initialized?
    bool initialized;

    /// ctor that respects scope
    this(T * rawPtr) @trusted return scope
    {
        this.ptr = rawPtr;
        this.refct = cast(int *) calloc(int.sizeof,1);
        (*this.refct) = 1;
        this.initialized = true;
    }
    
    /// postblit that respects scope
    this(this) @trusted return scope
    {
        if(initialized)(*this.refct)++;
    }

    /// allow SafeHtslibPtr to be used as 
    /// underlying ptr type
    alias getRef this;

    /// get underlying data pointer
    @property nothrow pure @nogc
    ref inout(T*) getRef() inout return
    {
        return ptr;
    }

    /// take ownership of underlying data pointer
    @property nothrow pure @nogc
    T* moveRef()
    {
        T * ptr;
        move(this.getRef, ptr);
        return ptr;
    }

    /// dtor that respects scope
    ~this() @trusted return scope
    {
        
        if(!this.initialized) return;
        if(--(*this.refct)) return;
        if(this.ptr){
            free(this.refct);
            /// if destroy function return is void 
            /// just destroy
            /// else if return is int
            /// destroy then check return value 
            /// else don't compile
            static if(is(ReturnType!destroyFun == void))
                destroyFun(this.ptr);
            else static if(is(ReturnType!destroyFun == int))
            {
                auto err = destroyFun(this.ptr);
                if(err != 0) 
                    hts_log_errorNoGC!__FUNCTION__("Couldn't destroy/close "~T.stringof~" * data using function "~__traits(identifier, destroyFun));
            }else{
                static assert(0, "HtslibPtr doesn't recognize destroy function return type");
            }
        }
    }
}

/// reference counted bam1_t wrapper
/// can be used directly as a bam1_t *
alias Bam1 = SafeHtslibPtr!(bam1_t, bam_destroy1);

/// reference counted bam_hdr_t wrapper
/// can be used directly as a bam_hdr_t *
alias BamHdr = SafeHtslibPtr!(bam_hdr_t, bam_hdr_destroy);

/// reference counted bcf1_t wrapper
/// can be used directly as a bcf1_t *
alias Bcf1 = SafeHtslibPtr!(bcf1_t, bcf_destroy);

/// reference counted bcf_hdr_t wrapper
/// can be used directly as a bcf_hdr_t *
alias BcfHdr = SafeHtslibPtr!(bcf_hdr_t, bcf_hdr_destroy);

/// reference counted htsFile wrapper
/// can be used directly as a htsFile *
alias HtsFile = SafeHtslibPtr!(htsFile, hts_close);

/// reference counted htsFile wrapper
/// can be used directly as a htsFile *
alias HtsIdx = SafeHtslibPtr!(hts_idx_t, hts_idx_destroy);

/// reference counted htsFile wrapper
/// can be used directly as a htsFile *
alias HtsItr = SafeHtslibPtr!(hts_itr_t, hts_itr_destroy);

/// reference counted htsFile wrapper
/// can be used directly as a htsFile *
alias VcfFile = HtsFile;

/// reference counted tbx_t wrapper
/// can be used directly as a tbx_t *
alias Tbx = SafeHtslibPtr!(tbx_t, tbx_destroy);

/// reference counted BGZF wrapper
/// can be used directly as a BGZF *
alias Bgzf = SafeHtslibPtr!(BGZF, bgzf_close);

/// reference counted faidx_t wrapper
/// can be used directly as a faidx_t *
alias Faidx = SafeHtslibPtr!(faidx_t, fai_destroy);

/// reference counted Kstring wrapper
/// can be used directly as a kstring_t *
alias Kstring = SafeHtslibPtr!(kstring_t, ks_free);

alias HtsItrMulti = HtsItr;

debug(dhtslib_unittest) unittest 
{
    auto rc1 = Bam1(bam_init1);
    assert(rc1.core.pos == 0);
    // No more allocation, add just one extra reference count
    auto rc2 = rc1;
    // Reference semantics
    rc2.core.pos = 42;
    assert(rc1.core.pos == 42);
}

static if(dip1000Enabled){
    debug(dhtslib_unittest) unittest
    {
        auto rc1 = Bam1(bam_init1);
        assert(rc1.core.pos == 0);
        // No more allocation, add just one extra reference count
        auto rc2 = rc1;
        // Reference semantics
        rc2.core.pos = 42;
        assert(rc1.core.pos == 42);
    }

    debug(dhtslib_unittest) unittest 
    {
        auto testfun(bool noScope = false)()
        {
            // auto rc = Bam1(bam_init);
            auto rc1 = Bam1(bam_init1);
            assert(rc1.core.pos == 0);
            // No more allocation, add just one extra reference count
            static if(noScope)
                auto rc2 = rc1;
            else
                scope rc2 = rc1;
            // Reference semantics
            rc2.core.pos = 42;
            assert(rc1.core.pos == 42);
            return rc2.getRef;
        }
        static assert(!__traits(compiles,testfun!true));
        static assert(!__traits(compiles,testfun!false));
    }

    debug(dhtslib_unittest) unittest 
    {
        auto testfun(bool noScope = false)()
        {
            struct T
            {
                Bam1 rc1;
            }
            T test;
            test.rc1 = Bam1(bam_init1);

            assert(test.rc1.core.pos == 0);
            // No more allocation, add just one extra reference count
            static if(noScope)
                auto rc2 = test.rc1;
            else
                scope rc2 = rc1;
            // Reference semantics
            rc2.core.pos = 42;
            assert(test.rc1.core.pos == 42);
            return rc2.getRef;
        }
        static assert(!__traits(compiles,testfun!true));
        static assert(!__traits(compiles,testfun!false));
    }
}