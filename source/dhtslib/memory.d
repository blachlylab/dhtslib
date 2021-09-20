module dhtslib.memory;

import std.traits : isPointer, isSomeFunction, ReturnType;
import core.lifetime : move;
import std.typecons : RefCounted, RefCountedAutoInitialize;
import htslib.sam : bam1_t, bam_hdr_t, bam_destroy1, bam_hdr_destroy;
import htslib.vcf : bcf1_t,bcf_hdr_t, bcf_destroy, bcf_hdr_destroy;
import htslib.hts : htsFile, hts_close;
import htslib.bgzf : BGZF, bgzf_close;
import htslib.tbx : tbx_t, tbx_destroy;
import htslib.faidx : faidx_t, fai_destroy;
import htslib.hts_log;

/// Template struct that performs reference
/// counting on htslib pointers and destroys with specified function
struct HtslibMemory(T, alias destroy)
if(!isPointer!T && isSomeFunction!destroy)
{
    @safe:
    /// Pointer Wrapper
    struct HtslibPtr
    {
        /// data pointer
        T * ptr;

        /// no copying this as that could result
        /// in premature destruction
        @disable this(this);

        /// destroy 
        ~this() @trusted
        {
            /// if destroy function return is void 
            /// just destroy
            /// else if int
            /// destroy then check return value 
            /// else don't compile
            if(this.ptr){
                static if(is(ReturnType!destroy == void))
                    destroy(this.ptr);
                else static if(is(ReturnType!destroy == int))
                {
                    auto err = destroy(this.ptr);
                    if(err != 0) 
                        hts_log_error(__FUNCTION__,"Couldn't destroy/close "~T.stringof~" * data using function "~__traits(identifier, destroy));
                }else{
                    static assert(0, "HtslibMemory doesn't recognize destroy function return type");
                }
            }
        }
    }

    /// reference counted HtslibPtr
    RefCounted!(HtslibPtr, RefCountedAutoInitialize.yes) rcPtr;

    /// get underlying data pointer
    @property nothrow pure @nogc
    ref inout(T*) getPtr() inout return
    {
        return rcPtr.refCountedPayload.ptr;
    }

    /// allow HtslibMemory to be used as 
    /// underlying ptr type
    alias getPtr this;

    /// ctor from raw pointer
    this(T * rawPtr) @trusted
    {
        auto wrapped = HtslibPtr(rawPtr);
        move(wrapped,this.rcPtr.refCountedPayload);
    }
}

/// reference counted bam1_t wrapper
/// can be used directly as a bam1_t *
alias Bam1_t = HtslibMemory!(bam1_t, bam_destroy1);

/// reference counted bam_hdr_t wrapper
/// can be used directly as a bam_hdr_t *
alias Bam_hdr_t = HtslibMemory!(bam_hdr_t, bam_hdr_destroy);

/// reference counted bcf1_t wrapper
/// can be used directly as a bcf1_t *
alias Bcf1_t = HtslibMemory!(bcf1_t, bcf_destroy);

/// reference counted bcf_hdr_t wrapper
/// can be used directly as a bcf_hdr_t *
alias Bcf_hdr_t = HtslibMemory!(bcf_hdr_t, bcf_hdr_destroy);

/// reference counted htsFile wrapper
/// can be used directly as a htsFile *
alias HtsFile = HtslibMemory!(htsFile, hts_close);

/// reference counted htsFile wrapper
/// can be used directly as a htsFile *
alias VcfFile = HtsFile;

/// reference counted tbx_t wrapper
/// can be used directly as a tbx_t *
alias Tbx_t = HtslibMemory!(tbx_t, tbx_destroy);

/// reference counted BGZF wrapper
/// can be used directly as a BGZF *
alias BgzfPtr = HtslibMemory!(BGZF, bgzf_close);

/// reference counted faidx_t wrapper
/// can be used directly as a faidx_t *
alias Faidx_t = HtslibMemory!(faidx_t, fai_destroy);


unittest
{
    import htslib.sam : bam_init1;
    auto rc1 = Bam1_t(bam_init1);
    assert(rc1.core.pos == 0);
    // No more allocation, add just one extra reference count
    auto rc2 = rc1;
    // Reference semantics
    rc2.core.pos = 42;
    assert(rc1.core.pos == 42);
}