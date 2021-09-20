module dhtslib.memory;

import std.traits : isPointer, isSomeFunction, ReturnType;
import core.lifetime : move;
import std.typecons : RefCounted, RefCountedAutoInitialize;
import htslib.sam;

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
            static if(is(ReturnType!destroy == void))
                destroy(this.ptr);
            else static if(is(ReturnType!destroy == int))
            {
                auto success = destroy(this.ptr);
                if(!success) hts_log_error(__FUNCTION__,"Couldn't destroy "~T.stringof~" * data using function "~destroy.stringof);
            }else{
                static assert(0, "HtslibMemory doesn't recognize destroy function return type");
            }
        }
    }

    /// reference counted HtslibPtr
    RefCounted!(HtslibPtr, RefCountedAutoInitialize.yes) rcPtr;

    /// get underlying data pointer
    @property nothrow pure @nogc
    ref inout(T*) truePtr() inout return
    {
        return rcPtr.refCountedPayload.ptr;
    }

    /// allow HtslibMemory to be used as 
    /// underlying ptr type
    alias truePtr this;

    /// ctor from raw pointer
    this(T * rawPtr) @trusted
    {
        auto wrapped = HtslibPtr(rawPtr);
        move(wrapped,this.rcPtr.refCountedPayload);
    }
}

alias Bam1_t = HtslibMemory!(bam1_t, bam_destroy1);

unittest
{
    auto rc1 = Bam1_t(bam_init1);
    assert(rc1.core.pos == 0);
    // No more allocation, add just one extra reference count
    auto rc2 = rc1;
    // Reference semantics
    rc2.core.pos = 42;
    assert(rc1.core.pos == 42);
}