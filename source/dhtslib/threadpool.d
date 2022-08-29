/**
This module provides multithreading support through htslib. It allows for mapping delagates,
and functions via a range interface with `parallelMap`. Crucially, the htslib task pool performs
all work in order. All functions and interfaces below should as well.
*/
module dhtslib.threadpool;
import dhtslib.memory;
import htslib.hts;
import htslib.thread_pool;
import std.algorithm : map;
import std.parallelism : totalCPUs;
import std.traits;
import std.range;
import std.meta;
import core.stdc.stdlib : calloc, free;
import core.atomic : atomicOp;

/// global thread pool
__gshared ThreadPool globalPool;

/// ensure global pool is initialized
/// NOTE: only sets threads if the global pool is not already initialized
void enforceGlobalThreadPool(int threads = -2) @trusted nothrow @nogc 
{
    if(globalPool.tpool && (globalPool.threads == threads || threads == -2)) {
        return;
    } else
        globalPool = ThreadPool(threads, true);   
}

private template noUnsharedAliasing(T)
{
    enum bool noUnsharedAliasing = !hasUnsharedAliasing!T;
}

// This template tests whether a function may be executed in parallel from
// @safe code via Task.executeInNewThread().  There is an additional
// requirement for executing it via a TaskPool.  (See isSafeReturn).
private template isSafeTask(F)
{
    enum bool isSafeTask =
        (functionAttributes!F & (FunctionAttribute.safe | FunctionAttribute.trusted)) != 0 &&
        (functionAttributes!F & FunctionAttribute.ref_) == 0 &&
        (isFunctionPointer!F || !hasUnsharedAliasing!F) &&
        allSatisfy!(noUnsharedAliasing, Parameters!F);
}

/// Thread pool instance. Wraps htsThreadPool and hts_tpool. 
struct ThreadPool 
{
    htsThreadPool htsPool;
    HtsTPool _tpool;
    int threads;
    bool isGlobal;
    bool owned;
    @nogc @trusted nothrow: 
    this(int threads, bool global = false) {
        this.threads = threads;
        if (threads == -1)
        {
            this.threads = totalCPUs;
        }
        this._tpool = HtsTPool(hts_tpool_init(cast(int)this.threads));
        this.htsPool.pool = this._tpool;
        this.htsPool.qsize = 0;
        this.isGlobal = global;
        this.owned = true;
    }

    this(this) {
        if(this.owned) this._tpool = _tpool;
    }

    @property tpool() {
        return this.htsPool.pool;
    }

}


/**
Below code is ported from std.parallelism to work with htslib's thread pool.
It is BSL1.0 licensed.

Source:    $(PHOBOSSRC std/parallelism.d)
Author:  David Simcha
Copyright:  Copyright (c) 2009-2011, David Simcha.
License:    $(HTTP boost.org/LICENSE_1_0.txt, Boost License 1.0)
*/

///
T* addressOf(T)(ref T val) { return &val; }

/**
`Task` represents the fundamental unit of work.  A `Task` may be
executed in parallel with any other `Task`.  Using this struct directly
allows future/promise parallelism.  In this paradigm, a function (or delegate
or other callable) is executed in a thread other than the one it was called
from.  The calling thread does not block while the function is being executed.

The $(REF task, std,parallelism) function can
be used to create an instance of this struct.  See `task` for usage examples.
If `fun` returns by ref, the reference will point
to the returned reference of `fun`.  Otherwise it will point to a
field in this struct.
Copying of this struct is disabled, since it would provide no useful semantics.
If you want to pass this struct around, you should do so by reference or
pointer.
Bugs:  Changes to `ref` and `out` arguments are not propagated to the
       call site, only to `args` in this struct.
*/
struct Task(alias fun, Args...)
{
    ThreadPool * pool;
    private bool isScoped;  // True if created with scopedTask.
    Args _args;

    extern(C) static void * runTask(void* myTask)
    {

        Task* myCastedTask = cast(typeof(this)*) myTask;
        static if (is(ReturnType == void))
        {
            fun(myCastedTask._args);
        }
        else static if (is(typeof(&(fun(myCastedTask._args)))))
        {
            myCastedTask.returnVal = addressOf(fun(myCastedTask._args));
        }
        else
        {
            myCastedTask.returnVal = fun(myCastedTask._args);
        }
        return myTask;
    }

    /**
    The arguments the function was called with.  Changes to `out` and
    `ref` arguments will be visible here.
    */
    static if (__traits(isSame, fun, run))
    {
        alias args = _args[1..$];
    }
    else
    {
        alias args = _args;
    }


    // The purpose of this code is to decide whether functions whose
    // return values have unshared aliasing can be executed via
    // TaskPool from @safe code.  See isSafeReturn.
    static if (__traits(isSame, fun, run))
    {
        static if (isFunctionPointer!(_args[0]))
        {
            private enum bool isPure =
            (functionAttributes!(Args[0]) & FunctionAttribute.pure_) != 0;
        }
        else
        {
            // BUG:  Should check this for delegates too, but std.traits
            //       apparently doesn't allow this.  isPure is irrelevant
            //       for delegates, at least for now since shared delegates
            //       don't work.
            private enum bool isPure = false;
        }

    }
    else
    {
        // We already know that we can't execute aliases in @safe code, so
        // just put a dummy value here.
        private enum bool isPure = false;
    }


    /**
    The return type of the function called by this `Task`.  This can be
    `void`.
    */
    alias ReturnType = typeof(fun(_args));

    static if (!is(ReturnType == void))
    {
        static if (is(typeof(&fun(_args))))
        {
            // Ref return.
            ReturnType* returnVal;

            ref ReturnType fixRef(ReturnType* val)
            {
                return *val;
            }

        }
        else
        {
            ReturnType returnVal;

            ref ReturnType fixRef(ref ReturnType val)
            {
                return val;
            }
        }
    }

    private void enforcePool()
    {
        import std.exception : enforce;
        enforce(this.pool !is null, "Job not submitted yet.");
    }

    static if (Args.length > 0)
    {
        private this(Args args)
        {
            _args = args;
        }
    }

    // Work around DMD bug https://issues.dlang.org/show_bug.cgi?id=6588,
    // allow immutable elements.
    static if (allSatisfy!(isAssignable, Args))
    {
        typeof(this) opAssign(typeof(this) rhs)
        {
            foreach (i, Type; typeof(this.tupleof))
            {
                this.tupleof[i] = rhs.tupleof[i];
            }
            return this;
        }
    }
    else
    {
        @disable typeof(this) opAssign(typeof(this) rhs);
    }
}

// Calls `fpOrDelegate` with `args`.  This is an
// adapter that makes `Task` work with delegates, function pointers and
// functors instead of just aliases.
ReturnType!F run(F, Args...)(F fpOrDelegate, ref Args args)
{
    return fpOrDelegate(args);
}

/**
Creates a malloc'd `Task` that calls an alias.  This may be executed
via `Task.executeInNewThread` or by submitting to a
$(REF TaskPool, std,parallelism).  A globally accessible instance of
`TaskPool` is provided by $(REF taskPool, std,parallelism).
Returns:  A pointer to the `Task`.
Example:
---
// Read two files into memory at the same time.
import std.file;
void main()
{
    // Create and execute a Task for reading
    // foo.txt.
    auto file1Task = task!read("foo.txt");
    file1Task.executeInNewThread();
    // Read bar.txt in parallel.
    auto file2Data = read("bar.txt");
    // Get the results of reading foo.txt.
    auto file1Data = file1Task.yieldForce;
}
---
---
// Sorts an array using a parallel quick sort algorithm.
// The first partition is done serially.  Both recursion
// branches are then executed in parallel.
//
// Timings for sorting an array of 1,000,000 doubles on
// an Athlon 64 X2 dual core machine:
//
// This implementation:               176 milliseconds.
// Equivalent serial implementation:  280 milliseconds
void parallelSort(T)(T[] data)
{
    // Sort small subarrays serially.
    if (data.length < 100)
    {
         std.algorithm.sort(data);
         return;
    }
    // Partition the array.
    swap(data[$ / 2], data[$ - 1]);
    auto pivot = data[$ - 1];
    bool lessThanPivot(T elem) { return elem < pivot; }
    auto greaterEqual = partition!lessThanPivot(data[0..$ - 1]);
    swap(data[$ - greaterEqual.length - 1], data[$ - 1]);
    auto less = data[0..$ - greaterEqual.length - 1];
    greaterEqual = data[$ - greaterEqual.length..$];
    // Execute both recursion branches in parallel.
    auto recurseTask = task!parallelSort(greaterEqual);
    taskPool.put(recurseTask);
    parallelSort(less);
    recurseTask.yieldForce;
}
---
*/
auto task(alias fun, Args...)(Args args)
{
    alias T = Task!(fun, Args);
    auto task = cast(T*) calloc(T.sizeof, 1);
    *task = T(args);
    return task;
}

/**
Creates a malloc'd `Task` that calls a function pointer, delegate, or
class/struct with overloaded opCall.
Example:
---
// Read two files in at the same time again,
// but this time use a function pointer instead
// of an alias to represent std.file.read.
import std.file;
void main()
{
    // Create and execute a Task for reading
    // foo.txt.
    auto file1Task = task(&read!string, "foo.txt", size_t.max);
    file1Task.executeInNewThread();
    // Read bar.txt in parallel.
    auto file2Data = read("bar.txt");
    // Get the results of reading foo.txt.
    auto file1Data = file1Task.yieldForce;
}
---
Notes: This function takes a non-scope delegate, meaning it can be
       used with closures.  If you can't allocate a closure due to objects
       on the stack that have scoped destruction, see `scopedTask`, which
       takes a scope delegate.
 */
auto task(F, Args...)(F delegateOrFp, Args args)
if (is(typeof(delegateOrFp(args))) && !isSafeTask!F)
{
    alias T = Task!(run, F, Args);
    auto task = cast(T*) calloc(T.sizeof, 1);
    *task = T(delegateOrFp, args);
    return task;
}

/**
Version of `task` usable from `@safe` code.  Usage mechanics are
identical to the non-@safe case, but safety introduces some restrictions:
1.  `fun` must be @safe or @trusted.
2.  `F` must not have any unshared aliasing as defined by
    $(REF hasUnsharedAliasing, std,traits).  This means it
    may not be an unshared delegate or a non-shared class or struct
    with overloaded `opCall`.  This also precludes accepting template
    alias parameters.
3.  `Args` must not have unshared aliasing.
4.  `fun` must not return by reference.
5.  The return type must not have unshared aliasing unless `fun` is
    `pure` or the `Task` is executed via `executeInNewThread` instead
    of using a `TaskPool`.
*/
@trusted auto task(F, Args...)(F fun, Args args)
if (is(typeof(fun(args))) && isSafeTask!F)
{
    alias T = Task!(run, F, Args);
    auto task = cast(T*) calloc(T.sizeof, 1);
    *task = T(fun, args);
    return task;
}

debug(dhtslib_unittest) unittest {
    auto add = (int a, int b) => a + b;
    auto t = add.task(3, 4);
    t = cast(typeof(t))t.runTask(t);
    assert(t.returnVal == 7);
    import std.exception : assertThrown;
    assertThrown(t.enforcePool);
}

/// End of Boost Licensed code

/// Wraps an htslib hts_tpool_process into a range interface. Acts as a parallel forward range.
struct Process(alias fun, R)
{
    ThreadPool * pool;

    HtsProcess proc;

    alias Args = ElementType!R;

    R range;

    alias TaskType = PointerTarget!(ReturnType!(task!(fun, Args)));
    

    TaskType frontVal;

    bool end;
    bool empty;
    
    this(R range, ThreadPool * pool) {
        this.pool = pool;
        this.range = range;
        this.proc = HtsProcess(hts_tpool_process_init(this.pool.tpool, cast(int)this.pool.threads * 2, 0));
        frontVal = frontVal.init;
        auto qsize =  hts_tpool_process_qsize(this.proc);
        if(this.range.empty) {
            this.empty = true;
        }
        for(auto i = 0; i < qsize; i++) {
            if(this.range.empty) {
                break;
            }
            this.dispatch(task!fun(this.range.front));
            this.range.popFront;
        }
        if(!this.empty) this.popFront;
    }

    this(this) {
        this.proc = proc;
    }

    /// dispatch task to thread pool
    /// non-blocking
    auto dispatch(TaskType * task) {
        task.pool = this.pool;
        auto ret = hts_tpool_dispatch(this.pool.tpool, this.proc, &TaskType.runTask, cast(void *) task);
        assert(!ret);
    }
    
    /// get front task result
    auto front() {
        return this.frontVal.returnVal;
    }

    /// get next result from task pool
    /// and submit another task if input range is not empty
    void popFront() {
        if(this.end){
            this.empty = true;
            return;
        }
        auto res = hts_tpool_next_result_wait(this.proc);
        frontVal = *cast(TaskType *)hts_tpool_result_data(res);
        hts_tpool_delete_result(res, 1);
        if(!this.range.empty) {
            this.dispatch(task!fun(range.front));
            this.range.popFront;
        }
        if(range.empty && cast(bool)hts_tpool_process_empty(this.proc)) 
            this.end = true;
    }
}

/// map a function to range in parallel
template parallelMap(fun...) {
    auto parallelMap(R)(R range, ThreadPool * pool = null){
        import std.functional : unaryFun;
        alias _fun = unaryFun!fun;
        alias _funs = AliasSeq!(_fun);
        if(pool is null) {
            assert(!(globalPool.htsPool.pool is null));
            pool = &globalPool;
        }
        return Process!(fun[0], R)(range, pool);
    }
}

debug(dhtslib_unittest) unittest {
    import std.stdio;
    import std.range : iota;
    import std.array : array;
    int[] input = iota(0, 1000).array;

    int[] result = iota(1, 1001).array;

    alias x = (x) => a + 1;
    enforceGlobalThreadPool(4);
    auto res = input.parallelMap!( x => x + 1).array;

    assert(res == result);
}

debug(dhtslib_unittest) unittest {
    int[] input;

    int[] result;

    /// test empty range
    alias x = (x) => a + 1;
    enforceGlobalThreadPool(4);
    auto res = input.parallelMap!( x => x + 1).array;

    assert(res == result);

    /// test range smaller than pool queue size
    input = [1, 2, 3];
    result = [2, 3, 4];
    
    res = input.parallelMap!( x => x + 1).array;

    assert(res == result);
}