/** Coordinates and Coordinate Systems

    Interval include `start` and `end`, but no reference sequence id (chr, contig, etc.)

    The `Coordinate` type is templated on `CoordSystem` enum, so that the actual coordinate
    system in use -- and by this, we mean zero- or one-based, and half-open vs. closed --
    is encoded within the type itself. e.g. `Coordinate!(CoordSystem.zbho)(0, 100)`

    In this way, dhtslib functions that take type `Coordinate` can enforce safety checks,
    or optionally, interconversions, to avoid off-by-one errors when dealing with different
    systems, as is common in different HTS/NGS file formats.

    In general, zero-based half-open (ZBHO) coordinates are easy to compute on, whereas
    many data sources intended to be read and understood by humans (including NCBI) are
    one-based closed systems. You may enjoy reading the following articles:

    https://www.biostars.org/p/84686/
    https://www.biostars.org/p/6373/
    http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/
    http://genomewiki.ucsc.edu/index.php/Coordinate_Transforms

    Zero-based, half open:  BED, BAM
    Zero-based, closed:     HTSlib function faidx_fetch_seq , strangely
    One-based, half open:   ?
    One-based, closed:      GFF3, SAM (text; not BAM), VCF (text; not BCF)
*/

/*

Translation table for my reference:
    zbho    zbc     obho    obc
zbho -       0,-1   +1,+1   +1, 0
zbc  0,+1   -       +1,+2   +1,+1
obho -1,-1  -1,-2   -        0,-1
obc  -1,0   -1,-1    0,+1    -

*/

module dhtslib.coordinates;
import std.traits : isIntegral;
import std.conv : to;

/// Represents 0-based vs 1-based coordinate types
enum Basis
{
    zero = 0,
    one
}

/// Represents whether a coordinate set's end 
/// coordinate is open or closed
enum End
{
    open = 0,
    closed
}

/// Coordinate type for a single position with in a Coordinate set,
/// where the type itself encodes the coordinate system details
/// (zero or one-based)
struct Coordinate(Basis bs)
{
    @safe @nogc nothrow pure:
    /// Coordinate value
    long pos;
    alias pos this;

    invariant
    {
        assert(pos >= 0);
        // in one based systems, pos cannot be zero
        static if(bs == Basis.one)
        {
            assert(pos != 0);
        }
    }
    
    /// Convert coordinate to another based system 
    auto to(Basis tobs)()
    {
        // return same Base type
        static if (bs == tobs) return Coordinate!bs(pos);
        // one to zero base, decrement
        else static if (tobs == Basis.zero) return Coordinate!bs(pos - 1);
        // zero to one base, increment
        else static if (tobs == Basis.one) return Coordinate!bs(pos + 1);
        else static assert(0, "Coordinate Type error");
    }

    /// Convert coordinate to another based system using shortcuts
    auto to(T: ZB)()
    {
        return this.to!(Basis.zero);
    }

    auto to(T: OB)()
    {
        return this.to!(Basis.one);
    }

    /// make a new coordinate with a value of this.pos + off
    /// this.pos + off
    auto offset(T)(T off)
    if(isIntegral!T)
    {
        return Coordinatse!(bs)(cast(long)(this.pos + off));
    }
}

alias ZB = Coordinate!(Basis.zero);
alias OB = Coordinate!(Basis.one);

alias ZeroBased = Coordinate!(Basis.zero);
alias OneBased = Coordinate!(Basis.one);

/// template to convert Basis, End enum combination to 
/// respective CoordinateSystem enum
template getCoordinateSystem(Basis bs, End es){
    static if(bs == Basis.zero){
        static if(es == End.open)
            enum CoordSystem getCoordinateSystem = CoordSystem.zbho;
        else static if(es == End.closed)
            enum CoordSystem getCoordinateSystem = CoordSystem.zbc;
    }
    else static if(bs == Basis.one){
        static if(es == End.open)
            enum CoordSystem getCoordinateSystem = CoordSystem.obho;
        else static if(es == End.closed)
            enum CoordSystem getCoordinateSystem = CoordSystem.obc;
    }
}

// just sanity checking
static assert(getCoordinateSystem!(Basis.zero, End.open) == CoordSystem.zbho);
static assert(getCoordinateSystem!(Basis.zero, End.closed) == CoordSystem.zbc);
static assert(getCoordinateSystem!(Basis.one, End.open) == CoordSystem.obho);
static assert(getCoordinateSystem!(Basis.one, End.closed) == CoordSystem.obc);

/// template to convert CoordinateSystem enum to 
/// respective Basis enum
template coordinateSystemToBasis(CoordSystem cs){
    static if(cs == CoordSystem.zbho){
        enum Basis coordinateSystemToBasis = Basis.zero;
    }
    else static if(cs == CoordSystem.zbc){
        enum Basis coordinateSystemToBasis = Basis.zero;
    }
    else static if(cs == CoordSystem.obho){
        enum Basis coordinateSystemToBasis = Basis.one;
    }
    else static if(cs == CoordSystem.obc){
        enum Basis coordinateSystemToBasis = Basis.one;
    }
}

// just sanity checking
static assert(coordinateSystemToBasis!(CoordSystem.zbho) == Basis.zero);
static assert(coordinateSystemToBasis!(CoordSystem.zbc) == Basis.zero);
static assert(coordinateSystemToBasis!(CoordSystem.obho) == Basis.one);
static assert(coordinateSystemToBasis!(CoordSystem.obc) == Basis.one);

/// template to convert CoordinateSystem enum to 
/// respective End enum
template coordinateSystemToEnd(CoordSystem cs){
    static if(cs == CoordSystem.zbho){
        enum End coordinateSystemToEnd = End.open;
    }
    else static if(cs == CoordSystem.zbc){
        enum End coordinateSystemToEnd = End.closed;
    }
    else static if(cs == CoordSystem.obho){
        enum End coordinateSystemToEnd = End.open;
    }
    else static if(cs == CoordSystem.obc){
        enum End coordinateSystemToEnd = End.closed;
    }
}

// just sanity checking
static assert(coordinateSystemToEnd!(CoordSystem.zbho) == End.open);
static assert(coordinateSystemToEnd!(CoordSystem.zbc) == End.closed);
static assert(coordinateSystemToEnd!(CoordSystem.obho) == End.open);
static assert(coordinateSystemToEnd!(CoordSystem.obc) == End.closed);

/// Coordsytem types
enum CoordSystem
{
    zbho = 0,   /// zero-based, half-open (BED, BAM)
    zbc,        /// zero-based, closed (HTSlib.faidx_fetch_seq)
    obho,       /// one-based, half-open
    obc         /// one-based, closed (GFF3, VCF [not BCF], SAM [not BAM])
}

/// Labels for each CoordSystem Type
enum CoordSystemLabels = __traits(allMembers, CoordSystem);

/// The (start, end) coordinates within a coordinate system,
/// where the type itself encodes the coordinate system details
/// (zero or one-based; half-open vs. closed)
struct Interval(CoordSystem cs)
{
    @safe nothrow pure:

    /// alias Basis and End enums for this Coordsystem type
    alias basetype = coordinateSystemToBasis!cs;
    alias endtype = coordinateSystemToEnd!cs;

    /// Starting coordinate
    Coordinate!basetype start;
    /// Ending coordinate
    Coordinate!basetype end;

    /// long constructor
    this(long start, long end) @nogc
    {
        this.start.pos = start;
        this.end.pos = end;
    }

    invariant
    {
        // In half-open systems, ensure start strictly less than end,
        // unless the zero-length range is truly desired, in which case use (0,0)
        static if (cs == CoordSystem.zbho || cs == CoordSystem.obho)
            assert(this.start < this.end || (this.start == 0 && this.end == 0));
        else
            assert(this.start <= this.end);
    }

    /// Return the size of the interval spanned by start and end
    long size() @nogc
    {
        static if (cs == CoordSystem.zbho || cs == CoordSystem.obho)
            return end - start;
        else static if (cs == CoordSystem.zbc || cs == CoordSystem.obc)
            return end - start + 1;
        else
            static assert(0, "Coordinate Type error");
    }
    
    /// Convert coordinates to another coordinate system 
    auto to(CoordSystem tocs)() @nogc
    {

        /// alias Basis and End enums for the tocs Coordsystem type
        alias tobasetype = coordinateSystemToBasis!tocs;
        alias toendtype = coordinateSystemToEnd!tocs;

        // new coordinates
        auto newStart = this.start;
        auto newEnd = this.end;

        // return same CoordinateSystem type
        static if (cs == tocs)
            return Interval!(tocs)(newStart, newEnd);
        else{
            // convert coordinates to new base type
            newStart = newStart.to!tobasetype;
            newEnd = newEnd.to!tobasetype;

            // if going between end types
            static if (endtype != toendtype){
                // open to closed end, decrement end
                // closed to open end, increment end
                static if(toendtype == End.closed){
                    newEnd--;
                }else{
                    newEnd++;
                }
            }
            return Interval!(tocs)(newStart, newEnd);
        }
    }

    /// Convert coordinate to another based system using shortcuts
    auto to(T: ZBHO)() @nogc
    {
        return this.to!(CoordSystem.zbho);
    }

    /// Convert coordinate to another based system using shortcuts
    auto to(T: OBHO)() @nogc
    {
        return this.to!(CoordSystem.obho);
    }

    /// Convert coordinate to another based system using shortcuts
    auto to(T: ZBC)() @nogc
    {
        return this.to!(CoordSystem.zbc);
    }
    
    /// Convert coordinate to another based system using shortcuts
    auto to(T: OBC)() @nogc
    {
        return this.to!(CoordSystem.obc);
    }

    /// make a new coordinate pair with a value of 
    /// this.start + off and this.end + off
    auto offset(T)(T off) @nogc
    if(isIntegral!T)
    {
        return Interval!(cs)(cast(long)(this.start + off), cast(long)(this.end + off));
    }

    /// intersection of two regions
    auto intersectImpl(Interval!cs other) @nogc
    {
        if(!this.isOverlap(other)){
            return Interval!(cs).init;
        }
        return Interval!(cs)(
            this.getMaxStart(other),
            this.getMinEnd(other)
            );
    }

    /// union of two regions
    auto unionImpl(Interval!cs other) @nogc
    {
        if(!this.isOverlap(other)){
            return Interval!(cs).init;
        }
        return Interval!(cs)(
            this.getMinStart(other),
            this.getMaxEnd(other)
            );
    }

    auto isOverlap(Interval!cs other) @nogc
    {
        static if(endtype == End.closed){
            return this.getMaxStart(other) <= this.getMinEnd(other);
        }else{
            return this.getMaxStart(other) < this.getMinEnd(other);
        }
    }

    auto getMinStart(Interval!cs other) @nogc
    {
        return this.start < other.start ? this.start : other.start;
    }

    auto getMaxStart(Interval!cs other) @nogc
    {
        return this.start > other.start ? this.start : other.start;
    }

    auto getMinEnd(Interval!cs other) @nogc
    {
        return this.end < other.end ? this.end : other.end;
    }
    
    auto getMaxEnd(Interval!cs other) @nogc
    {
        return this.end > other.end ? this.end : other.end;
    }

    /// set operators for intersect, union, and difference
    auto opBinary(string op)(Interval!cs other) @nogc
    {
        static if(op == "|") return unionImpl(other);
        else static if(op == "&") return intersectImpl(other);
        else static assert(0,"Operator "~op~" not implemented");
    }

    /// Get string representation for printing
    string toString() const{
        return "[" ~ CoordSystemLabels[cs] ~ "] " ~ this.start.pos.to!string ~ "-" ~ this.end.pos.to!string;
    }
}

alias ZBHO = Interval!(CoordSystem.zbho);
alias OBHO = Interval!(CoordSystem.obho);
alias ZBC = Interval!(CoordSystem.zbc);
alias OBC = Interval!(CoordSystem.obc);

alias ZeroBasedHalfOpen = Interval!(CoordSystem.zbho);
alias OneBasedHalfOpen = Interval!(CoordSystem.obho);
alias ZeroBasedClosed = Interval!(CoordSystem.zbc);
alias OneBasedClosed = Interval!(CoordSystem.obc);


debug(dhtslib_unittest) unittest
{
    auto c0 = Interval!(CoordSystem.zbho)(0, 100);
    assert(c0.size == 100);

    auto c1 = c0.to!(CoordSystem.zbc);
    auto c2 = c0.to!(CoordSystem.obc);
    auto c3 = c0.to!(CoordSystem.obho);
    auto c4 = c0.to!(CoordSystem.zbho);
    
    assert(c1 == ZBC(0, 99));
    assert(c2 == OBC(1, 100));
    assert(c3 == OBHO(1, 101));
    assert(c4 == ZBHO(0, 100));
    // ...
}

debug(dhtslib_unittest) unittest
{
    auto c0 = Interval!(CoordSystem.zbho)(0, 100);
    assert(c0.size == 100);

    auto c1 = c0.to!(CoordSystem.zbc);
    auto c2 = c0.to!(CoordSystem.obc);
    auto c3 = c0.to!(CoordSystem.obho);
    auto c4 = c0.to!(CoordSystem.zbho);
    
    assert(c1 == ZBC(0, 99));
    assert(c2 == OBC(1, 100));
    assert(c3 == OBHO(1, 101));
    assert(c4 == ZBHO(0, 100));
}

debug(dhtslib_unittest) unittest
{
    auto c0 = ZBHO(0, 100);
    assert(c0.size == 100);

    auto c1 = c0.to!ZBC;
    auto c2 = c0.to!OBC;
    auto c3 = c0.to!OBHO;
    auto c4 = c0.to!ZBHO;
    
    assert(c1 == ZBC(0, 99));
    assert(c2 == OBC(1, 100));
    assert(c3 == OBHO(1, 101));
    assert(c4 == ZBHO(0, 100));
}

debug(dhtslib_unittest) unittest
{
    auto c0 = Interval!(CoordSystem.obho)(1, 101);
    assert(c0.size == 100);

    auto c1 = c0.to!(CoordSystem.zbc);
    auto c2 = c0.to!(CoordSystem.obc);
    auto c3 = c0.to!(CoordSystem.obho);
    auto c4 = c0.to!(CoordSystem.zbho);
    
    assert(c1 == ZBC(0, 99));
    assert(c2 == OBC(1, 100));
    assert(c3 == OBHO(1, 101));
    assert(c4 == ZBHO(0, 100));
    // ...
}

debug(dhtslib_unittest) unittest
{
    ZBHO c0 = ZBHO(0, 100);
    assert(c0.size == 100);

    auto c1 = c0.offset(50);
    assert(c1 == ZBHO(50, 150));
    assert((c0 & c1) == ZBHO(50, 100));
    assert((c0 | c1) == ZBHO(0, 150));

    c1 = c0.offset(99);
    assert(c1 == ZBHO(99, 199));

    assert((c0 & c1) == ZBHO(99, 100));
    assert((c0 | c1) == ZBHO(0, 199));
}

debug(dhtslib_unittest) unittest
{
    OBC c0 = OBC(1, 100);
    assert(c0.size == 100);

    auto c1 = c0.offset(50);
    assert(c1 == OBC(51, 150));

    assert((c0 & c1) == OBC(51, 100));
    assert((c0 | c1) == OBC(1, 150));

    c1 = c0.offset(99);
    assert(c1 == OBC(100, 199));
    import std.stdio;
    writeln((c0 & c1));
    assert((c0 & c1) == OBC(100, 100));
    assert((c0 | c1) == OBC(1, 199));
}

///
// debug(dhtslib_unittest) unittest
// {
//     import dhtslib.sam;
//     import htslib.hts_log;
//     import std.path : buildPath, dirName;
//     import std.string : fromStringz;
//     import std.array : array; 

//     hts_set_log_level(htsLogLevel.HTS_LOG_WARNING);
//     hts_log_info(__FUNCTION__, "Testing SAMFile & SAMRecord");
//     hts_log_info(__FUNCTION__, "Loading test file");
//     auto sam = SAMFile(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","auxf#values.sam"), 0);
    
//     auto reg = getIntervalFromString("sheila:1-10",sam.header);
//     assert(reg.tid == 0);
//     assert(reg.interval == ZBHO(0,11));
// }