/** Coordinates and Coordinate Systems

    STATUS: Experimental

    Coordinates include `start` and `end`, but no reference sequence id (chr, contig, etc.)

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

enum Based
{
    zero = 0,
    one
}

enum End
{
    open = 0,
    closed
}

/// Coordinate type for a single position with in a Coordinate set,
/// where the type itself encodes the coordinate system details
/// (zero or one-based)
struct Coordinate(Based bs)
{
    long pos;
    alias pos this;

    invariant
    {
        // in one based systems, pos cannot be zero
        static if(bs == Based.one)
        {
            assert(pos != 0);
        }
    }
    
    /// Convert coordinate to another based system 
    auto to(Based tobs)()
    {
        // return same Base type
        static if (bs == tobs) return Coordinate!bs(pos);
        // one to zero base, decrement
        else static if (tobs == Based.zero) return Coordinate!bs(pos - 1);
        // zero to one base, increment
        else static if (tobs == Based.one) return Coordinate!bs(pos + 1);
        else static assert(0, "Coordinate Type error");
    }
}

/// template to convert Based, End enum combination to 
/// respective CoordinateSystem enum
template getCoordinateSystem(Based bs, End es){
    static if(bs == Based.zero){
        static if(es == End.open)
            enum CoordSystem getCoordinateSystem = CoordSystem.zbho;
        else static if(es == End.closed)
            enum CoordSystem getCoordinateSystem = CoordSystem.zbc;
    }
    else static if(bs == Based.one){
        static if(es == End.open)
            enum CoordSystem getCoordinateSystem = CoordSystem.obho;
        else static if(es == End.closed)
            enum CoordSystem getCoordinateSystem = CoordSystem.obc;
    }
}

// just sanity checking
static assert(getCoordinateSystem!(Based.zero, End.open) == CoordSystem.zbho);
static assert(getCoordinateSystem!(Based.zero, End.closed) == CoordSystem.zbc);
static assert(getCoordinateSystem!(Based.one, End.open) == CoordSystem.obho);
static assert(getCoordinateSystem!(Based.one, End.closed) == CoordSystem.obc);

/// template to convert CoordinateSystem enum to 
/// respective Based enum
template coordinateSystemToBased(CoordSystem cs){
    static if(cs == CoordSystem.zbho){
        enum Based coordinateSystemToBased = Based.zero;
    }
    else static if(cs == CoordSystem.zbc){
        enum Based coordinateSystemToBased = Based.zero;
    }
    else static if(cs == CoordSystem.obho){
        enum Based coordinateSystemToBased = Based.one;
    }
    else static if(cs == CoordSystem.obc){
        enum Based coordinateSystemToBased = Based.one;
    }
}

// just sanity checking
static assert(coordinateSystemToBased!(CoordSystem.zbho) == Based.zero);
static assert(coordinateSystemToBased!(CoordSystem.zbc) == Based.zero);
static assert(coordinateSystemToBased!(CoordSystem.obho) == Based.one);
static assert(coordinateSystemToBased!(CoordSystem.obc) == Based.one);

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


enum CoordSystem
{
    zbho = 0,   /// zero-based, half-open (BED, BAM)
    zbc,        /// zero-based, closed (HTSlib.faidx_fetch_seq)
    obho,       /// one-based, half-open
    obc         /// one-based, closed (GFF3, VCF [not BCF], SAM [not BAM])
}

/// The (start, end) coordinates within a coordinate system,
/// where the type itself encodes the coordinate system details
/// (zero or one-based; half-open vs. closed)
struct Coordinates(CoordSystem cs)
{
    /// alias Based and End enums for this Coordsystem type
    alias basetype = coordinateSystemToBased!cs;
    alias endtype = coordinateSystemToEnd!cs;

    Coordinate!basetype start;
    Coordinate!basetype end;

    this(long start, long end)
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

    // Return the size of the interval spanned by start and end
    long size()
    {
        static if (cs == CoordSystem.zbho || cs == CoordSystem.obho)
            return end - start;
        else static if (cs == CoordSystem.zbc || cs == CoordSystem.obc)
            return end - start + 1;
        else
            static assert(0, "Coordinate Type error");
    }
    
    /// Convert coordinates to another coordinate system 
    auto to(CoordSystem tocs)()
    {

        /// alias Based and End enums for the tocs Coordsystem type
        alias tobasetype = coordinateSystemToBased!tocs;
        alias toendtype = coordinateSystemToEnd!tocs;

        // new coordinates
        auto newStart = this.start;
        auto newEnd = this.end;

        // return same CoordinateSystem type
        static if (cs == tocs)
            return Coordinates!(tocs)(newStart, newEnd);
        
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
        return Coordinates!(tocs)(newStart, newEnd);
    }
}
unittest
{
    import std.stdio;
    auto c0 = Coordinates!(CoordSystem.zbho)(0, 100);
    assert(c0.size == 100);

    auto c1 = c0.to!(CoordSystem.zbc);
    auto c2 = c0.to!(CoordSystem.obc);
    auto c3 = c0.to!(CoordSystem.obho);
    
    assert(c1 == Coordinates!(CoordSystem.zbc)(0, 99));
    assert(c2 == Coordinates!(CoordSystem.obc)(1, 100));
    assert(c3 == Coordinates!(CoordSystem.obho)(1, 101));
    
    // ...
}
