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
    /// Coordinate value
    long pos;
    alias pos this;

    invariant
    {
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
static immutable CoordSystemLabels = ["zbho", "zbc", "obho", "obc"];

/// The (start, end) coordinates within a coordinate system,
/// where the type itself encodes the coordinate system details
/// (zero or one-based; half-open vs. closed)
struct Coordinates(CoordSystem cs)
{
    /// alias Basis and End enums for this Coordsystem type
    alias basetype = coordinateSystemToBasis!cs;
    alias endtype = coordinateSystemToEnd!cs;

    /// Starting coordinate
    Coordinate!basetype start;
    /// Ending coordinate
    Coordinate!basetype end;

    /// long constructor
    this(long start, long end)
    {
        this.start.pos = start;
        this.end.pos = end;
    }

    /// Coordinate constructor
    this(Coordinate!basetype start, Coordinate!basetype end){
        this.start = start;
        this.end = end;
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

        /// alias Basis and End enums for the tocs Coordsystem type
        alias tobasetype = coordinateSystemToBasis!tocs;
        alias toendtype = coordinateSystemToEnd!tocs;

        // new coordinates
        auto newStart = this.start;
        auto newEnd = this.end;

        // return same CoordinateSystem type
        static if (cs == tocs)
            return Coordinates!(tocs)(newStart, newEnd);
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
            return Coordinates!(tocs)(newStart, newEnd);
        }
    }

    /// Get string representation for printing
    string toString() const{
        return "[" ~ CoordSystemLabels[cs] ~ "] " ~ this.start.pos.to!string ~ "-" ~ this.end.pos.to!string;
    }
}

/// Represents a string version of coordinates
/// takes a string of form: chrom:start-end
struct ChromCoordinates(CoordSystem cs)
{
    /// chromosome or contig
    string chrom;

    /// Coordinates
    Coordinates!cs coords;

    /// string constructor
    this(string chrom, Coordinates!cs coords){
        this.chrom = chrom;
        this.coords = coords;
    }
    /// string constructor
    this(string region){
        import std.array : split;

        auto first = region.split(":");
        assert(first.length == 2);

        this.chrom = first[0];

        auto strcoords = first[1].split("-");
        assert(strcoords.length == 2);

        auto start = strcoords[0].to!long;
        auto end = strcoords[1].to!long;

        this.coords = Coordinates!cs(start,end);
    }

    /// Convert coordinates to another coordinate system 
    auto to(CoordSystem tocs)()
    {
        ChromCoordinates!tocs c;
        c.chrom = this.chrom;
        c.coords = this.coords.to!tocs;
        return c;
    }

    /// Get string representation for printing
    string toString() const{
        return "[" ~ CoordSystemLabels[cs] ~ "] " ~ 
            chrom ~ ":" ~ this.coords.start.pos.to!string ~ 
            "-" ~ this.coords.end.pos.to!string;
    } 
}

alias ZBHO = Coordinates!(CoordSystem.zbho);
alias OBHO = Coordinates!(CoordSystem.obho);
alias ZBC = Coordinates!(CoordSystem.zbc);
alias OBC = Coordinates!(CoordSystem.obc);

alias ZeroBasedHalfOpen = Coordinates!(CoordSystem.zbho);
alias OneBasedHalfOpen = Coordinates!(CoordSystem.obho);
alias ZeroBasedClosed = Coordinates!(CoordSystem.zbc);
alias OneBasedClosed = Coordinates!(CoordSystem.obc);

alias ChromZBHO = ChromCoordinates!(CoordSystem.zbho);
alias ChromOBHO = ChromCoordinates!(CoordSystem.obho);
alias ChromZBC = ChromCoordinates!(CoordSystem.zbc);
alias ChromOBC = ChromCoordinates!(CoordSystem.obc);

alias ChromZeroBasedHalfOpen = ChromCoordinates!(CoordSystem.zbho);
alias ChromOneBasedHalfOpen = ChromCoordinates!(CoordSystem.obho);
alias ChromZeroBasedClosed = ChromCoordinates!(CoordSystem.zbc);
alias ChromOneBasedClosed = ChromCoordinates!(CoordSystem.obc);


debug(dhtslib_unittest) unittest
{
    import std.stdio;
    auto c0 = Coordinates!(CoordSystem.zbho)(0, 100);
    assert(c0.size == 100);

    auto c1 = c0.to!(CoordSystem.zbc);
    auto c2 = c0.to!(CoordSystem.obc);
    auto c3 = c0.to!(CoordSystem.obho);
    auto c4 = c0.to!(CoordSystem.zbho);
    
    assert(c1 == ZBC(0, 99));
    assert(c2 == OBC(1, 100));
    assert(c3 == OBHO(1, 101));
    assert(c4 == ZBHO(0, 100));
    
    writeln(c0);
    writeln(c1);
    writeln(c2);
    writeln(c3);
    // ...
}

debug(dhtslib_unittest) unittest
{
    import std.stdio;
    auto c0 = ChromCoordinates!(CoordSystem.zbho)("chrom1:0-100");
    assert(c0.coords.size == 100);

    auto c1 = c0.to!(CoordSystem.zbc);
    auto c2 = c0.to!(CoordSystem.obc);
    auto c3 = c0.to!(CoordSystem.obho);
    auto c4 = c0.to!(CoordSystem.zbho);
    
    assert(c1 == ChromZBC("chrom1:0-99"));
    assert(c2 == ChromOBC("chrom1:1-100"));
    assert(c3 == ChromOBHO("chrom1:1-101"));
    assert(c4 == ChromZBHO("chrom1:0-100"));
    
    writeln(c0);
    writeln(c1);
    writeln(c2);
    writeln(c3);
    // ...
}

debug(dhtslib_unittest) unittest
{
    import std.stdio;
    auto c0 = Coordinates!(CoordSystem.zbho)(ZB(0), ZB(100));
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
    auto c0 = Coordinates!(CoordSystem.obho)(OB(1), OB(101));
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