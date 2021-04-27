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
    long start;
    long end;

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
        static if (cs == tocs)
            return Coordinates!(tocs)(this.start, this.end);

        // from zero-based, half-open
        else static if (cs == CoordSystem.zbho && tocs == CoordSystem.zbc) {
            assert(start <= end);
            return Coordinates!tocs(this.start, this.end - 1);
        }
        else static if (cs == CoordSystem.zbho && tocs == CoordSystem.obho)
            return Coordinates!tocs(this.start + 1, this.end + 1);
        else static if (cs == CoordSystem.zbho && tocs == CoordSystem.obc) {
            assert(start < end);
            return Coordinates!tocs(this.start + 1, this.end);
        }

        // from zero-based, closed
        else static if (cs == CoordSystem.zbc && tocs == CoordSystem.zbho)
            return Coordinates!tocs(this.start, this.end + 1);
        else static if (cs == CoordSystem.zbc && tocs == CoordSystem.obho)
            return Coordinates!tocs(this.start + 1, this.end + 2);
        else static if (cs == CoordSystem.zbc && tocs == CoordSystem.obc)
            return Coordiantes!tocs(this.start + 1, this.end + 1);

        // from one-based, half-open
        else static if (cs == CoordSystem.obho && tocs == CoordSystem.zbho)
            return Coordinates!tocs(this.start - 1, this.end - 1);
        else static if (cs == CoordSystem.obho && tocs == CoordSystem.zbc)
            return Coordinates!tocs(this.start - 1, this.end - 2);
        else static if (cs == CoordSystem.obho && tocs == CoordSystem.obc)
            return Coordiantes!tocs(this.start, this.end - 1);

        // from one-based, closed
        else static if (cs == CoordSystem.obc && tocs == CoordSystem.zbho)
            return Coordinates!tocs(this.start - 1, this.end);
        else static if (cs == CoordSystem.obc && tocs == CoordSystem.zbc)
            return Coordinates!tocs(this.start - 1, this.end - 1);
        else static if (cs == CoordSystem.obc && tocs == CoordSystem.obho)
            return Coordinates!tocs(this.start, this.end + 1);

        else static assert(0, "Coordinate Type error");
    }
}
unittest
{
    auto c0 = Coordinates!(CoordSystem.zbho)(0, 100);
    assert(c0.size == 100);

    auto c1 = c0.to!(CoordSystem.obc);
    assert(c1 == Coordinates!(CoordSystem.obc)(1, 100));

    // ...
}
