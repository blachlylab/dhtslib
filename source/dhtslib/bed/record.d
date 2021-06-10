module dhtslib.bed.record;

import std.range : inputRangeObject, InputRangeObject;
import std.algorithm: splitter, map;
import std.range: drop, enumerate;
import std.conv: to;
import std.array: join, split;

import dhtslib.coordinates;
import dhtslib.bgzf;
import dhtslib.tabix;

struct RGB
{
    ubyte red;
    ubyte green;
    ubyte blue;
}

/// Represents a record in a bed file.
/// Based on UCSC format and methods are derived from
/// UCSC's specifications.
struct BedRecord
{
    private string line;

    private bool unpacked = true;
    private string[] fields;

    /// string ctor
    this(string line)
    {
        this.line = line;
        this.unpacked = false;
    }

    /// unpack fields of bed line for mutability
    private void unpack()
    {
        if(unpacked) return;
        this.fields = this.line.split("\t");
        this.unpacked = true;
    }

    /// column 1: The name of the chromosome or scaffold.
    /// getter
    @property contig() const
    {
        if(unpacked) return this.fields[0];
        return this.line.splitter('\t').front;
    }

    /// column 1: The name of the chromosome or scaffold.
    /// setter
    @property contig(string val)
    {
        unpack;
        if(this.fields.length == 0) this.fields.length = 3;  
        this.fields[0] = val;
    }

    /// Columns 2 & 3 as coordinate pair, Zero-based half-open.
    /// column 2: The starting position of the feature in the chromosome or scaffold.
    /// column 3: The ending position of the feature in the chromosome or scaffold. 
    /// getter
    @property coordinates() const
    {
        if(unpacked)
            return ZBHO(this.fields[1].to!long, this.fields[2].to!long);
        auto start = (cast(string)this.line.splitter('\t').drop(1).front).to!long;
        auto end = (cast(string)this.line.splitter('\t').drop(2).front).to!long;
        return ZBHO(start, end);
    }

    /// Columns 2 & 3 as coordinate pair, Zero-based half-open.
    /// column 2: The starting position of the feature in the chromosome or scaffold.
    /// column 3: The ending position of the feature in the chromosome or scaffold. 
    /// setter
    @property coordinates(CoordSystem cs)(Coordinates!cs coords)
    {
        unpack;
        if(this.fields.length < 3) this.fields.length = 3;
        auto newCoords = coords.to!(CoordSystem.zbho);
        this.fields[1] = newCoords.start.pos.to!string;
        this.fields[2] = newCoords.end.pos.to!string; 
    }
    
    /// column 2: The starting position of the feature in the chromosome or scaffold.
    @property start() const { return this.coordinates.start; }
    
    /// column 3: The ending position of the feature in the chromosome or scaffold. 
    @property end() const { return this.coordinates.end; }

    /// column 4: Defines the name of the BED line.
    /// getter
    @property name() const
    {
        if(unpacked) return this.fields[3];
        return this.line.splitter('\t').drop(3).front; 
    }

    /// column 4: Defines the name of the BED line.
    /// setter
    @property name(string val)
    {
        unpack;
        if(this.fields.length < 4) this.fields.length = 4;
        this.fields[3] = val;
    }

    /// column 5: A score between 0 and 1000.
    /// getter
    @property score() const
    {
        if(unpacked) return this.fields[4].to!int;
        return this.line.splitter('\t').drop(4).front.to!int;
    }

    /// column 5: A score between 0 and 1000.
    /// setter
    @property score(int val)
    {
        unpack;
        if(this.fields.length < 5) this.fields.length = 5;
        this.fields[4] = val.to!string;
    }    
    
    /// column 6: Defines the strand. Either "." (=no strand) or "+" or "-".
    /// getter
    @property strand() const 
    {
        if(unpacked) return this.fields[5][0];
        return cast(char)this.line.splitter('\t').drop(5).front[0];
    }

    /// column 6: Defines the strand. Either "." (=no strand) or "+" or "-".
    /// setter
    @property strand(char val)
    {
        assert(val == '.' || val == '+'|| val == '-');
        unpack;
        if(this.fields.length < 6) this.fields.length = 6;
        this.fields[5] = [val].idup;
    }

    /// column 7: The starting position at which the feature is drawn thickly;
    /// getter
    @property thickStart() const
    {
        if(unpacked) return this.fields[6].to!int;
        return this.line.splitter('\t').drop(6).front.to!int;
    }

    /// column 7: The starting position at which the feature is drawn thickly;
    /// setter
    @property thickStart(int val)
    {
        unpack;
        if(this.fields.length < 7) this.fields.length = 7;
        this.fields[6] = val.to!string;
    }

    /// column 8: The ending position at which the feature is drawn thickly
    /// getter
    @property thickEnd() const
    {
        if(unpacked) return this.fields[7].to!int;
        return this.line.splitter('\t').drop(7).front.to!int;
    }

    /// column 8: The ending position at which the feature is drawn thickly
    /// setter
    @property thickEnd(int val)
    {
        unpack;
        if(this.fields.length < 8) this.fields.length = 8;
        this.fields[7] = val.to!string;
    }
    
    /// column 9: An RGB value of the form R,G,B (e.g. 255,0,0).
    /// getter
    @property itemRGB() const 
    {
        string str;
        if(unpacked)
            str = this.fields[8];
        else
            str = this.line.splitter('\t').drop(8).front;
        if(str == "0") return RGB.init;
        return RGB(
            str.splitter(',').front.to!ubyte,
            str.splitter(',').drop(1).front.to!ubyte,
            str.splitter(',').drop(2).front.to!ubyte
            );
    }

    /// column 9: An RGB value of the form R,G,B (e.g. 255,0,0).
    /// Setter
    @property itemRGB(RGB rgb)
    {
        unpack;
        if(this.fields.length < 9) this.fields.length = 9;
        this.fields[8] = [rgb.red, rgb.green, rgb.blue].map!(x => x.to!string).join(",");
    }

    /// column 10: The number of blocks (exons) in the BED line.
    /// getter
    @property blockCount() const
    {
        if(unpacked) return this.fields[9].to!int;
        return this.line.splitter('\t').drop(9).front.to!int;
    }

    /// column 10: The number of blocks (exons) in the BED line.
    /// setter
    @property blockCount(int count)
    {
        unpack;
        if(this.fields.length < 10) this.fields.length = 10;
        this.fields[9] = count.to!string;
    }
    
    /// column 11: A comma-separated list of the block sizes. 
    /// The number of items in this list should correspond to blockCount.
    /// getter
    @property blockSizes() const 
    {
        string str;
        if(unpacked)
            str = this.fields[10];
        else
            str = this.line.splitter('\t').drop(10).front;
        int[] arr = new int[this.blockCount];
        foreach (i, key; str.splitter(',').enumerate)
        {
            arr[i] = key.to!int;
        }
        return arr;
    }
    
    /// column 11: A comma-separated list of the block sizes. 
    /// The number of items in this list should correspond to blockCount.
    /// setter
    @property blockSizes(int[] vals)
    {
        unpack;
        if(this.fields.length < 11) this.fields.length = 11;
        this.fields[10] = vals.map!(x => x.to!string).join(",");
    }


    /// column 12: A comma-separated list of block starts. 
    /// All of the blockStart positions should be calculated relative to chromStart. 
    /// The number of items in this list should correspond to blockCount.
    /// getter
    @property blockStarts() const 
    {
        string str;
        if(unpacked)
            str = this.fields[11];
        else
            str = this.line.splitter('\t').drop(11).front;
        int[] arr = new int[this.blockCount];
        foreach (i, key; str.splitter(',').enumerate)
        {
            arr[i] = key.to!int;
        }
        return arr;
    }

    /// column 12: A comma-separated list of block starts. 
    /// All of the blockStart positions should be calculated relative to chromStart. 
    /// The number of items in this list should correspond to blockCount.
    /// setter
    @property blockStarts(int[] vals)
    {
        unpack;
        if(this.fields.length < 12) this.fields.length = 12;
        this.fields[11] = vals.map!(x => x.to!string).join(",");
    }

    /// get column idx as a string
    /// helps if your bed file isn't a ucsc bed file
    auto opIndex(ulong idx)
    {
        if(unpacked) return this.fields[idx];
        return this.line.splitter('\t').drop(idx).front;
    }

    /// set column idx as a string
    /// helps if your bed file isn't a ucsc bed file
    auto opIndexAssign(string val, ulong idx)
    {
        unpack;
        if(this.fields.length < idx) this.fields.length = idx;
        return this.fields[idx] = val;
    }

    /// return bed line
    string toString() const
    {
        if(unpacked) return this.fields.join("\t");
        return this.line;
    }
}

