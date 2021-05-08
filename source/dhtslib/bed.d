module dhtslib.bed;

import std.range : inputRangeObject, InputRangeObject;
import std.algorithm.iteration: splitter;
import std.range: drop;
import std.conv: to;

import dhtslib.coordinates;

struct BedRecord
{
    string line;
    
    @property contig() const { return this.raw.splitter('\t').front; }

    @property coordinates() const
    {
        auto start = (cast(string)this.raw.splitter('\t').drop(1).front).to!long;
        auto end = (cast(string)this.raw.splitter('\t').drop(2).front).to!long;
        return Zbho(start, end);
    }
    
    @property start() const { return this.coordinates.start; }
    
    @property end() const { return this.coordinates.end; }

    @property name() const { return this.raw.splitter('\t').drop(3).front; }

    @property score() const { return this.raw.splitter('\t').drop(4).front.to!int; }
    
    @property strand() const {
        return cast(char)this.raw.splitter('\t').drop(5).front[0];
    }

    @property thickStart() const { return this.raw.splitter('\t').drop(6).front.to!int; }

    @property thickEnd() const { return this.raw.splitter('\t').drop(7).front.to!int; }

    @property itemRGB() const { return this.raw.splitter('\t').drop(8).front.to!int; }

    @property blockCount() const { return this.raw.splitter('\t').drop(9).front.to!int; }

    @property blockSizes() const { return this.raw.splitter('\t').drop(10).front.to!int; }

    @property blockStarts() const { return this.raw.splitter('\t').drop(11).front.to!int; }

}

struct BedReader
{
    BGZFile file;
    InputRangeObject range;

    this(string fn)
    {
        this.file = BGZFile(fn);
        this.range = this.file.byLineCopy.inputRangeObject;
    }

    auto front()
    {
        return ChromCoordinates;
    }

    void popFront()
    {
        this.range.popFront;
    }

    auto empty()
    {
        return this.range.empty;
    }
}