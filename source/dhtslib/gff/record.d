module dhtslib.gff.record;

import std.stdio;
import std.algorithm.iteration: splitter;
import std.stdio: chunks, File;
import std.algorithm: filter,map;
import std.traits: ReturnType;
import std.array: replace;
import std.conv: to;
import std.string;
import std.file: mkdir;
import std.algorithm.searching: countUntil;
import std.range: drop, ForwardRange, chain;
import std.format : format;
import std.array : split;
import std.math : isNaN;

import dhtslib.coordinates;

enum GFFVersion
{
    GFF2 = 0,
    GTF = 0,
    GFF3 = 1
}

alias GTFRecord = GFFRecord!(GFFVersion.GTF);
alias GFF2Record = GFFRecord!(GFFVersion.GFF2);
alias GFF3Record = GFFRecord!(GFFVersion.GFF3);

/** GFF Generic Record

Format Documentation:
 *  http://gmod.org/wiki/GFF2#The_GFF2_File_Format
 *  https://useast.ensembl.org/info/website/upload/gff.html
 *   
 *  Represents a GFF record. Version selected only changes formatting of 
 *  String representations and how the record is parsed. Implementation  
 *  is mostly the same. This record has getters and setters and 
 *  should be able to be built from scratch.
 */
struct GFFRecord(GFFVersion ver)
{
    private ubyte[] raw;

    private bool unpackedFields = true;
    private bool unpackedAttributes = true;
    private string[9] fields;
    private string[string] kvmap;

    /// constructor (raw ubytes)
    this(ubyte[]data){
        this.raw = data;
        this.unpackedFields = false;
        this.unpackedAttributes = false;
    }

    /// constructor (string)
    this(string data){
        this(cast(ubyte[])data);
    }

    /// unpack fields of gff line for mutability
    private void unpack()
    {
        if(this.unpackedFields) return;
        this.fields = (cast(string) this.raw).split('\t');
        assert(this.fields.length == 9);
        this.unpackedFields = true;
    }

    /// unpack attributes of gff3 line for mutability
    private void unpackGFF3Attributes(){
        if(this.unpackedAttributes) return;
        unpack;
        auto kvpairs = this.fields[8].split(';');
        foreach (string kv; kvpairs)
        {
            auto kvfields = kv.split('=');
            assert(kvfields.length == 2);
            this.kvmap[kvfields[0]] = kvfields[0];
        }
        this.unpackedAttributes = true;
    }

    /// unpack attributes of gff2 line for mutability
    private void unpackGFF2Attributes(){
        if(this.unpackedAttributes) return;
        unpack;
        auto kvpairs = this.fields[8].split(';');
        foreach (string kv; kvpairs)
        {
            auto kvfields = kv.strip.split(' ');
            assert(kvfields.length == 2);
            this.kvmap[kvfields[0].strip] = kvfields[0].strip;
        }
        this.unpackedAttributes = true;
    }

    /// TODO: Not implemented; (almost) always true
    @property bool isValid()
    {
        //hts_log_trace(__FUNCTION__, format("raw.length %d", raw.length));
        return (raw.length >= 0 ? true : false);
    }

    /// Column 1: seqid (aka contig); basis for the coordinate system
    /// getter
    @property seqid() const
    {
        if(unpackedFields) return this.fields[0];
        return cast(string)this.raw.splitter('\t').front; 
    }

    /// Column 1: seqid (aka contig); basis for the coordinate system
    /// setter
    @property seqid(string chr)
    {
        unpack;
        this.fields[0] = chr; 
    }

    /// ditto
    @property contig() const
    {
        return seqid;
    }

    /// ditto
    @property contig(string chr)
    {
        seqid(chr); 
    }

    /// Column 2: source; software, procedure, or database originating the record
    /// getter
    @property source() const
    {
        if(unpackedFields) return this.fields[1];
        return cast(string)this.raw.splitter('\t').drop(1).front;
    }

    /// Column 2: source; software, procedure, or database originating the record
    /// setter
    @property source(string src)
    {
        unpack;
        this.fields[1] = src; 
    }

    /// Column 3: feature type; sequence ontology (SO) defined type, or SO accession number
    /// getter
    @property type() const
    {
        if(unpackedFields) return this.fields[2];
        return cast(string)this.raw.splitter('\t').drop(2).front;
    }
    
    /// Column 3: feature type; sequence ontology (SO) defined type, or SO accession number
    /// setter
    @property type(string typ)
    {
        unpack;
        this.fields[2] = typ; 
    }

    /// Columns 4 & 5: returns Coordinate set: OBC format
    /// getter
    @property coordinates() const
    {
        long start, end;
        if(unpackedFields){
            start = this.fields[3].to!long;
            end = this.fields[4].to!long;
        }else{
            start = (cast(string)this.raw.splitter('\t').drop(3).front).to!long;
            end = (cast(string)this.raw.splitter('\t').drop(4).front).to!long;
        }
        return OBC(start, end);
    }

    /// Columns 4 & 5: returns Coordinate set: OBC format
    /// setter
    @property coordinates(CoordSystem cs)(Interval!cs coords)
    {
        unpack;
        auto newCoords = coords.to!(CoordSystem.obc);
        this.fields[3] = newCoords.start.pos.to!string;
        this.fields[4] = newCoords.end.pos.to!string; 
    }

    /// Columns 4: start; 1-based integer start position of the feature
    @property start() const { return this.coordinates.start; }
    /// Column 5: end; closed coordinate integer ending nucleotide position of the feature
    @property end() const { return this.coordinates.end; }

    /// Column 6: score; float. From the standard: "the semantics of the score are ill-defined."
    /// Tragically, score can be either a float, or not present (".")
    /// Totally arbitrarily, we will represent absent as -1
    /// getter
    @property score() const
    {
        string val;
        if(unpackedFields) val = this.fields[5];
        else val = cast(string)this.raw.splitter('\t').drop(5).front;
        if(val=="."){
            return -1.0;
        }
        return val.to!float;
    }

    /// Column 6: score; float. From the standard: "the semantics of the score are ill-defined."
    /// Tragically, score can be either a float, or not present (".")
    /// Totally arbitrarily, we will represent absent as -1
    /// setter
    @property score(float s)
    {
        unpack;
        if(s.isNaN) this.fields[5] = ".";
        else this.fields[5] = s.to!string;
    }

    /// Column 7: strand; '+', '-', or '.' (or '?' for relevant but unknown)
    // getter
    @property strand() const
    {
        if(unpackedFields) return cast(char)this.fields[6][0];
        return cast(char)this.raw.splitter('\t').drop(6).front[0];
    }

    /// Column 7: strand; '+', '-', or '.' (or '?' for relevant but unknown)
    //setter
    @property strand(char s)
    {
        unpack;
        assert(s == '+' || s == '-' || s == '.');
        this.fields[6] = [s].idup;
    }

    /** Column 8: phase;
    For features of type "CDS", the phase indicates where the feature begins with
    reference to the reading frame. The phase is one of the integers 0, 1, or 2,
    indicating the number of bases that should be removed from the beginning of
    this feature to reach the first base of the next codon. In other words, a
    phase of "0" indicates that the next codon begins at the first base of the
    region described by the current line, a phase of "1" indicates that the next
    codon begins at the second base of this region, and a phase of "2" indicates
    that the codon begins at the third base of this region. This is NOT to be
    confused with the frame, which is simply start modulo 3.

    For forward strand features, phase is counted from the start field.
    For reverse strand features, phase is counted from the end field.

    The phase is REQUIRED for all CDS features.
    
    Tragically, phase can be either an integer (0, 1, 2), or not present (".")
    Totally arbitrarily, we will represent absent as -1
    **/
    @property phase() const
    {
        string val;
        if(unpackedFields) val = this.fields[7];
        else val = cast(string)this.raw.splitter('\t').drop(7).front;
        if(val=="."){
            return -1;
        }
        return val.to!int;
    }

    @property phase(long p)
    {
        unpack;
        assert(p == 0 || p == 1 || p == 2 || p == -1);
        if(p == -1) this.fields[7] = ".";
        else this.fields[7] = p.to!string;
    }

    /// Column 9: backwards compatible GFF2 group field
    @property group() const
    {
        if(unpackedAttributes){
            static if(ver == GFFVersion.GFF3) 
                return ("\t" ~ this.kvmap.byKeyValue.map!(a => a.key ~ "=" ~ a.value).join(";"));
            else
                return ("\t" ~ (this.kvmap.byKeyValue.map!(a => " " ~ a.key ~ " " ~ a.value ~ " ").join(";"))[1 .. $]);
        }
        if(unpackedFields) return this.fields[8];
        return cast(string)this.raw.splitter('\t').drop(8).front;
    }

    /// Column 9: backwards compatible GFF2 group field
    @property group(string g)
    {
        unpack;
        this.fields[8] = g;

        /// clear hashmap as values have been overwritten
        if(this.unpackedAttributes) this.kvmap.clear;
        this.unpackedAttributes = false;
    }

    /// Column 9: attributes; A list of ;-separated feature attributes in key=value form
    string attributes(const string field){ return this.opIndex(field); }

    /// Provides map key lookup semantics for column 9 attributes
    string opIndex(string field) const
    {
        if(unpackedAttributes) return this.kvmap[field];
        static if(ver != GFFVersion.GFF3){
            auto attrs=this.raw.splitter('\t').drop(8).front.splitter(";");
            auto val = attrs    // actualy a Range of key=val
                .filter!(kv => ((cast(string) kv).strip[0 .. (cast(string) kv).strip.countUntil(' ')]) == field);
                //.front; // -- AssertError if range is empty
            if (!val.empty) return ((cast(string) val.front).strip[(cast(string) val.front).strip.countUntil(' ') + 1..$]).strip.strip("\"");
            else return "";
        }
        else
        {
            auto attrs=this.raw.splitter('\t').drop(8).front.splitter(";");
            auto val = attrs    // actualy a Range of key=val
                .filter!(kv => cast(string)(kv[0 .. kv.countUntil('=')]) == field);
                //.front; // -- AssertError if range is empty
            if (!val.empty) return cast(string) (val.front[val.front.countUntil('=')+1..$]);
            else return "";
        }
        
        /+ Alternative impl -- benchmark (also pull field ~ "=" out of the filter and combine it once upfront)
        auto vals = attrs
                    .filter!(kv => kv.startsWith(field ~ "="))
                    .map!(kv => kv[kv.countUntil('=')+1 .. $]);
        if (!vals.empty) return cast(string) val.front;
        else return "";
        +/
    }

    /// Provides map key assignment semantics for column 9 attributes
    void opIndexAssign(string val, string field)
    {
        static if(ver == GFFVersion.GFF3) unpackGFF3Attributes;
        else unpackGFF2Attributes;
        static if(ver != GFFVersion.GFF3) val = "\"" ~ val ~ "\"";
        this.kvmap[field] = val;
    }

    /// Column 9 attributes may also include a comma-sep list of tags: (key:tag)={t1,t2,t3,...}
    bool hasTag(string tagName)()
    {
        return this["tag"].splitter(",").filter!(t => t == tagName).count > 0;
    }

    /// Computed feature length
    @property length() const { return this.end.pos - (this.start.pos - 1); }
    /// Relative start === 1
    @property relativeStart() const { return OB(1); }
    /// Relative start === the feature length
    @property relativeEnd() const { return OB(this.length); }

    /// Genomic coordinate at offset into feature, taking strandedness into account
    @property coordinateAtOffset(long offset) const
    {
        // GFF3 features are 1-based coordinates
        assert(offset > 0);
        offset--;   // necessary in a 1-based closed coordinate system
        
        // for - strand count from end; for +, ., and ? strand count from start
        immutable auto begin = (this.strand == '-' ? this.end : this.start);

        // for - strand count backward; for +, ., and ? strand count forward
        immutable auto direction = (this.strand == '-' ? -1 : 1);

        return OB(begin.pos + (direction * offset));
    }
    /// Genomic coordinate at beginning of feature, taking strandedness into account
    @property coordinateAtBegin() const
    {
        return this.coordinateAtOffset(1);
    }
    /// Genomic coordinate at end of feature, taking strandedness into account
    @property coordinateAtEnd() const
    {
        return this.coordinateAtOffset(this.length);
    }

    string toString() const
    {
        if(!unpackedFields) 
            return cast(string) this.raw;
        string ret;
        if(unpackedAttributes){
            ret = this.fields[0..8].join("\t");
            static if(ver == GFFVersion.GFF3) 
                ret ~= ("\t" ~ this.kvmap.byKeyValue.map!(a => a.key ~ "=" ~ a.value).join(";"));
            else
                ret ~= ("\t" ~ (this.kvmap.byKeyValue.map!(a => " " ~ a.key ~ " " ~ a.value ~ " ").join(";"))[1 .. $]);
        }else{
            ret = this.fields[].join("\t");
        }
        return ret;
    }

    /// Returns a string with the canonical "chr:start-end" representation
    @property string canonicalRepresentation() const
    {
        return this.seqid~":"~this.start.pos.to!string~"-"~this.end.pos.to!string;
    }
    /// Return the seqURI representation
    @property string seqURI() const {
        return format("seq:unk/%s", this.canonicalRepresentation);
    }

}