module dhtslib.gff3d.gtfrecord;

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

import dhtslib.coordinates;

//import dhtslib.htslib.hts_log;

/** GTF Record

Format Documentation:
 *  http://gmod.org/wiki/GFF2#The_GFF2_File_Format
 *  https://useast.ensembl.org/info/website/upload/gff.html
 *  
 *  TODO: make sortable
 *  TODO: Make record builder (i.e. start with blank record and add attrs) to prep for writing
 */
struct GTF_Record
{
    private ubyte[] raw;

    /// constructor (raw ubytes)
    this(ubyte[]data){
        this.raw = data;
    }
    /// constructor (string)
    this(string data){
        this.raw = cast(ubyte[])data;
    }

    /// TODO: Not implemented; (almost) always true
    @property bool isValid()
    {
        //hts_log_trace(__FUNCTION__, format("raw.length %d", raw.length));
        return (raw.length >= 0 ? true : false);
    }

    /// Column 1: seqid (aka contig); basis for the coordinate system
    @property seqid() const { return cast(string)this.raw.splitter('\t').front; }
    /// ditto
    @property contig() const{ return cast(string)this.raw.splitter('\t').front; }

    /// Column 2: source; software, procedure, or database originating the record
    @property source() const { return cast(string)this.raw.splitter('\t').drop(1).front; }

    /// Column 3: feature type; sequence ontology (SO) defined type, or SO accession number
    @property type() const { return cast(string)this.raw.splitter('\t').drop(2).front; }

    /// Columns 4 & 5: returns Coordinate set: Obc format
    @property coordinates() const
    {
        auto start = (cast(string)this.raw.splitter('\t').drop(3).front).to!long;
        auto end = (cast(string)this.raw.splitter('\t').drop(4).front).to!long;
        return Obc(start, end);
    }
    /// Columns 4: start; 1-based integer start position of the feature
    @property start() const { return this.coordinates.start; }
    /// Column 5: end; closed coordinate integer ending nucleotide position of the feature
    @property end() const { return this.coordinates.end; }

    /// Column 6: score; float. From the standard: "the semantics of the score are ill-defined."
    /// Tragically, score can be either a float, or not present (".")
    /// Totally arbitrarily, we will represent absent as -1
    @property score() const {
        if(cast(string)this.raw.splitter('\t').drop(5).front=="."){
            return -1.0;
        }
        return (cast(string)this.raw.splitter('\t').drop(5).front).to!float;
    }

    /// Column 7: strand; '+', '-', or '.' (or '?' for relevant but unknown)
    @property strand() const {
        return cast(char)this.raw.splitter('\t').drop(6).front[0];
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
    @property phase() const {
        if(cast(string)this.raw.splitter('\t').drop(7).front=="."){
            return -1;
        }
        return (cast(string)this.raw.splitter('\t').drop(7).front).to!int;
    }

    /// Column 9: attributes; A list of ;-separated feature attributes in key=value form
    string attributes(const string field) const { return this.opIndex(field); }

    /// Provides map key lookup semantics for column 9 attributes
    string opIndex(string field) const {
        const attrs=this.raw.splitter('\t').drop(8).front.splitter(";");
        auto val = attrs    // actualy a Range of key=val
            .filter!(kv => ((cast(string) kv).strip[0 .. (cast(string) kv).strip.countUntil(' ')]) == field);
            //.front; // -- AssertError if range is empty
        if (!val.empty) return ((cast(string) val.front).strip[(cast(string) val.front).strip.countUntil(' ') + 1..$]).strip.strip("\"");
        else return "";

        /+ Alternative impl -- benchmark (also pull field ~ "=" out of the filter and combine it once upfront)
        auto vals = attrs
                    .filter!(kv => kv.startsWith(field ~ "="))
                    .map!(kv => kv[kv.countUntil('=')+1 .. $]);
        if (!vals.empty) return cast(string) val.front;
        else return "";
        +/
    }

    /// Column 9 attributes may also include a comma-sep list of tags: (key:tag)={t1,t2,t3,...}
    bool hasTag(string tagName)()
    {
        return this["tag"].splitter(",").filter!(t => t == tagName).count > 0;
    }

    /// Computed feature length
    @property length() const { return this.end - (this.start-1); }
    /// Relative start === 1
    @property relativeStart() const { return Ob(1); }
    /// Relative start === the feature length
    @property relativeEnd() const  { return Ob(this.length); }

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

        return Ob(begin + (direction * offset));
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

    string toString() const {
        return cast(string) this.raw;
    }
    /// Returns a string with the canonical "chr:start-end" representation
    @property string canonicalRepresentation() const {
        return this.seqid~":"~this.start.pos.to!string~"-"~this.end.pos.to!string;
    }
    /// Return the seqURI representation
    @property string seqURI() const {
        return format("seq:unk/%s", this.canonicalRepresentation);
    }

}
unittest{
    auto rec    = GTF_Record("chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\tID \"ENSG00000223972.5\" ; gene_id ENSG00000223972.5 ; gene_id ENSG00000223972.5 ; gene_type transcribed_unprocessed_pseudogene ; gene_name DDX11L1 ; level 2 ; havana_gene OTTHUMG00000000961.2"); // @suppress(dscanner.style.long_line)
    auto rec_neg= GTF_Record("chr1\tHAVANA\tgene\t11869\t14409\t.\t-\t.\tID \"ENSG00000223972.5\" ; gene_id ENSG00000223972.5 ; gene_id ENSG00000223972.5 ; gene_type transcribed_unprocessed_pseudogene ; gene_name DDX11L1 ; level 2 ; havana_gene OTTHUMG00000000961.2"); // @suppress(dscanner.style.long_line)

    assert(rec.seqid=="chr1");
    assert(rec.source=="HAVANA");
    assert(rec.type=="gene");
    assert(rec.start==11_869);
    assert(rec.end==14_409);
    assert(rec.score==-1.0);
    assert(rec.strand()=='+');
    assert(rec.phase==-1);
    assert(rec["ID"] == "ENSG00000223972.5");
    assert(rec["gene_id"] == "ENSG00000223972.5");
    assert(rec["gene_type"] == "transcribed_unprocessed_pseudogene");
    assert(rec["gene_name"] == "DDX11L1");
    assert(rec["level"] == "2");
    assert(rec["havana_gene"] == "OTTHUMG00000000961.2");

    assert(rec.length == 2541);
    assert(rec.relativeStart == 1);
    assert(rec.relativeEnd == 2541);

    // Test forward and backward offsets
    assert(rec.coordinateAtOffset(2) == 11_870);
    assert(rec_neg.coordinateAtOffset(2) == 14_408);

    assert(rec.coordinateAtBegin == 11_869);
    assert(rec.coordinateAtEnd   == 14_409);

    assert(rec_neg.coordinateAtBegin == 14_409);
    assert(rec_neg.coordinateAtEnd   == 11_869);

    // TODO validator
    assert(rec.isValid);
}
