module dhtslib.sam.header;

import htslib.sam;
import htslib.kstring;

import core.stdc.stdlib : free;
import core.stdc.string : memcpy;

import std.conv : to;
import std.string : toStringz;
import std.traits : isSomeString;

/// SAM specifications Section 1.3
/// Each header line begins with the character '@' followed by one of the
/// two-letter header record type codes defined in this section.
enum RecordType : immutable(char)[2]
{
    HD = "HD",
    SQ = "SQ",
    RG = "RG",
    PG = "PG",
    CO = "CO",
}

/** SAMHeader encapsulates `sam_hdr_t*`
    and provides convenience wrappers to manipulate the header metadata/records.
*/
struct SAMHeader
{
    private sam_hdr_t* h;

    this(sam_hdr_t* h)
    {
        this.h = h;
    }
    // no destructor

    /* Array-like indexing */

    /// 'in' membership operator.
    /// usage: RecordType.SQ in hdr; => <bool>
    bool opBinaryRight(string op)(RecordType lhs)
    if (op == "in")
    {
        if (numRecords(lhs)) return true;
        return false;
    }

    /// For position-based lookups of key,
    /// e.g. a sample-name lookup in Pysam is ["RG"][0]["SM"] ,
    /// while in dhtslib:
    /// [RecordType.RG, 0, "SN"]
    const(char)[] opIndex(RecordType rt, size_t pos, const(char)[] key)
    {
        return this.valueByPos(rt, pos, key);
    }

    /// number of records (lines) of type e.g. SQ, RG, etc.
    size_t numRecords(RecordType rt)
    {
        return sam_hdr_count_lines(this.h, rt.ptr);
    }

    /* ==== Line level methods ==== */

    /// add multiple \n-terminated full SAM header records, eg "@SQ\tSN:foo\tLN:100"
    /// (passed line does not require \n)
    auto addLines(const(char)[] lines)
    {
        import std.algorithm.searching : maxElement;
        import std.algorithm.iteration : map;
        auto maxlen = lines.map!(x => x.length).maxElement;
        return sam_hdr_add_lines(this.h, lines.ptr, maxlen);
    }

    /// Add a single line to an existing header
    auto addLine(T...)(RecordType type, T kvargs)
    if(kvargs.length > 0 && isSomeString!(T[0]))
    {
        static assert (kvargs.length %2 == 0);   // K-V pairs => even number of variadic args
/*
        // NOTE: both (runtime) type safe variadic params, and compile-time variadic templates
        // use dynamic arrays, which cannot be passed to C variadic functions no matter what.
        // complicating this, we also need to convert to char*. The below won't work period;
        // the analogous variadic template won't work either without the mixin foolishness below.
        const(char)*[] varargs;
        varargs.length = kvargs.length + 1;
        for(int i=0; i < kvargs.length; i++)
            varargs[i] = kvargs[i].ptr;
        varargs[$-1] = null;  // last vararg param null signals to sam_hdr_add_line end of list
        
        return sam_hdr_add_line(this.h, type.ptr, varargs.ptr);
*/
        string varargMagic(size_t len)
        {
            string args = "sam_hdr_add_line(this.h, type.ptr, ";
            for(int i=0; i<len; i++)
                args ~= "toStringz(kvargs[" ~ i.to!string ~ "]), ";
            args ~= "null)";
            return args;
        }

        // if mixin result is "toStringz(kvargs[0], ..." error is:
        // Error: Using the result of a comma expression is not allowed
        //return sam_hdr_add_line(this.h, type.ptr, mixin(varargMagic(kvargs.length)) );
        return mixin(varargMagic(kvargs.length));
    }

    /// Return a complete line of formatted text for a given type and ID,
    /// or if no ID, first line matching type.
    ///
    /// Parameters:
    ///     * type      - enum
    ///     * id_key    - may be empty, in which case the first line matching type is returned
    ///     * id_val    - may be empty IFF id_key empty; otherwise must be value for key
    const(char)[] lineById(RecordType type, string id_key = "", string id_val = "")
    in (id_key.length == 0 ? id_val.length == 0 : id_val.length > 0)
    {

        kstring_t ks_line;
        
        // looking at samtools header.c sam_hrecs_Find_type_id (called by sam_hdr_find_line_id),
        // passing non-null terminated two-char char* appears safe
        auto res = sam_hdr_find_line_id(this.h, type.ptr,
                                        id_key == "" ? null : id_key.ptr,
                                        id_val == "" ? null : id_val.ptr,
                                        &ks_line);

        // 0: success, -1: no match found, -2: error
        if (res < 0)
            return "";

        char[] line;
        line.length = ks_line.l;
        memcpy(line.ptr, ks_line.s, ks_line.l);
        free(ks_line.s);
        return line;
    }

    /*
int sam_hdr_find_line_pos(sam_hdr_t *h, const char *type,
                          int pos, kstring_t *ks);
    */

    /* int sam_hdr_remove_line_id(sam_hdr_t *h, const char *type, const char *ID_key, const char *ID_value); */

    /* int sam_hdr_remove_line_pos(sam_hdr_t *h, const char *type, int position); */

    /* int sam_hdr_update_line(sam_hdr_t *h, const char *type,
        const char *ID_key, const char *ID_value, ...); */

    /* int sam_hdr_remove_except(sam_hdr_t *h, const char *type, const char *ID_key, const char *ID_value); */

    /* int sam_hdr_remove_lines(sam_hdr_t *h, const char *type, const char *id, void *rh); */

    /+
    int sam_hdr_count_lines(sam_hdr_t *h, const char *type);
    int sam_hdr_line_index(sam_hdr_t *bh, const char *type, const char *key);
    const char *sam_hdr_line_name(sam_hdr_t *bh, const char *type, int pos);
    +/

    //// //// int sam_hdr_find_tag_id(sam_hdr_t *h, const char *type, const char *ID_key, const char *ID_value, const char *key, kstring_t *ks);


    /// Return the value associated with a key for a header line identified by position
    const(char)[] valueByPos(RecordType type, size_t pos, const(char)[] key)
    in (pos <= int.max)
    in (key.length > 0)
    {
        kstring_t ks;
        auto res = sam_hdr_find_tag_pos(this.h, type.ptr, cast(int)pos, toStringz(key), &ks);
        
        // 0: success, -1: tag DNE, -2: error
        if (res < 0)
            return "";

        char[] ret;
        ret.length = ks.l;
        memcpy(ret.ptr, ks.s, ks.l);
        free(ks.s);
        return ret;
    }
}
/// Example
debug (dhtslib_unittest) unittest
{
    import std;

    auto h = sam_hdr_init();
    auto hdr = SAMHeader(h);

    assert(!(RecordType.RG in hdr));

    //sam_hdr_add_line(h, RecordType.RG.ptr, "ID".ptr, "001".ptr, "SM".ptr, "sample1".ptr, null);
    hdr.addLine(RecordType.RG, "ID", "001", "SM", "sample1");

    assert(RecordType.RG in hdr);

    auto line = hdr.lineById(RecordType.RG, "ID", "001");
    assert(line == "@RG	ID:001	SM:sample1");

    auto val = hdr.valueByPos(RecordType.RG, 0, "SM");
    assert(val == "sample1");
    assert(hdr[RecordType.RG, 0, "SM"] == "sample1");
    assert(hdr[RecordType.RG, 0, "XX"] == "");
}
