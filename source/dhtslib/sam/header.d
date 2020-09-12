module dhtslib.sam.header;

import htslib.sam;
import htslib.kstring;

import core.stdc.stdlib : free;
import core.stdc.string : memcpy;

import std.string : toStringz;
import std.traits : isSomeString;

/// SAM specifications Section 1.3
/// Each header line begins with the character '@' followed by one of the
/// two-letter header record type codes defined in this section.
///
/// (leftover from when was char(3) -- the allocator evidently always puts \0 at end of char[] also, not just string
/// Note that it could also be char[2] because the compiler inserts implicit \0
/// in statically defined strings, but char[3] lets us be safer and surer.
enum RecordType : immutable(char)[2]
{
    HD = "HD",
    SQ = "SQ",
    RG = "RG",
    PG = "PG",
    CO = "CO",
}

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
    /// e.g. a sample-name lookup in Pysam is ["RG"][0]["SN"] , while in dhtslib:
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
        return sam_hdr_add_lines(this.h, lines.ptr, lines.length);
    }

    /// Add a single line to an existing header
    /// parameter kvargs is a typesafe variadic arg 
    /// https://dlang.org/spec/function.html#typesafe_variadic_functions
    auto addLine(RecordType type, const(char)[][] kvargs ...)
    {
        assert (kvargs.length %2 == 0);   // K-V pairs => even number of variadic args

        const(char)*[] varargs;
        varargs.length = kvargs.length + 1;
        for(int i=0; i < kvargs.length; i++)
            varargs[i] = kvargs[i].ptr;
        varargs[$-1] = null;  // last vararg param null signals to sam_hdr_add_line end of list

        return sam_hdr_add_line(this.h, type.ptr, varargs.ptr);
    }

    /// Return a complete line of formatted text for a given type and ID,
    /// or if no ID, first line matching type.
    ///
    /// Parameters:
    ///     * type      - enum
    ///     * id_key    - may be empty, in which case the first line matching type is returned
    ///     * id_val    - may be empty IFF id_key empty; otherwise must be value for key
    deprecated("rewrite  like valueByPos")
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

unittest
{
    import std;

    auto h = sam_hdr_init();
    auto hdr = SAMHeader(h);

    assert(!(RecordType.RG in hdr));

    //sam_hdr_add_line(h, RecordType.RG.ptr, "ID".ptr, "001".ptr, "SM".ptr, "sample1".ptr, null);
    hdr.addLine(RecordType.RG, "ID", "001", "SM", "sample1");

    assert(RecordType.RG in hdr);
}
