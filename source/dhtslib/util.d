module dhtslib.util;

import std.string : toStringz;
import std.typecons : tuple;

import core.stdc.stdlib : malloc;
import core.stdc.stdio : SEEK_SET;
import dhtslib.coordinates;
import dhtslib.memory;
import htslib;

/// Returns tuple of String and ZBHO coordinates
/// representing input string. Supports htslib coordinate strings.
/// i.e chr1:1-10
auto getIntervalFromString(string region){
    ZBHO coords;
    auto regStr = toStringz(region);
    auto ptr = hts_parse_reg64(regStr,&coords.start.pos,&coords.end.pos);
    if(!ptr){
        throw new Exception("Region string could not be parsed");
    }
    auto contig = region[0..ptr - regStr];
    return tuple!("contig","interval")(contig,coords);
}

debug(dhtslib_unittest) unittest
{
    auto reg = getIntervalFromString("chrom1:0-100");
    assert(reg.contig == "chrom1");
    auto c0 = reg.interval;
    assert(c0 == Interval!(CoordSystem.zbho)(0, 100));
}

/// TODO: complete getIntervalFromString with the ability to check headers

// auto getIntervalFromString(Header)(string region, Header h)
//     if(is(Header == SAMHeader) || is(Header == VCFHeader))
// {
//     ZBHO coords;
//     int tid;
//     auto regStr = toStringz(region);
//     auto flag =  HTS_PARSE_FLAGS.HTS_PARSE_THOUSANDS_SEP;
//     static if(is(Header == SAMHeader)){
//         auto ptr = hts_parse_region(regStr,&tid, &coords.start.pos, &coords.end.pos, cast(hts_name2id_f * )&bam_name2id, h.h, flag);
//     } else static if(is(Header == VCFHeader)){
//         auto ptr = hts_parse_region(regStr,&tid, &coords.start.pos, &coords.end.pos, &bcf_hdr_id2name, h.hdr, flag);
//     }
//     if(!ptr){
//         throw new Exception("Region string could not be parsed");
//     }
//     return tuple!("tid","interval")(tid, coords);
// }


/**
    All functions below take a SafeHtslibPtr alias as input but return
    a raw htslib pointer. The returned pointer should be wrapped in
    a SafeHtslibPtr to maintain safety.

    (We cannot return a new SafeHtslibPtr from these functions because of 
    scope rules with dip1000)
*/

pragma(inline, true){

    /// deep copy an hts_itr_t
    /// htslib really needs a copy/dup function for this
    auto copyHtsItr(HtsItr itr) 
    {
        /// init
        auto newItr = cast(hts_itr_t *) malloc(hts_itr_t.sizeof);
        /// bulk copy data
        *newItr = *itr;

        /// initialize and copy region list
        /// if it is not null
        if(newItr.reg_list){
            /// initialize the region lists
            newItr.reg_list = cast(hts_reglist_t *) malloc(itr.n_reg * hts_reglist_t.sizeof);
            newItr.reg_list[0 .. newItr.n_reg] = itr.reg_list[0 .. itr.n_reg];
            /// for each list
            for(auto i=0; i < newItr.n_reg; i++)
            {
                /// copy all intervals
                newItr.reg_list[i].intervals = cast(hts_pair_pos_t *) malloc(itr.reg_list[i].count * hts_pair_pos_t.sizeof);
                newItr.reg_list[i].intervals[0 .. newItr.reg_list[i].count] = itr.reg_list[i].intervals[0 .. itr.reg_list[i].count];
            }
        }

        /// initialize and copy bins list
        newItr.bins.a = cast(int *) malloc(itr.bins.m * int.sizeof);
        assert(newItr.bins.m >= newItr.bins.n);
        newItr.bins.a[0 .. newItr.bins.m] = itr.bins.a[0 .. itr.bins.m];

        /// initialize and copy off list
        newItr.off = cast(hts_pair64_max_t *) malloc(itr.n_off * hts_pair64_max_t.sizeof);
        newItr.off[0 .. newItr.n_off] = itr.off[0 .. itr.n_off];
        return newItr;
    }

    /// copy a htsFile
    auto copyHtsFile(HtsFile fp){
        // collect offset of the file
        // we would like to copy
        off_t offset;
        if(fp.is_bgzf) offset = bgzf_tell(fp.fp.bgzf);
        else offset = htell(fp.fp.hfile);
        return copyHtsFile(fp, offset);
    }

    /// copy a htsFile with custom offset
    auto copyHtsFile(HtsFile fp, off_t offset)
    {
        auto newfp = hts_open(fp.fn, cast(immutable(char)*) "r");
        // Need to seek to current offset
        if (newfp.is_bgzf)
            bgzf_seek(newfp.fp.bgzf, offset, SEEK_SET);
        else
            hseek(newfp.fp.hfile, offset, SEEK_SET);
        return newfp;
    }

    /// copy a BGZF
    auto copyBgzf(Bgzf fp, const(char) * fn)
    {
        // collect offset of the file
        // we would like to copy
        auto offset = bgzf_tell(fp);
        return copyBgzf(fp, offset, fn);
    }

    /// copy a BGZF with custom offset
    auto copyBgzf(Bgzf fp, off_t offset, const(char)* fn)
    {
        auto newfp = bgzf_open(fn, cast(immutable(char)*) "r");
        // Need to seek to current offset
        bgzf_seek(newfp, offset, SEEK_SET);
        return newfp;
    }


    /// Copy a kstring_t
    auto copyKstring(Kstring str)
    {
        auto newStr = cast(kstring_t*)malloc(kstring_t.sizeof);
        newStr.l = str.l;
        newStr.m = str.m;
        newStr.s = cast(char*)malloc(str.m);
        newStr.s[0..newStr.m] = str.s[0..str.m];
        return newStr;
    }

    /// create a kstring_t *
    auto initKstring()
    {
        return cast(kstring_t*)malloc(kstring_t.sizeof);
    }
}