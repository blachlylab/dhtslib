module dhtslib.util;

import std.string : toStringz;
import std.typecons : tuple;

import dhtslib.coordinates;
import htslib.sam : bam_name2id;
import htslib.vcf : bcf_hdr_id2name;
import htslib.hts : hts_parse_reg64;

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