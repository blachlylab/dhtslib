module dhtslib.md;
import dhtslib.sam: SAMRecord;
import dhtslib.htslib.hts_log;
import std.regex;

/// regex to extract MD string groups
/// ex: "11A3^G" -> [(11, "A"), (3, "^G")] 
auto MDREGEX = regex(`(\d+)([\^ATGCN]*)`, "g");

enum MDOps
{
    MATCH,
    MISMATCH,
    DEL
}

struct MDPair 
{
    int numMatches;
    string __mismatch;
    string mismatch(){
        if(isDel)
            return __mismatch[1..$];
        else
            return __mismatch;
    }
    bool isDel(){ return __mismatch[0]=='^';}
    bool isLast(){ return __mismatch.length == 0; }
    invariant 
    {
        assert(isDel && (__mismatch[0] == '^'));
        assert(!isDel && (__mismatch.length <= 1));
    }
}

auto getMDPairs(SAMRecord rec){
    struct MDPairs
    {
        
    }
    auto tag = rec["MD"];
    debug if(!tag.exists) hts_log_error(__FUNCTION__,"MD tag not present");


}

// Parse MD String into ReadAllele array
ReadAllele[] getReadAllelesFromMd(SAMRecord rec)
{
    int rpos = 0;
    int qpos = 0;
    
    // get md and cigar
    auto mdtag = rec["MD"];
    auto cigar = rec.cigar.ops.dup; // if we don't duplicate, we edit bam record

    // return empty if no cigar or md
    if (rec.b.core.n_cigar == 0)
        return [];

    if (!mdtag.exists)
        return []; //throw new Exception("MD tag not present");

    // Skip hard clips
    if (cigar[0].op == Ops.HARD_CLIP)
        cigar = cigar[1 .. $];

    // Skip soft clips
    if (cigar[0].op == Ops.SOFT_CLIP)
        qpos += cigar[0].length, cigar = cigar[1 .. $];
        
    ReadAllele[] refAlleles;
    ubyte[] seq;

    // reserve array large enough to hold all unpacked nibbles
    // if odd reserve enough for null nibble as well
    seq.reserve(rec.b.core.l_qseq + (rec.b.core.l_qseq & 1));

    // loop over all packed nibbles in bam seq and unpack 
    foreach (ubyte key; bam_get_seq(rec.b)[0 .. (rec.b.core.l_qseq >> 1) + (rec.b.core.l_qseq & 1)])
    {
        seq ~= (key & 0xF0) >> 4;
        seq ~= key & 0x0F;
    }
    
    // do regex match
    auto regexMatches = match(mdtag.toString, MDREGEX).array;

    // reserve alleles
    refAlleles.reserve(rec.b.core.l_qseq + (rec.b.core.l_qseq & 1) + regexMatches.length - 1);
    
    int mod;
    foreach (m; regexMatches)
    {
        // get allele and number of matching bases
        mod = m[1].to!int;
        auto mismatch = m[2].to!string;
        
        // keep cigar in sync with md string as we parse
        while (mod != 0)
        {
            // if matching bases greater than cigar op length
            // skip to next cigar op
            // else decrement cigar op
            if (mod >= cigar[0].length)
            {
                // record reference alleles
                refAlleles ~= iota(cigar[0].length).map!(x => 
                    ReadAllele(
                        rpos + x, 
                        qpos + x, 
                        false, 
                        false, 
                        seq[qpos + x..qpos + x + 1], 
                        bam_get_qual(rec.b)[qpos + x..qpos + x + 1])
                    ).array;
                qpos += cigar[0].length;
                rpos += cigar[0].length;
                mod -= cigar[0].length;
                cigar = cigar[1 .. $];
            }else{
                // record reference alleles
                refAlleles ~= iota(mod).map!(x => 
                    ReadAllele(
                        rpos + x, 
                        qpos + x, 
                        false, 
                        false, 
                        seq[qpos + x..qpos + x + 1], 
                        bam_get_qual(rec.b)[qpos + x..qpos + x + 1])
                    ).array;
                qpos += mod;
                rpos += mod;
                cigar[0].length = cigar[0].length - mod;
                break;
            }
            if(cigar.length == 0) break;
            // if we have an INS we must add new INS mismatch as 
            // MD string doesn't contain insertions
            if (cigar[0].op == Ops.INS)
            {
                ubyte[] ins = seq[qpos - 1 .. qpos + cigar[0].length];
                ubyte[] insbq = bam_get_qual(rec.b)[qpos - 1 .. qpos + cigar[0].length];

                refAlleles ~= ReadAllele(rpos - 1, qpos - 1, false, true, ins, insbq);

                // skip to next op
                qpos += ins.length - 1;
                cigar = cigar[1 .. $];
            }
        }
        if(mismatch.length == 0){
            break;
        }
        // if DEL
        if (mismatch.length != 1)
        {
            // create nibble allele for deletion 
            auto allele = m[2].to!string
                .map!(x => cast(ubyte) seq_nt16_table[x])
                .array;
            refAlleles ~= ReadAllele(rpos - 1, qpos - 1, true, false,
                    cast(ubyte) seq[qpos - 1] ~ allele[1 .. $], // append allele to previous reference base
                    bam_get_qual(rec.b)[qpos - 1 .. qpos + 1]);
            rpos += allele.length - 1;
            debug if(cigar[0].op != Ops.DEL) hts_log_warning(__FUNCTION__,"Expected DEL op but got " ~ cigar[0].op.to!string);
            assert(refAlleles[$ - 1].bq.length == 2);
            cigar = cigar[1 .. $];
        }
        // else if SNP
        else if (mismatch.length == 1)
        {
            refAlleles ~= ReadAllele(
                    rpos, 
                    qpos, 
                    false, 
                    false, 
                    seq[qpos..qpos + 1], 
                    bam_get_qual(rec.b)[qpos..qpos + 1]);
            rpos++;
            qpos++;
            cigar[0].length = cigar[0].length - 1;
        }
        // if cigar op length is 0, skip to next
        if (cigar[0].length == 0)
            cigar = cigar[1 .. $];
    }
    return refAlleles;
}