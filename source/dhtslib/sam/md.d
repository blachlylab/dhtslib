/**
Module to deal with SAM records' `MD` auxillary tag.

This tag is a string encoding mismatched and deleted reference bases,
used in conjunction with CIGAR and SEQ fields to reconstruct the bases
of the reference sequence interval to which the alignment has been mapped.
This can enable variant calling without requiring access to the entire original reference.

"For example,  a string `10A5^AC6` means from the leftmost reference base in the alignment, there are 10 matches followed by an A on the reference which is different from the aligned read base; the next 5 reference bases are matches followed by a 2bp deletion from the reference; the deleted sequence is AC; the last 6 bases are matches."

Reference: https://samtools.github.io/hts-specs/SAMtags.pdf
*/
module dhtslib.sam.md;

import dhtslib.sam.record : SAMRecord;
import dhtslib.sam.cigar;
import htslib.hts_log;
import std.regex;
import std.traits : ReturnType;
import std.conv : to;

/// regex to extract MD string groups
/// ex: "11A3^G" -> [(11, "A"), (3, "^G")] 
auto MDREGEX = regex(`(\d+)([\^ATGCN]*)`);

struct MDPair
{
    int numMatches;
    string __mismatch;

    string mismatch() const
    {
        if (isDel)
            return __mismatch[1 .. $];
        else
            return __mismatch;
    }

    bool isDel() const
    {
        return !isLast && __mismatch[0] == '^';
    }

    bool isLast() const
    {
        return __mismatch.length == 0;
    }

    invariant
    {
        if (__mismatch.length > 1)
            assert(__mismatch[0] == '^');
        if (__mismatch.length == 1)
            assert(__mismatch[0] != '^');
    }

    string toString()
    {
        return numMatches.to!string ~ __mismatch;
    }
}

/// (?) For SAM record `rec`, return ForwardRange over read's MD tag data
auto getMDPairs(SAMRecord rec)
{
    struct MDPairs
    {
        string md_string;
        ReturnType!generateMDmatches matches;

        this(string md)
        {
            md_string = md;
            matches = generateMDmatches();
        }

        auto generateMDmatches()
        {
            return matchAll(md_string, MDREGEX);
        }

        MDPair front()
        {
            return MDPair(matches.front[1].to!int, matches.front[2]);
        }

        bool empty()
        {
            return matches.empty;
        }

        void popFront()
        {
            matches.popFront;
        }

    }

    auto tag = rec["MD"];
    import std.stdio;
    debug if (!tag.exists)
        hts_log_error(__FUNCTION__, "MD tag not present");
    return MDPairs(tag.toString);
}

/// (?) Iterator yielding mismatched base in query; or '=' if equal to reference
struct MDItr
{
    ReturnType!getMDPairs mdPairs;
    int currentlen;
    string currentSeq;
    bool isLast;

    this(SAMRecord rec)
    {
        mdPairs = getMDPairs(rec);
        currentlen = mdPairs.front.numMatches;
        currentSeq = mdPairs.front.mismatch;
        isLast = mdPairs.front.isLast;
        mdPairs.popFront;
    }

    char front()
    {
        if (currentlen)
            return '=';
        else
            return currentSeq[0];
    }

    void popFront()
    {
        if (currentlen)
            currentlen--;
        else if (currentSeq.length > 1)
            currentSeq = currentSeq[1 .. $];
        else
        {
            currentlen = mdPairs.front.numMatches;
            currentSeq = mdPairs.front.mismatch;
            isLast = mdPairs.front.isLast;
            mdPairs.popFront;
        }
    }

    bool empty()
    {
        return isLast && currentlen == 0;
    }

}

debug (dhtslib_unittest) unittest
{
    import std.stdio;
    import dhtslib.sam;
    import std.array : array;
    import std.path : buildPath, dirName;

    auto bam = SAMFile(buildPath(dirName(dirName(dirName(dirName(__FILE__)))), "htslib",
            "test", "range.bam"), 0);
    auto ar = bam.allRecords;
    assert(ar.empty == false);
    auto read = ar.front;
    read["MD"] = "2G11^GATC7T6^A11";
    assert(MDItr(read).array == "==G===========GATC=======T======A===========");
}
