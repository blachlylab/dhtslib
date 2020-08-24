module dhtslib.md;

import dhtslib.sam: SAMRecord;
import dhtslib.cigar;
import dhtslib.htslib.hts_log;
import std.regex;
import std.traits: ReturnType;
import std.conv: to;

/// regex to extract MD string groups
/// ex: "11A3^G" -> [(11, "A"), (3, "^G")] 
auto MDREGEX = regex(`(\d+)([\^ATGCN]*)`);

struct MDPair 
{
    int numMatches;
    string __mismatch;

    string mismatch() const {
        if(isDel)
            return __mismatch[1..$];
        else
            return __mismatch;
    }
    bool isDel() const { 
        return !isLast && __mismatch[0]=='^';
    }
    
    bool isLast() const { 
        return __mismatch.length == 0; 
    }

    invariant 
    {
        if(__mismatch.length > 1)
            assert(__mismatch[0] == '^');
        if(__mismatch.length == 1)
            assert(__mismatch[0] != '^');
    }
    string toString(){
        return numMatches.to!string~__mismatch;
    }
}

auto getMDPairs(SAMRecord rec){
    struct MDPairs
    {
        string md_string;
        ReturnType!generateMDmatches matches;

        this(string md){
            md_string = md;
            matches = generateMDmatches();
        }
        
        auto generateMDmatches(){
            return matchAll(md_string, MDREGEX);
        }

        MDPair front(){
            return MDPair(matches.front[1].to!int,matches.front[2]);
        }

        bool empty(){
            return matches.empty;
        }

        void popFront(){
            matches.popFront;
        }

    }
    auto tag = rec["MD"];
    debug if(!tag.exists) hts_log_error(__FUNCTION__,"MD tag not present");
    return MDPairs(tag.toString);
}

// debug(dhtslib_unittest) 
unittest
{
    import std.stdio;
    import dhtslib.sam;
    import std.path:buildPath,dirName;
    auto bam = SAMFile(buildPath(dirName(dirName(dirName(__FILE__))),"htslib","test","range.bam"), 0);
    auto read=bam.all_records.front;
    read["MD"] = "2G11^GATC7T6^A11";
    writeln(getMDPairs(read));
}

struct MDItr {
    ReturnType!getMDPairs mdPairs;
    int currentlen;
    string currentSeq;
    bool isLast;

    this(SAMRecord rec){
        mdPairs = getMDPairs(rec);
        currentlen = mdPairs.front.numMatches;
        currentSeq = mdPairs.front.mismatch;
        isLast = mdPairs.front.isLast;
        mdPairs.popFront;
    }

    char front(){
        if(currentlen) return '=';
        else return currentSeq[0];
    }

    void popFront(){
        if(currentlen) currentlen--;
        else if(currentSeq.length > 1) currentSeq = currentSeq[1..$];
        else
        {
            currentlen = mdPairs.front.numMatches;
            currentSeq = mdPairs.front.mismatch;
            isLast = mdPairs.front.isLast;
            mdPairs.popFront;
        }
    }

    bool empty(){
        return isLast && currentlen == 0;
    }

}

unittest
{
    import std.stdio;
    import dhtslib.sam;
    import std.path:buildPath,dirName;
    auto bam = SAMFile(buildPath(dirName(dirName(dirName(__FILE__))),"htslib","test","range.bam"), 0);
    auto read=bam.all_records.front;
    read["MD"] = "2G11^GATC7T6^A11";
    auto mdPairs = MDItr(read);
    writeln(mdPairs);
}