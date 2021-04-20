/**

SAMRecord and SAMFile are wrappers for htslib functions relating to SAM/BAM/CRAM* files

SAMRecord is a structured representation of a SAM/BAM/CRAM* record,
backed internally by htslib's bam1_t, but with convenient getters and setters
for record attributes, including functions not included in vanilla htslib 
like returning the sequence or the qscore string (NB: actually char*)

SAMFile is a structured representation of SAM/BAM/CRAM* file,
backed internally by htslib's htsFile and bam_hdr_t,
but with convenient getters and setters for file and header attributes,
as well as query functions accessible explicitly (`query("chr1:999-9999"`)
and by indexing (`samfile["chr1", 999 .. 9999]`).
The file object can be iterated as an InputRange to obtain every record in the file.

Authors: James S Blachly, MD <james.blachly@gmail.com> ; Thomas Gregory <charles.gregory@osumc.edu>

Bugs: SAMRecord and SAMFile function only as readers, rather than writers (i.e, cannot build SAMFile)
Bugs: (*CRAM functionality is limited and untested)

Date: 2020-09-12

License: Apache 2.0

Standards: Sequence Alignment/Map Format Specification v1 14 Dec 2018 http://samtools.github.io/hts-specs/

*/
module dhtslib.sam;

public import dhtslib.sam.header;
public import dhtslib.sam.record;
public import dhtslib.sam.cigar;
public import dhtslib.sam.reader;
public import dhtslib.sam.writer;

import htslib.sam;
import htslib.hts : hts_reglist_t, hts_pair32_t;


/// Used in sorting
bool cmpInterval(hts_pair32_t a, hts_pair32_t b)
{
    if (a.beg < b.beg)
    {
        return true;
    }
    if (a.end < b.end)
    {
        return true;
    }
    return false;
}

/// Used in sorting
bool cmpRegList(hts_reglist_t a, hts_reglist_t b)
{
    if (a.tid < b.tid)
    {
        return true;
    }
    return false;
}

/// Parse text line of SAM; Used in unittest
private int parseSam(string line, bam_hdr_t* header, bam1_t* b)
{
    import htslib.kstring : kstring_t;
    import std.utf : toUTFz;

    kstring_t k;
    k.s = toUTFz!(char*)(line.dup);
    k.m = line.length + 1;
    k.l = line.length + 1;
    return sam_parse1(&k, header, b);
}

/// Nucleotide complement table; from samtools/sam_view.c
private const(char)[16] seq_comp_table = [0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15];

/// Reverse a string in place; from samtools/sam_view.c
//  TODO: Could be sped up, and potentially made safer? by passing strlen since already known
private char* reverse(char* str)
{
    import core.stdc.string : strlen;

    auto i = strlen(str) - 1, j = 0;
    char ch;
    while (i > j)
    {
        ch = str[i];
        str[i] = str[j];
        str[j] = ch;
        i--;
        j++;
    }
    return str;
}

