module dhtslib.vcf;

public import dhtslib.vcf.record;
public import dhtslib.vcf.reader;
public import dhtslib.vcf.header;
public import dhtslib.vcf.writer;

import std.meta : AliasSeq;

/// Replacement for htslib BCF_HL_*
enum HeaderRecordType
{
    NULL = -1,
    FILTER =    0, /// header line: FILTER
    INFO =      1, /// header line: INFO
    FORMAT =    2, /// header line: FORMAT
    CONTIG =    3, /// header line: contig
    STRUCT =    4, /// header line: structured header line TAG=<A=..,B=..>
    GENERIC =   5, /// header line: generic header line
}

/// Strings for HDR_LINE
static immutable HeaderRecordTypeStrings = ["FILTER","INFO","FORMAT","contig", "struct", "generic"];

/// Replacement for htslib BCF_HT_*
enum HeaderTypes
{
    NULL = -1,
    FLAG =  0, /// header type: FLAG// header type
    INT =   1, /// header type: INTEGER
    REAL =  2, /// header type: REAL,
    FLOAT =  2, /// copy of REAL
    STR =   3, /// header type: STRING
    CHAR = 4,
    LONG =  INT | 0x100, // BCF_HT_INT, but for int64_t values; VCF only!
}

/// Strings for HDR_TYPE
static immutable HeaderTypesStrings = ["Flag","Integer","Float","String","Character","Long"];

/// Replacement for htslib BCF_VL_*
enum HeaderLengths
{
    NULL = -1,
    FIXED = 0, /// variable length: fixed length
    VAR =   1, /// variable length: variable
    A =     2, /// variable length: one field per alt allele
    G =     3, /// variable length: one field per genotype
    R =     4, /// variable length: one field per allele including ref
}

/// Strings for HDR_LENGTH
static immutable  HeaderLengthsStrings = ["FIXED",".","A","G","R"];

/// Replacement for htslib BCF_DT_*
enum HeaderDictTypes
{
    ID =     0, /// dictionary type: ID
    CTG =    1, /// dictionary type: CONTIG
    SAMPLE = 2, /// dictionary type: SAMPLE
}

/// Replacement for htslib BCF_BT_*
enum RecordType{
    NULL =   0,  /// null
    INT8 =   1,  /// int8
    INT16 =  2,  /// int16
    INT32 =  3,  /// int32
    INT64 =  4,  /// Unofficial, for internal use only per htslib headers 
    FLOAT =  5,  /// float (32?)
    CHAR =   7  /// char (8 bit)
}

/// Byte sizes for RecordType
static immutable RecordTypeSizes = [0, 1, 2, 4, 8, 4, 1];
alias RecordTypeToDType = AliasSeq!(null, byte, short, int, long, float, null, string);

/// Replacement for htslib VCF_*
enum VariantType
{
    REF =       0,  /// ref (e.g. in a gVCF)
    SNP =       1,  /// SNP 
    MNP =       2,  /// MNP
    INDEL =     4,  /// INDEL
    OTHER =     8,  /// other (e.g. SV)
    BND =       16, /// breakend
    OVERLAP =   32, /// overlapping deletion, ALT=* 
}

/// Replacement for htslib BCF_UN_*
enum UNPACK
{
    ALT =   1, // up to ALT inclusive
    FILT =   2, // up to FILTER
    INFO =  4, // up to INFO
    SHR =   ALT | FILT | INFO, // all shared information
    FMT =   8, // unpack format and each sample
    IND =   FMT, // a synonym of UNPACK.FMT
    ALL =   SHR | FMT, // everything
}
