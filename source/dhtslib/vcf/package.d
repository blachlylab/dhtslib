module dhtslib.vcf;

public import dhtslib.vcf.record;
public import dhtslib.vcf.reader;
public import dhtslib.vcf.header;
public import dhtslib.vcf.writer;

import std.meta : AliasSeq;
import htslib.vcf;

/// Replacement for htslib BCF_HL_*
enum HeaderRecordType
{
    NULL = -1,
    FILTER =    BCF_HL_FLT, /// header line: FILTER
    INFO =      BCF_HL_INFO, /// header line: INFO
    FORMAT =    BCF_HL_FMT, /// header line: FORMAT
    CONTIG =    BCF_HL_CTG, /// header line: contig
    STRUCT =    BCF_HL_STR, /// header line: structured header line TAG=<A=..,B=..>
    GENERIC =   BCF_HL_GEN, /// header line: generic header line
}

/// Strings for HeaderRecordType
static immutable HeaderRecordTypeStrings = ["FILTER","INFO","FORMAT","contig", "struct", "generic"];

/// Replacement for htslib BCF_HT_*
enum HeaderTypes
{
    NULL = -1,
    FLAG =  BCF_HT_FLAG, /// header type: FLAG// header type
    INT =   BCF_HT_INT, /// header type: INTEGER
    REAL =  BCF_HT_REAL, /// header type: REAL,
    FLOAT =  BCF_HT_REAL, /// copy of REAL
    STR =   BCF_HT_STR, /// header type: STRING
    CHAR = BCF_HT_STR,
    LONG =  BCF_HT_LONG, // BCF_HT_INT, but for int64_t values; VCF only!
}

/// Strings for HeaderTypes
static immutable HeaderTypesStrings = ["Flag","Integer","Float","String","Character","Long"];

/// Replacement for htslib BCF_VL_*
enum HeaderLengths
{
    NULL = -1,
    FIXED = BCF_VL_FIXED, /// variable length: fixed length
    VAR =   BCF_VL_VAR, /// variable length: variable
    A =     BCF_VL_A, /// variable length: one field per alt allele
    G =     BCF_VL_G, /// variable length: one field per genotype
    R =     BCF_VL_R, /// variable length: one field per allele including ref
}

/// Strings for HDR_LENGTH
static immutable  HeaderLengthsStrings = ["FIXED",".","A","G","R"];

/// Replacement for htslib BCF_DT_*
enum HeaderDictTypes
{
    ID =     BCF_DT_ID, /// dictionary type: ID
    CTG =    BCF_DT_CTG, /// dictionary type: CONTIG
    SAMPLE = BCF_DT_SAMPLE, /// dictionary type: SAMPLE
}

/// Replacement for htslib BCF_BT_*
enum RecordType
{
    NULL =   0,  /// null
    INT8 =   BCF_BT_INT8,  /// int8
    INT16 =  BCF_BT_INT16,  /// int16
    INT32 =  BCF_BT_INT32,  /// int32
    INT64 =  BCF_BT_INT64,  /// Unofficial, for internal use only per htslib headers 
    FLOAT =  BCF_BT_FLOAT,  /// float (32?)
    CHAR =   BCF_BT_CHAR  /// char (8 bit)
}

/// Byte sizes for RecordType
static immutable RecordTypeSizes = [0, 1, 2, 4, 8, 4, 1];
alias RecordTypeToDType = AliasSeq!(null, byte, short, int, long, float, null, string);

/// Replacement for htslib VCF_*
enum VariantType
{
    REF =       VCF_REF,  /// ref (e.g. in a gVCF)
    SNP =       VCF_SNP,  /// SNP 
    MNP =       VCF_MNP,  /// MNP
    INDEL =     VCF_INDEL,  /// INDEL
    OTHER =     VCF_OTHER,  /// other (e.g. SV)
    BND =       VCF_BND, /// breakend
    OVERLAP =   VCF_OVERLAP, /// overlapping deletion, ALT=* 
}

/// Replacement for htslib BCF_UN_*
enum UNPACK
{
    ALT =   BCF_UN_STR, // up to ALT inclusive
    FILT =   BCF_UN_FLT, // up to FILTER
    INFO =  BCF_UN_INFO, // up to INFO
    SHR =   BCF_UN_SHR, // all shared information
    FMT =   BCF_UN_FMT, // unpack format and each sample
    IND =   BCF_UN_IND, // a synonym of UNPACK.FMT
    ALL =   BCF_UN_ALL, // everything
}
