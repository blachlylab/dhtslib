module dhtslib.vcf;

public import dhtslib.vcf.record;
public import dhtslib.vcf.reader;
public import dhtslib.vcf.header;
public import dhtslib.vcf.writer;

import std.meta : AliasSeq;
import htslib.vcf;

/// Represents the classification of a headerline
///
/// ##INFO=<...>
///   ====
/// 
/// Replacement for htslib BCF_HL_*
enum HeaderRecordType
{
    None = -1,
    Filter =    BCF_HL_FLT, /// header line: FILTER
    Info =      BCF_HL_INFO,/// header line: INFO
    Format =    BCF_HL_FMT, /// header line: FORMAT
    Contig =    BCF_HL_CTG, /// header line: contig
    Struct =    BCF_HL_STR, /// header line: structured header line TAG=<A=..,B=..>
    Generic =   BCF_HL_GEN, /// header line: generic header line
}

/// Strings for HeaderRecordType
static immutable HeaderRecordTypeStrings = ["FILTER","INFO","FORMAT","contig", "struct", "generic"];

/// Represents the classification of a headerline
///
/// ##INFO=<Number=A, Type=Integer>
///                       =======
/// 
/// Replacement for htslib BCF_HT_*
enum HeaderTypes
{
    None =      -1,
    Flag =      BCF_HT_FLAG, /// header type: FLAG// header type
    Integer =   BCF_HT_INT, /// header type: INTEGER
    Float =      BCF_HT_REAL, /// header type: REAL,
    String =    BCF_HT_STR, /// header type: STRING
    Character = 4,
    Long =      BCF_HT_LONG, // BCF_HT_INT, but for int64_t values; VCF only!
}

/// Strings for HeaderTypes
/// 
/// doesn't work as needs compile time
/// enum HeaderTypesStrings = __traits(allMembers, HeaderTypes);
///
/// works but includes "None" which is throwing off indexing
/// enum HeaderTypesStrings = [__traits(allMembers, HeaderTypes)];
enum HeaderTypesStrings = [__traits(allMembers, HeaderTypes)][1..$];


/// Represents the classification of a headerline
///
/// ##INFO=<Number=A, Type=Integer>
///               =
///
/// if FIXED
/// ##INFO=<Number=2, Type=Integer>
///               =
/// 
/// Replacement for htslib BCF_VL_*
enum HeaderLengths
{
    None =              -1,
    Fixed =             BCF_VL_FIXED, /// variable length: fixed length
    Variable =          BCF_VL_VAR, /// variable length: variable
    OnePerAltAllele =   BCF_VL_A, /// variable length: one field per alt allele
    OnePerGenotype =    BCF_VL_G, /// variable length: one field per genotype
    OnePerAllele =      BCF_VL_R, /// variable length: one field per allele including ref
}

/// Strings for HDR_LENGTH
static immutable  HeaderLengthsStrings = ["FIXED",".","A","G","R"];

/// Used to index into bcf_hdr_t's id field of type bcf_idpair_t*[3]
///
/// i.e as used from VCFRecord where this.vcfheader is a VCFHeader:
/// this.vcfheader.hdr.id[HeaderDictTypes.Id]
///
/// Replacement for htslib BCF_DT_*
enum HeaderDictTypes
{
    Id =        BCF_DT_ID, /// dictionary type: ID
    Contig =    BCF_DT_CTG, /// dictionary type: Contig
    Sample =    BCF_DT_SAMPLE, /// dictionary type: SAMPLE
}

/// Used by InfoField (bcf_info_t) and FormatField (bcf_fmt_t) 
/// to identify the underlying htslib/bcf1_t info and format data
/// type and size. This data is stored in ubyte arrays.
///
/// 
/// Replacement for htslib BCF_BT_*
enum RecordType
{
    Null =   0,  /// null
    Int8 =   BCF_BT_INT8,  /// int8
    Int16 =  BCF_BT_INT16,  /// int16
    Int32 =  BCF_BT_INT32,  /// int32
    Int64 =  BCF_BT_INT64,  /// Unofficial, for internal use only per htslib headers 
    Float =  BCF_BT_FLOAT,  /// float (32?)
    Char =   BCF_BT_CHAR  /// char (8 bit)
}

/// Byte sizes for RecordType
static immutable RecordTypeSizes = [0, byte.sizeof, short.sizeof, int.sizeof, long.sizeof, float.sizeof, 0, char.sizeof];
alias RecordTypeToDType = AliasSeq!(null, byte, short, int, long, float, null, string);

/// Replacement for htslib VCF_*
enum VariantType
{
    Ref =           VCF_REF,  /// ref (e.g. in a gVCF)
    Snp =           VCF_SNP,  /// SNP 
    Mnp =           VCF_MNP,  /// MNP
    Indel =         VCF_INDEL,  /// INDEL
    Other =         VCF_OTHER,  /// other (e.g. SV)
    Breakend =      VCF_BND, /// breakend
    Overlap =       VCF_OVERLAP, /// overlapping deletion, ALT=* 
}

/// Levels identifiers for unpacking the underlying variable length
/// data in the bcf1_t. Values are inclusive 
/// i.e UnpackLevels.AltAllele unpacks all data before and including the ALT allele
/// Replacement for htslib BCF_UN_*
enum UnpackLevels
{
    None =              0,
    AltAllele =         BCF_UN_STR, // up to ALT inclusive
    Filter =            BCF_UN_FLT, // up to Filter
    Info =              BCF_UN_INFO, // up to Info
    SharedFields =      BCF_UN_SHR, // all shared information
    Format =            BCF_UN_FMT, // unpack format and each sample
    IndividualFields =  BCF_UN_IND, // a synonym of UNPACK.FMT
    All =               BCF_UN_ALL, // everything
}
