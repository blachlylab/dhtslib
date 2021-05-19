/// @file htslib/vcfutils.h
/// Allele-related utility functions.
/*
    Copyright (C) 2012, 2013, 2015-2016 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

extern (C):

struct kbitset_t;

/**
 *  bcf_trim_alleles() - remove ALT alleles unused in genotype fields
 *  @header:  for access to BCF_DT_ID dictionary
 *  @line:    VCF line obtain from vcf_parse1
 *
 *  Returns the number of removed alleles on success or negative
 *  on error:
 *      -1 .. some allele index is out of bounds
 *      -2 .. could not remove alleles
 */
int bcf_trim_alleles (const(bcf_hdr_t)* header, bcf1_t* line);

/**
 *  bcf_remove_alleles() - remove ALT alleles according to bitmask @mask
 *  @header:  for access to BCF_DT_ID dictionary
 *  @line:    VCF line obtained from vcf_parse1
 *  @mask:    alleles to remove
 *
 *  If you have more than 31 alleles, then the integer bit mask will
 *  overflow, so use bcf_remove_allele_set instead
 *  Returns 0 on success, <0 on error
 */
int bcf_remove_alleles (const(bcf_hdr_t)* header, bcf1_t* line, int mask);

/**
 *  bcf_remove_allele_set() - remove ALT alleles according to bitset @rm_set
 *  @header:  for access to BCF_DT_ID dictionary
 *  @line:    VCF line obtained from vcf_parse1
 *  @rm_set:  pointer to kbitset_t object with bits set for allele
 *            indexes to remove
 *
 *  Returns 0 on success or -1 on failure
 *
 *  Number=A,R,G INFO and FORMAT fields will be updated accordingly.
 */
int bcf_remove_allele_set (
    const(bcf_hdr_t)* header,
    bcf1_t* line,
    const(kbitset_t)* rm_set);

/**
 *  bcf_calc_ac() - calculate the number of REF and ALT alleles
 *  @header:  for access to BCF_DT_ID dictionary
 *  @line:    VCF line obtained from vcf_parse1
 *  @ac:      array of length line->n_allele
 *  @which:   determine if INFO/AN,AC and indv fields be used
 *
 *  Returns 1 if the call succeeded, or 0 if the value could not
 *  be determined.
 *
 *  The value of @which determines if existing INFO/AC,AN can be
 *  used (BCF_UN_INFO) and and if indv fields can be split (BCF_UN_FMT).
 */
int bcf_calc_ac (const(bcf_hdr_t)* header, bcf1_t* line, int* ac, int which);

/**
 * bcf_gt_type() - determines type of the genotype
 * @fmt_ptr:  the GT format field as set for example by set_fmt_ptr
 * @isample:  sample index (starting from 0)
 * @ial:      index of the 1st non-reference allele (starting from 1)
 * @jal:      index of the 2nd non-reference allele (starting from 1)
 *
 * Returns the type of the genotype (one of GT_HOM_RR, GT_HET_RA,
 * GT_HOM_AA, GT_HET_AA, GT_HAPL_R, GT_HAPL_A or GT_UNKN). If $ial
 * is not NULL and the genotype has one or more non-reference
 * alleles, $ial will be set. In case of GT_HET_AA, $ial is the
 * position of the allele which appeared first in ALT. If $jal is
 * not null and the genotype is GT_HET_AA, $jal will be set and is
 * the position of the second allele in ALT.
 */
enum GT_HOM_RR = 0; // note: the actual value of GT_* matters, used in dosage r2 calculation
enum GT_HOM_AA = 1;
enum GT_HET_RA = 2;
enum GT_HET_AA = 3;
enum GT_HAPL_R = 4;
enum GT_HAPL_A = 5;
enum GT_UNKN = 6;
int bcf_gt_type (bcf_fmt_t* fmt_ptr, int isample, int* ial, int* jal);

int bcf_acgt2int (char c);

extern (D) auto bcf_int2acgt(T)(auto ref T i)
{
    return "ACGT"[i];
}

/**
  * bcf_ij2G() - common task: allele indexes to Number=G index (diploid)
  * @i,j:  allele indexes, 0-based, i<=j
  *
  * Returns index to the Number=G diploid array
  */
extern (D) auto bcf_ij2G(T0, T1)(auto ref T0 i, auto ref T1 j)
{
    return j * (j + 1) / 2 + i;
}
