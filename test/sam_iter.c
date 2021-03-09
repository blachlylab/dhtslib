#include <stdio.h>

#include "htslib/htslib/hts.h"
#include "htslib/htslib/hts_log.h"
#include "htslib/htslib/sam.h"

int main(void)
{
    printf("Query c raw htslib\n");
    int j = 0;
    char  *fn = "../htslib/test/range.bam";

    samFile *fp = sam_open(fn, "r");
    sam_hdr_t* h = sam_hdr_read(fp);

    printf("fp      %x\n", (int)fp);
    printf("fp->fn  %s\n", fp->fn);
    printf("fp->fp  %x\n", (int)fp->fp.bgzf);
    printf("fp->format.category %d\n", fp->format.category);
    printf("fp->format.compression %d\n", fp->format.compression);

    hts_idx_t *idx = sam_index_load(fp, fn);
    bam1_t *b = bam_init1();
    hts_itr_t *iter;
    int r = 0;
    if ((iter = sam_itr_queryi(idx, 0, 1000, 2000)) == 0) {
        hts_log_error("Failed to parse region");
        return 1;
    }
    printf("iter == %x\n", (int)iter);
    printf("fp->fp  %x\n", (int)fp->fp.bgzf);
    
    while ((r = sam_itr_next(fp, iter, b)) >= 0) {
        j++;
    }

    printf("Processed %d records with raw iter\n", j);

    printf("Now sam_read1\n");
    sam_read1(fp, h, b);
    printf("%s\n", bam_get_qname(b));

    hts_itr_destroy(iter);
    bam_destroy1(b);
    sam_hdr_destroy(h);
    hts_close(fp);

    return 0;
}
