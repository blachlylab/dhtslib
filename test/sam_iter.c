#include <stdio.h>

#include "dhtslib/htslib/include/hts.h"
#include "dhtslib/htslib/include/hts_log.h"
#include "dhtslib/htslib/include/sam.h"

int main(void)
{
    printf("Query c raw htslib\n");
    int j = 0;
    char  *fn = "/Users/james/Documents/Development/blachlylab/funmap/ENCFF399AWI.bam";

    samFile *fp = sam_open(fn, "r");
    printf("fp      %x\n", (int)fp);
    printf("fp->fn  %s\n", fp->fn);
    printf("fp->fp  %x\n", (int)fp->fp.bgzf);
    printf("fp->format.category %d\n", fp->format.category);
    printf("fp->format.compression %d\n", fp->format.compression);

    hts_idx_t *idx = sam_index_load(fp, fn);
    bam1_t *b = bam_init1();
    hts_itr_t *iter;
    int r = 0;
    if ((iter = sam_itr_queryi(idx, 0, 1000000, 2000000)) == 0) {
        hts_log_error("Failed to parse region");
        return 1;
    }
    printf("iter == %x\n", (int)iter);
    printf("fp->fp  %x\n", (int)fp->fp.bgzf);
    
    while ((r = sam_itr_next(fp, iter, b)) >= 0) {
        j++;
    }

    printf("Processed %d records with raw iter\n", j);

    hts_itr_destroy(iter);
    bam_destroy1(b);
    hts_close(fp);

    return 0;
}