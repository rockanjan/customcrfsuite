/*
 *      CRF1d feature generator (dyad features).
 *
 * Copyright (c) 2007-2010, Naoaki Okazaki
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the names of the authors nor the names of its contributors
 *       may be used to endorse or promote products derived from this
 *       software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/* $Id$ */


#ifdef    HAVE_CONFIG_H
#include <config.h>
#endif/*HAVE_CONFIG_H*/

#include <os.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <crfsuite.h>

#include "logging.h"
#include "crf1d.h"
#include "rumavl.h"    /* AVL tree library necessary for feature generation. */

/**
 * Feature set.
 */
typedef struct {
    RUMAVL* avl;    /**< Root node of the AVL tree. */
    int num;        /**< Number of features in the AVL tree. */
} featureset_t;


#define    COMP(a, b)    ((a)>(b))-((a)<(b))

int ANJAN_NUM_LABEL = -1; 
long* label_counts;
float* chisquare_table_005;

void display_observed_table(const long** observed){
    //display the table
    fprintf(stdout, "OBSERVED TABLE:\n");
    for(int x=0; x< ANJAN_NUM_LABEL+1; x++){
        fprintf(stdout, "%d\t%d\t%d\n", observed[x][0], observed[x][1], observed[x][2]);
    }
    fprintf(stdout, "\n");
}
void display_expected_table(const float** expected){
    //display the table
    fprintf(stdout, "EXPECTED TABLE:\n");
    for(int x=0; x< ANJAN_NUM_LABEL; x++){
        fprintf(stdout, "%f\t%f\n", expected[x][0], expected[x][1]);
    }
    fprintf(stdout, "\n");
}

static float compute_chisquared(long** observed, int DEBUG_TRANSITIONS){
    float chisquared = 0;
    //calculate totals for expected
    long sum_col1 = 0;
    long sum_col2 = 0;
    for(int x=0; x< ANJAN_NUM_LABEL; x++){
        observed[x][2] = observed[x][0] + observed[x][1];
        sum_col1 += observed[x][0];
        sum_col2 += observed[x][1];
    }
    observed[ANJAN_NUM_LABEL][0]=sum_col1;
    observed[ANJAN_NUM_LABEL][1]=sum_col2;
    observed[ANJAN_NUM_LABEL][2]=sum_col1 + sum_col2;

    //Expected
    float** expected = malloc( ANJAN_NUM_LABEL* sizeof(float*));
    if(expected != NULL){
        for(int x=0; x<ANJAN_NUM_LABEL; x++){
            expected[x] = (float*) malloc(2 * sizeof(float));
        }
    } else {
        fprintf(stderr, "Error in creating expected counts");
        exit -1;
    }
    
    //compute table
    for(int x=0; x<ANJAN_NUM_LABEL; x++){
        expected[x][0] = (float) 1.0 * observed[x][2] * observed[ANJAN_NUM_LABEL][0] / observed[ANJAN_NUM_LABEL][2];
        expected[x][1] = (float) 1.0 * observed[x][2] * observed[ANJAN_NUM_LABEL][1] / observed[ANJAN_NUM_LABEL][2];
    }
    
    if(DEBUG_TRANSITIONS)
        display_expected_table(expected);
    
    //calculate value
    for(int x=0; x<ANJAN_NUM_LABEL; x++){
        for(int y=0; y<2; y++){
            chisquared += (observed[x][y] - expected[x][y]) * (observed[x][y] - expected[x][y]) / expected[x][y];
        }
    }    
    //free expected
    for(int x=0; x<ANJAN_NUM_LABEL; x++){
        free(expected[x]);
    }
    free(expected);
    return chisquared;
}

static int featureset_comp(const void *x, const void *y, size_t n, void *udata)
{
    int ret = 0;
    const crf1df_feature_t* f1 = (const crf1df_feature_t*)x;
    const crf1df_feature_t* f2 = (const crf1df_feature_t*)y;

    ret = COMP(f1->type, f2->type);
    if (ret == 0) {
        ret = COMP(f1->src, f2->src);
        if (ret == 0) {
            ret = COMP(f1->dst, f2->dst);
        }
    }
    return ret;
}

static featureset_t* featureset_new()
{
    featureset_t* set = NULL;
    set = (featureset_t*)calloc(1, sizeof(featureset_t));
    if (set != NULL) {
        set->num = 0;
        set->avl = rumavl_new(
            sizeof(crf1df_feature_t), featureset_comp, NULL, NULL);
        if (set->avl == NULL) {
            free(set);
            set = NULL;
        }
    }
    return set;
}

static void featureset_delete(featureset_t* set)
{
    if (set != NULL) {
        rumavl_destroy(set->avl);
        free(set);
    }
}

static int featureset_add(featureset_t* set, const crf1df_feature_t* f)
{
    /* Check whether if the feature already exists. */
    crf1df_feature_t *p = (crf1df_feature_t*)rumavl_find(set->avl, f);
    if (p == NULL) {
        /* Insert the feature to the feature set. */
        rumavl_insert(set->avl, f);
        ++set->num;
    } else {
        /* An existing feature: add the observation expectation. */
        p->freq += f->freq;
    }
    return 0;
}

//converts feature set (in AVL) into feature array
static crf1df_feature_t*
featureset_generate(
    int *ptr_num_features,
    featureset_t* set,
    floatval_t minfreq
    )
{
    chisquare_table_005 = (float*) malloc(11 * sizeof(float)); //for upto 20 degree of freedoms, 0 ignored
    chisquare_table_005[1] = 3.84;
    chisquare_table_005[2] = 5.99;
    chisquare_table_005[3] = 7.82;
    chisquare_table_005[4] = 9.49;
    chisquare_table_005[5] = 11.07;
    chisquare_table_005[6] = 12.59;
    chisquare_table_005[7] = 14.07;
    chisquare_table_005[8] = 15.51;
    chisquare_table_005[9] = 16.92;
    chisquare_table_005[10] = 18.31;
    
    float chisquare_table_001_dof2 = 9.21;
    float chisquare_table_0001_dof2 = 13.82;
    
    float chisquare_table_005_value = chisquare_table_005[ANJAN_NUM_LABEL - 1];
    //free chisquared table
    free(chisquare_table_005);
    
    int n = 0, k = 0;
    RUMAVL_NODE *node = NULL;
    crf1df_feature_t *f = NULL;
    crf1df_feature_t *features = NULL;
    
    fprintf(stdout, "LABEL COUNTS: \n");
    for(int x=0; x<ANJAN_NUM_LABEL; x++){
        fprintf(stdout, "%d = %ld\n", x, label_counts[x]);
    }
    node = NULL;
    
    while ((node = rumavl_node_next(set->avl, node, 1, (void**)&f)) != NULL) {
        //observed table
        long** observed = malloc( (ANJAN_NUM_LABEL+1)* sizeof(long*));
        if(observed){
            for(int x=0; x<ANJAN_NUM_LABEL+1; x++){
                observed[x] = malloc(3 * sizeof(long));
            }
        }
        for(int x=0; x<ANJAN_NUM_LABEL+1; x++){
            observed[x][0] = 0; observed[x][1] = 0; observed[x][2] = 0;
        }
        
                
        for(int x=0; x<ANJAN_NUM_LABEL; x++){
            if(f->dst == x){
                observed[x][1] = f->freq;
            } else {
                observed[x][0] = label_counts[x];
            }
        }
        observed[f->dst][0] = label_counts[f->dst] - observed[f->dst][1]; //remaining
        int DEBUG_TRANSITIONS = 0;
        if(f->type == 1){
            DEBUG_TRANSITIONS = 1;
        }
        float chisquared_observed_value = compute_chisquared(observed, DEBUG_TRANSITIONS);
        if(f->type == 1)
            display_observed_table(observed);
        f->chisquared_value = chisquared_observed_value;
        //free
        for(int x=0; x<ANJAN_NUM_LABEL+1; x++){
            free(observed[x]);
        }
        free(observed);
    }    
    /*Anjan code complete*/
    
    //number of features without selection
    long number_of_features_without_selection = 0;
    
    /* The first pass: count the number of valid features. */
    while ((node = rumavl_node_next(set->avl, node, 1, (void**)&f)) != NULL) {
        if(f->type == 1){
            //transition feature
            fprintf(stdout, "Transition chisquared_obs = %f, table = %f\n", f->chisquared_value, chisquare_table_005_value);
        }
        if (minfreq <= f->freq &&  f->chisquared_value > chisquare_table_005_value ) {
            ++n;
        }
        ++number_of_features_without_selection;
    }
    fprintf(stdout, "Feature before selection: %ld \t after selection: %ld\n", number_of_features_without_selection, n);

    /* The second path: copy the valid features to the feature array. */
    features = (crf1df_feature_t*)calloc(n, sizeof(crf1df_feature_t));
    if (features != NULL) {
        node = NULL;
        while ((node = rumavl_node_next(set->avl, node, 1, (void**)&f)) != NULL) {
            //fprintf(stdout, "Feature: type:%d, %d --> %d ____ freq:%f", f->type, f->src, f->dst, f->freq);
            //fprintf(stdout, "\n");
            if (minfreq <= f->freq &&  f->chisquared_value > chisquare_table_005_value) {
                memcpy(&features[k], f, sizeof(crf1df_feature_t));
                ++k;
            }
        }
        *ptr_num_features = n;
        return features;
    } else {
        *ptr_num_features = 0;
        return NULL;
    }
}


/* creates all possible features and stores the counts of them
 */
//state feature and transition feature generation
crf1df_feature_t* crf1df_generate(
    int *ptr_num_features,
    dataset_t *ds,
    int num_labels,
    int num_attributes,
    int connect_all_attrs,
    int connect_all_edges,
    floatval_t minfreq,
    crfsuite_logging_callback func,
    void *instance
    )
{
    ANJAN_NUM_LABEL = num_labels;
    int c, i, j, s, t;
    crf1df_feature_t f;
    crf1df_feature_t *features = NULL;
    featureset_t* set = NULL;
    const int N = ds->num_instances;
    const int L = num_labels;
    logging_t lg;

    lg.func = func;
    lg.instance = instance;
    lg.percent = 0;

    /* Create an instance of feature set. */
    set = featureset_new();

    /* Loop over the sequences in the training data. */
    logging_progress_start(&lg);
    
    /*Anjan code*/
    label_counts = (long*) calloc(ANJAN_NUM_LABEL, sizeof(long)); 
    
    for (s = 0;s < N;++s) { //each sentence (or sequence)
        int prev = L, cur = 0;
        const crfsuite_item_t* item = NULL;
        const crfsuite_instance_t* seq = dataset_get(ds, s);
        const int T = seq->num_items;
        /* Loop over the items in the sequence. */
        for (t = 0;t < T;++t) { //for each item in a line
            item = &seq->items[t];
            cur = seq->labels[t];
            
            /* Transition feature: label #prev -> label #(item->yid).
               Features with previous label #L are transition BOS. */
            
            /*****Modify here if you want to include BOS to other label transition *****/
            if (prev != L) {
                f.type = FT_TRANS;
                f.src = prev;
                f.dst = cur;
                f.freq = 1;
                featureset_add(set, &f);
            } 
           
            //state feature generation
            for (c = 0;c < item->num_contents;++c) {
                /* State feature: attribute #a -> state #(item->yid). */
                f.type = FT_STATE;
                f.src = item->contents[c].aid;
                f.dst = cur;
                f.freq = item->contents[c].value;
                featureset_add(set, &f);

                /* Generate state features connecting attributes with all
                   output labels. These features are not unobserved in the
                   training data (zero expexcations). */
                if (connect_all_attrs) {
                    for (i = 0;i < L;++i) {
                        f.type = FT_STATE;
                        f.src = item->contents[c].aid;
                        f.dst = i;
                        f.freq = 0;
                        featureset_add(set, &f);
                    }
                }
            }
            label_counts[cur]++;
            prev = cur;
        }
        
        logging_progress(&lg, s * 100 / N);
    }
    logging_progress_end(&lg);

    /* Generate edge features representing all pairs of labels.
       These features are not unobserved in the training data
       (zero expexcations). */
    if (connect_all_edges) {
        for (i = 0;i < L;++i) {
            for (j = 0;j < L;++j) {
                f.type = FT_TRANS;
                f.src = i;
                f.dst = j;
                f.freq = 0;
                featureset_add(set, &f);
            }
        }
    }

    /* Convert the feature set to an feature array. */
    features = featureset_generate(ptr_num_features, set, minfreq);//generates features

    /* Delete the feature set. */
    featureset_delete(set);

    return features;
}

int crf1df_init_references(
    feature_refs_t **ptr_attributes,
    feature_refs_t **ptr_trans,
    const crf1df_feature_t *features,
    const int K,
    const int A,
    const int L
    )
{
    int i, k;
    feature_refs_t *fl = NULL;
    feature_refs_t *attributes = NULL;
    feature_refs_t *trans = NULL;

    /*
        The purpose of this routine is to collect references (indices) of:
        - state features fired by each attribute (attributes)
        - transition features pointing from each label (trans)
    */

    /* Allocate arrays for feature references. */
    attributes = (feature_refs_t*)calloc(A, sizeof(feature_refs_t));
    if (attributes == NULL) goto error_exit;
    trans = (feature_refs_t*)calloc(L, sizeof(feature_refs_t));
    if (trans == NULL) goto error_exit;

    /*
        Firstly, loop over the features to count the number of references.
        We don't use realloc() to avoid memory fragmentation.
     */
    for (k = 0;k < K;++k) {
        const crf1df_feature_t *f = &features[k];
        switch (f->type) {
        case FT_STATE:
            attributes[f->src].num_features++;
            break;
        case FT_TRANS:
            trans[f->src].num_features++;
            break;
        }
    }

    /*
        Secondarily, allocate memory blocks to store the feature references.
        We also clear fl->num_features fields, which will be used as indices
        in the next phase.
     */
    for (i = 0;i < A;++i) {
        fl = &attributes[i];
        fl->fids = (int*)calloc(fl->num_features, sizeof(int));
        if (fl->fids == NULL) goto error_exit;
        fl->num_features = 0;
    }
    for (i = 0;i < L;++i) {
        fl = &trans[i];
        fl->fids = (int*)calloc(fl->num_features, sizeof(int));
        if (fl->fids == NULL) goto error_exit;
        fl->num_features = 0;
    }

    /*
        Finally, store the feature indices.
     */
    for (k = 0;k < K;++k) {
        const crf1df_feature_t *f = &features[k];
        switch (f->type) {
        case FT_STATE:
            fl = &attributes[f->src];
            fl->fids[fl->num_features++] = k;
            break;
        case FT_TRANS:
            fl = &trans[f->src];
            fl->fids[fl->num_features++] = k;
            break;
        }
    }

    *ptr_attributes = attributes;
    *ptr_trans = trans;
    return 0;

error_exit:
    if (attributes != NULL) {
        for (i = 0;i < A;++i) free(attributes[i].fids);
        free(attributes);
    }
    if (trans != NULL) {
        for (i = 0;i < L;++i) free(trans[i].fids);
        free(trans);
    }
    *ptr_attributes = NULL;
    *ptr_trans = NULL;
    return -1;
}
