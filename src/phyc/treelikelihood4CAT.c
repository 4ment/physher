// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

// Partial likelihood kernels for the empirical CAT site model, nucleotides only.
//
// CAT assigns every site pattern to a single rate category instead of averaging the
// likelihood over all of them. So the partials hold one block per node
// (tlk->cat_count == 1) while the transition matrices still hold sm->cat_count blocks
// per branch, and each pattern reads the block its category points at. That is the
// whole difference from treelikelihood4.c: the category loop is gone and the matrix
// offset is picked per pattern instead of per category.
//
// These live in their own file rather than as a branch inside the shared kernels
// because those kernels are on the hot path of every other model; threading a
// per-pattern indirection through their inner loops would put every model at risk to
// serve this one. Nothing here is reachable unless the site model is CAT.
//
// Only the lower partials need a CAT version. The transition matrices are already
// built for every category by _calculate_partials, and the root likelihood reads the
// single partials block directly (CAT sets sm->integrate to false), so
// node_log_likelihoods_4, the scaling factors and the matrix allocation are shared and
// already correct.

#include "treelikelihood4CAT.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef SSE3_ENABLED
#if defined(__aarch64__)
#include "neon2sse.h"
#else
#include <pmmintrin.h>  // SSE3
#include <xmmintrin.h>  // SSE
#endif
#endif

#include "treelikelihood.h"

// Offset of the transition matrix block holding the rate category this pattern was
// assigned to. Read straight from the array rather than through
// sm->get_site_category so the inner loop keeps no indirect call; these kernels only
// ever run for CAT, where site_category is always allocated.
static inline int _cat_offset(const SingleTreeLikelihood *tlk, int pattern) {
    return tlk->sm->site_category[pattern] * 16;
}

#pragma mark -
#pragma mark Lower Likelihood

static void _partials_states_and_states_4_cat(const SingleTreeLikelihood *tlk, int idx1,
                                              const double *matrices1, int idx2,
                                              const double *matrices2,
                                              double *partials) {
    double *pPartials = partials;

    for (int k = 0; k < tlk->pattern_count; k++) {
        int state1 = tlk->sp->patterns[idx1][k];
        int state2 = tlk->sp->patterns[idx2][k];
        int w = _cat_offset(tlk, k);

        if (state1 < 4 && state2 < 4) {
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2];
            w += 4;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2];
            w += 4;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2];
            w += 4;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2];
        } else if (state1 < 4) {
            // child 1 has a gap or unknown state so treat it as unknown
            *pPartials++ = matrices1[w + state1];
            w += 4;
            *pPartials++ = matrices1[w + state1];
            w += 4;
            *pPartials++ = matrices1[w + state1];
            w += 4;
            *pPartials++ = matrices1[w + state1];
        } else if (state2 < 4) {
            // child 2 has a gap or unknown state so treat it as unknown
            *pPartials++ = matrices2[w + state2];
            w += 4;
            *pPartials++ = matrices2[w + state2];
            w += 4;
            *pPartials++ = matrices2[w + state2];
            w += 4;
            *pPartials++ = matrices2[w + state2];
        } else {
            // both children have a gap or unknown state so set partials to 1
            memcpy(pPartials, TWENTY_DOUBLE_ONES, sizeof(double) * 4);
            pPartials += 4;
        }
    }
}

// compute e^{tQ}.partials
static void _partials_states_4_cat(const SingleTreeLikelihood *tlk, int idx1,
                                   const double *matrices1, double *partials) {
    double *pPartials = partials;

    for (int k = 0; k < tlk->pattern_count; k++) {
        int state1 = tlk->sp->patterns[idx1][k];
        int w = _cat_offset(tlk, k);

        if (state1 < 4) {
            *pPartials++ = matrices1[w + state1];
            w += 4;
            *pPartials++ = matrices1[w + state1];
            w += 4;
            *pPartials++ = matrices1[w + state1];
            w += 4;
            *pPartials++ = matrices1[w + state1];
        } else {
            // P.[1 1 1 1]^T = [1 1 1 1] when P is a probability matrix.
            // When derivatives are calculated P is NOT a probability matrix so we do
            // the calculation.
            *pPartials++ =
                matrices1[w] + matrices1[w + 1] + matrices1[w + 2] + matrices1[w + 3];
            *pPartials++ = matrices1[w + 4] + matrices1[w + 5] + matrices1[w + 6] +
                           matrices1[w + 7];
            *pPartials++ = matrices1[w + 8] + matrices1[w + 9] + matrices1[w + 10] +
                           matrices1[w + 11];
            *pPartials++ = matrices1[w + 12] + matrices1[w + 13] + matrices1[w + 14] +
                           matrices1[w + 15];
        }
    }
}

static void _partials_states_and_undefined_4_cat(const SingleTreeLikelihood *tlk,
                                                 int idx1, const double *matrices1,
                                                 const double *partials2,
                                                 const double *matrices2,
                                                 double *partials3) {
    double sum;
    int v = 0;
    double *pPartials = partials3;

    for (int k = 0; k < tlk->pattern_count; k++) {
        int state1 = tlk->sp->patterns[idx1][k];
        int w = _cat_offset(tlk, k);

        if (state1 < 4) {
            sum = matrices2[w] * partials2[v];
            sum += matrices2[w + 1] * partials2[v + 1];
            sum += matrices2[w + 2] * partials2[v + 2];
            sum += matrices2[w + 3] * partials2[v + 3];

            *pPartials++ = matrices1[w + state1] * sum;

            sum = matrices2[w + 4] * partials2[v];
            sum += matrices2[w + 5] * partials2[v + 1];
            sum += matrices2[w + 6] * partials2[v + 2];
            sum += matrices2[w + 7] * partials2[v + 3];

            *pPartials++ = matrices1[w + 4 + state1] * sum;

            sum = matrices2[w + 8] * partials2[v];
            sum += matrices2[w + 9] * partials2[v + 1];
            sum += matrices2[w + 10] * partials2[v + 2];
            sum += matrices2[w + 11] * partials2[v + 3];

            *pPartials++ = matrices1[w + 8 + state1] * sum;

            sum = matrices2[w + 12] * partials2[v];
            sum += matrices2[w + 13] * partials2[v + 1];
            sum += matrices2[w + 14] * partials2[v + 2];
            sum += matrices2[w + 15] * partials2[v + 3];

            *pPartials++ = matrices1[w + 12 + state1] * sum;
        } else {
            // Child 1 has a gap or unknown state so don't use it
            *pPartials = matrices2[w] * partials2[v];
            *pPartials += matrices2[w + 1] * partials2[v + 1];
            *pPartials += matrices2[w + 2] * partials2[v + 2];
            *pPartials += matrices2[w + 3] * partials2[v + 3];

            pPartials++;

            *pPartials = matrices2[w + 4] * partials2[v];
            *pPartials += matrices2[w + 5] * partials2[v + 1];
            *pPartials += matrices2[w + 6] * partials2[v + 2];
            *pPartials += matrices2[w + 7] * partials2[v + 3];

            pPartials++;

            *pPartials = matrices2[w + 8] * partials2[v];
            *pPartials += matrices2[w + 9] * partials2[v + 1];
            *pPartials += matrices2[w + 10] * partials2[v + 2];
            *pPartials += matrices2[w + 11] * partials2[v + 3];

            pPartials++;

            *pPartials = matrices2[w + 12] * partials2[v];
            *pPartials += matrices2[w + 13] * partials2[v + 1];
            *pPartials += matrices2[w + 14] * partials2[v + 2];
            *pPartials += matrices2[w + 15] * partials2[v + 3];

            pPartials++;
        }
        v += 4;
    }
}

static void _partials_undefined_and_undefined_4_cat(
    const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1,
    const double *partials2, const double *matrices2, double *partials3) {
    double sum1, sum2;
    int v = 0;
    double *pPartials = partials3;

    for (int k = 0; k < tlk->pattern_count; k++) {
        int w = _cat_offset(tlk, k);

        sum1 = matrices1[w] * partials1[v];
        sum2 = matrices2[w] * partials2[v];

        sum1 += matrices1[w + 1] * partials1[v + 1];
        sum2 += matrices2[w + 1] * partials2[v + 1];

        sum1 += matrices1[w + 2] * partials1[v + 2];
        sum2 += matrices2[w + 2] * partials2[v + 2];

        sum1 += matrices1[w + 3] * partials1[v + 3];
        sum2 += matrices2[w + 3] * partials2[v + 3];

        *pPartials++ = sum1 * sum2;

        sum1 = matrices1[w + 4] * partials1[v];
        sum2 = matrices2[w + 4] * partials2[v];

        sum1 += matrices1[w + 5] * partials1[v + 1];
        sum2 += matrices2[w + 5] * partials2[v + 1];

        sum1 += matrices1[w + 6] * partials1[v + 2];
        sum2 += matrices2[w + 6] * partials2[v + 2];

        sum1 += matrices1[w + 7] * partials1[v + 3];
        sum2 += matrices2[w + 7] * partials2[v + 3];

        *pPartials++ = sum1 * sum2;

        sum1 = matrices1[w + 8] * partials1[v];
        sum2 = matrices2[w + 8] * partials2[v];

        sum1 += matrices1[w + 9] * partials1[v + 1];
        sum2 += matrices2[w + 9] * partials2[v + 1];

        sum1 += matrices1[w + 10] * partials1[v + 2];
        sum2 += matrices2[w + 10] * partials2[v + 2];

        sum1 += matrices1[w + 11] * partials1[v + 3];
        sum2 += matrices2[w + 11] * partials2[v + 3];

        *pPartials++ = sum1 * sum2;

        sum1 = matrices1[w + 12] * partials1[v];
        sum2 = matrices2[w + 12] * partials2[v];

        sum1 += matrices1[w + 13] * partials1[v + 1];
        sum2 += matrices2[w + 13] * partials2[v + 1];

        sum1 += matrices1[w + 14] * partials1[v + 2];
        sum2 += matrices2[w + 14] * partials2[v + 2];

        sum1 += matrices1[w + 15] * partials1[v + 3];
        sum2 += matrices2[w + 15] * partials2[v + 3];

        *pPartials++ = sum1 * sum2;

        v += 4;
    }
}

static void _partials_undefined_4_cat(const SingleTreeLikelihood *tlk,
                                      const double *partials1, const double *matrices1,
                                      double *partials3) {
    double sum1;
    int v = 0;
    double *pPartials = partials3;

    for (int k = 0; k < tlk->pattern_count; k++) {
        int w = _cat_offset(tlk, k);

        sum1 = matrices1[w] * partials1[v];
        sum1 += matrices1[w + 1] * partials1[v + 1];
        sum1 += matrices1[w + 2] * partials1[v + 2];
        sum1 += matrices1[w + 3] * partials1[v + 3];

        *pPartials++ = sum1;

        sum1 = matrices1[w + 4] * partials1[v];
        sum1 += matrices1[w + 5] * partials1[v + 1];
        sum1 += matrices1[w + 6] * partials1[v + 2];
        sum1 += matrices1[w + 7] * partials1[v + 3];

        *pPartials++ = sum1;

        sum1 = matrices1[w + 8] * partials1[v];
        sum1 += matrices1[w + 9] * partials1[v + 1];
        sum1 += matrices1[w + 10] * partials1[v + 2];
        sum1 += matrices1[w + 11] * partials1[v + 3];

        *pPartials++ = sum1;

        sum1 = matrices1[w + 12] * partials1[v];
        sum1 += matrices1[w + 13] * partials1[v + 1];
        sum1 += matrices1[w + 14] * partials1[v + 2];
        sum1 += matrices1[w + 15] * partials1[v + 3];

        *pPartials++ = sum1;

        v += 4;
    }
}

void update_partials_4_cat(SingleTreeLikelihood *tlk, int partialsIndex,
                           int partialsIndex1, int matrixIndex1, int partialsIndex2,
                           int matrixIndex2) {
    // The negative index selects partials_undefined_t_and_undefined_4, a transposed
    // variant with no live caller in the lower likelihood. Refuse rather than fall
    // through to a kernel that would read the wrong matrix block.
    if (partialsIndex1 < 0) {
        fprintf(stderr,
                "update_partials_4_cat: the transposed path is not implemented for "
                "the CAT site model\n");
        exit(2);
    }

    if (partialsIndex2 < 0) {
        if (tlk->partials[0][partialsIndex1] != NULL) {
            _partials_undefined_4_cat(
                tlk,
                tlk->partials[tlk->current_partials_indexes[partialsIndex1]]
                             [partialsIndex1],
                tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]]
                             [matrixIndex1],
                tlk->partials[tlk->current_partials_indexes[partialsIndex]]
                             [partialsIndex]);
        } else {
            _partials_states_4_cat(
                tlk, tlk->mapping[partialsIndex1],
                tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]]
                             [matrixIndex1],
                tlk->partials[tlk->current_partials_indexes[partialsIndex]]
                             [partialsIndex]);
        }
    } else if (tlk->partials[0][partialsIndex1] != NULL) {
        if (tlk->partials[0][partialsIndex2] != NULL) {
            _partials_undefined_and_undefined_4_cat(
                tlk,
                tlk->partials[tlk->current_partials_indexes[partialsIndex1]]
                             [partialsIndex1],
                tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]]
                             [matrixIndex1],
                tlk->partials[tlk->current_partials_indexes[partialsIndex2]]
                             [partialsIndex2],
                tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]]
                             [matrixIndex2],
                tlk->partials[tlk->current_partials_indexes[partialsIndex]]
                             [partialsIndex]);
        } else {
            _partials_states_and_undefined_4_cat(
                tlk, tlk->mapping[partialsIndex2],
                tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]]
                             [matrixIndex2],
                tlk->partials[tlk->current_partials_indexes[partialsIndex1]]
                             [partialsIndex1],
                tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]]
                             [matrixIndex1],
                tlk->partials[tlk->current_partials_indexes[partialsIndex]]
                             [partialsIndex]);
        }
    } else {
        if (tlk->partials[0][partialsIndex2] != NULL) {
            _partials_states_and_undefined_4_cat(
                tlk, tlk->mapping[partialsIndex1],
                tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]]
                             [matrixIndex1],
                tlk->partials[tlk->current_partials_indexes[partialsIndex2]]
                             [partialsIndex2],
                tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]]
                             [matrixIndex2],
                tlk->partials[tlk->current_partials_indexes[partialsIndex]]
                             [partialsIndex]);
        } else {
            _partials_states_and_states_4_cat(
                tlk, tlk->mapping[partialsIndex1],
                tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]]
                             [matrixIndex1],
                tlk->mapping[partialsIndex2],
                tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]]
                             [matrixIndex2],
                tlk->partials[tlk->current_partials_indexes[partialsIndex]]
                             [partialsIndex]);
        }
    }

    if (tlk->scale) {
        SingleTreeLikelihood_scalePartials(tlk, partialsIndex, partialsIndex1,
                                           partialsIndex2);
    }
}

#pragma mark -
#pragma mark SSE

#ifdef SSE3_ENABLED

// The matrices of a node whose children are tip states are stored transposed (see the
// p_t_transpose calls in _calculate_partials), which is why the two kernels reading a
// state index below address them as matrices1[w + 4*state] instead of
// matrices1[w + state].

static void _partials_states_and_states_4_SSE_cat(const SingleTreeLikelihood *tlk,
                                                  int idx1, const double *matrices1,
                                                  int idx2, const double *matrices2,
                                                  double *partials) {
    const double *m1;
    const double *m2;
    double *pPartials = partials;
    __m128d m1v0, m1v2, m2v0, m2v2;

    for (int k = 0; k < tlk->pattern_count; k++) {
        int state1 = tlk->sp->patterns[idx1][k];
        int state2 = tlk->sp->patterns[idx2][k];
        int w = _cat_offset(tlk, k);

        if (state1 < 4 && state2 < 4) {
            m1 = &matrices1[w + 4 * state1];
            m2 = &matrices2[w + 4 * state2];

            m1v0 = _mm_load_pd(&m1[0]);
            m1v2 = _mm_load_pd(&m1[2]);

            m2v0 = _mm_load_pd(&m2[0]);
            m2v2 = _mm_load_pd(&m2[2]);

            _mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
            pPartials += 2;
            _mm_store_pd(pPartials, _mm_mul_pd(m1v2, m2v2));
            pPartials += 2;
        } else if (state1 < 4) {
            // child 1 has a gap or unknown state so treat it as unknown
            m1 = &matrices1[w + 4 * state1];

            *pPartials++ = *m1++;
            *pPartials++ = *m1++;
            *pPartials++ = *m1++;
            *pPartials++ = *m1;
        } else if (state2 < 4) {
            // child 2 has a gap or unknown state so treat it as unknown
            m2 = &matrices2[w + 4 * state2];

            *pPartials++ = *m2++;
            *pPartials++ = *m2++;
            *pPartials++ = *m2++;
            *pPartials++ = *m2;
        } else {
            // both children have a gap or unknown state so set partials to 1
            *pPartials++ = 1.0;
            *pPartials++ = 1.0;
            *pPartials++ = 1.0;
            *pPartials++ = 1.0;
        }
    }
}

static void _partials_states_4_SSE_cat(const SingleTreeLikelihood *tlk, int idx1,
                                       const double *matrices1, double *partials) {
    double *pPartials = partials;

    for (int k = 0; k < tlk->pattern_count; k++) {
        int state1 = tlk->sp->patterns[idx1][k];
        int w = _cat_offset(tlk, k);

        if (state1 < 4) {
            memcpy(pPartials, &matrices1[w + 4 * state1], 4 * sizeof(double));
            pPartials += 4;
        } else {
            // P.[1 1 1 1]^T = [1 1 1 1] when P is a probability matrix.
            // When derivatives are calculated P is NOT a probability matrix so we do
            // the calculation. P is transposed here, hence the stride of 4.
            *pPartials++ =
                matrices1[w] + matrices1[w + 4] + matrices1[w + 8] + matrices1[w + 12];
            *pPartials++ = matrices1[w + 1] + matrices1[w + 5] + matrices1[w + 9] +
                           matrices1[w + 13];
            *pPartials++ = matrices1[w + 2] + matrices1[w + 6] + matrices1[w + 10] +
                           matrices1[w + 14];
            *pPartials++ = matrices1[w + 3] + matrices1[w + 7] + matrices1[w + 11] +
                           matrices1[w + 15];
        }
    }
}

static void _partials_states_and_undefined_4_SSE_cat(const SingleTreeLikelihood *tlk,
                                                     int idx1, const double *matrices1,
                                                     const double *partials2,
                                                     const double *matrices2,
                                                     double *partials3) {
    int v = 0;
    __m128d p2v0, p2v2, m2v0, m2v2, *m1, temp;
    __m128d *pPartials = (__m128d *)partials3;

    for (int k = 0; k < tlk->pattern_count; k++) {
        int state1 = tlk->sp->patterns[idx1][k];
        int w = _cat_offset(tlk, k);

        p2v0 = _mm_load_pd(&partials2[v]);
        p2v2 = _mm_load_pd(&partials2[v + 2]);

        if (state1 < 4) {
            m1 = (__m128d *)&matrices1[w + state1 * 4];

            m2v0 = _mm_load_pd(&matrices2[w]);
            m2v2 = _mm_load_pd(&matrices2[w + 2]);

            m2v0 = _mm_mul_pd(m2v0, p2v0);
            m2v2 = _mm_mul_pd(m2v2, p2v2);
            temp = _mm_add_pd(m2v0, m2v2);

            m2v0 = _mm_load_pd(&matrices2[w + 4]);
            m2v2 = _mm_load_pd(&matrices2[w + 6]);

            m2v0 = _mm_mul_pd(m2v0, p2v0);
            m2v2 = _mm_mul_pd(m2v2, p2v2);

            m2v0 = _mm_add_pd(m2v0, m2v2);

            m2v2 = _mm_unpacklo_pd(temp, m2v0);
            temp = _mm_unpackhi_pd(temp, m2v0);
            *pPartials++ = _mm_mul_pd(*m1, _mm_add_pd(m2v2, temp));
            m1++;

            m2v0 = _mm_load_pd(&matrices2[w + 8]);
            m2v2 = _mm_load_pd(&matrices2[w + 10]);

            m2v0 = _mm_mul_pd(m2v0, p2v0);
            m2v2 = _mm_mul_pd(m2v2, p2v2);
            temp = _mm_add_pd(m2v0, m2v2);

            m2v0 = _mm_load_pd(&matrices2[w + 12]);
            m2v2 = _mm_load_pd(&matrices2[w + 14]);

            m2v0 = _mm_mul_pd(m2v0, p2v0);
            m2v2 = _mm_mul_pd(m2v2, p2v2);

            m2v0 = _mm_add_pd(m2v0, m2v2);

            m2v2 = _mm_unpacklo_pd(temp, m2v0);
            temp = _mm_unpackhi_pd(temp, m2v0);
            *pPartials++ = _mm_mul_pd(*m1, _mm_add_pd(m2v2, temp));
        } else {
            // Child 1 has a gap or unknown state so don't use it

            m2v0 = _mm_load_pd(&matrices2[w]);
            m2v2 = _mm_load_pd(&matrices2[w + 2]);

            m2v0 = _mm_mul_pd(m2v0, p2v0);
            m2v2 = _mm_mul_pd(m2v2, p2v2);
            temp = _mm_add_pd(m2v0, m2v2);

            m2v0 = _mm_load_pd(&matrices2[w + 4]);
            m2v2 = _mm_load_pd(&matrices2[w + 6]);

            m2v0 = _mm_mul_pd(m2v0, p2v0);
            m2v2 = _mm_mul_pd(m2v2, p2v2);

            m2v0 = _mm_add_pd(m2v0, m2v2);

            m2v2 = _mm_unpacklo_pd(temp, m2v0);
            temp = _mm_unpackhi_pd(temp, m2v0);
            *pPartials++ = _mm_add_pd(m2v2, temp);

            m2v0 = _mm_load_pd(&matrices2[w + 8]);
            m2v2 = _mm_load_pd(&matrices2[w + 10]);

            m2v0 = _mm_mul_pd(m2v0, p2v0);
            m2v2 = _mm_mul_pd(m2v2, p2v2);
            temp = _mm_add_pd(m2v0, m2v2);

            m2v0 = _mm_load_pd(&matrices2[w + 12]);
            m2v2 = _mm_load_pd(&matrices2[w + 14]);

            m2v0 = _mm_mul_pd(m2v0, p2v0);
            m2v2 = _mm_mul_pd(m2v2, p2v2);

            m2v0 = _mm_add_pd(m2v0, m2v2);

            m2v2 = _mm_unpacklo_pd(temp, m2v0);
            temp = _mm_unpackhi_pd(temp, m2v0);
            *pPartials++ = _mm_add_pd(m2v2, temp);
        }
        v += 4;
    }
}

static void _partials_undefined_and_undefined_4_SSE_cat(
    const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1,
    const double *partials2, const double *matrices2, double *partials3) {
    int v = 0;
    __m128d *pPartials = (__m128d *)partials3;
    __m128d m1v0, m1v2, m2v0, m2v2, p1v0, p1v2, p2v0, p2v2, temp;

    for (int k = 0; k < tlk->pattern_count; k++) {
        int w = _cat_offset(tlk, k);

        p1v0 = _mm_load_pd(&partials1[v]);
        p1v2 = _mm_load_pd(&partials1[v + 2]);

        p2v0 = _mm_load_pd(&partials2[v]);
        p2v2 = _mm_load_pd(&partials2[v + 2]);

        m1v0 = _mm_load_pd(&matrices1[w]);
        m1v2 = _mm_load_pd(&matrices1[w + 2]);

        m2v0 = _mm_load_pd(&matrices2[w]);
        m2v2 = _mm_load_pd(&matrices2[w + 2]);

        m1v0 = _mm_mul_pd(m1v0, p1v0);
        m1v2 = _mm_mul_pd(m1v2, p1v2);

        m1v0 = _mm_add_pd(m1v0, m1v2);

        m2v0 = _mm_mul_pd(m2v0, p2v0);
        m2v2 = _mm_mul_pd(m2v2, p2v2);

        m2v0 = _mm_add_pd(m2v0, m2v2);

        temp = _mm_hadd_pd(m1v0, m2v0);

        m1v0 = _mm_load_pd(&matrices1[w + 4]);
        m1v2 = _mm_load_pd(&matrices1[w + 6]);

        m2v0 = _mm_load_pd(&matrices2[w + 4]);
        m2v2 = _mm_load_pd(&matrices2[w + 6]);

        m1v0 = _mm_mul_pd(m1v0, p1v0);
        m1v2 = _mm_mul_pd(m1v2, p1v2);

        m1v0 = _mm_add_pd(m1v0, m1v2);

        m2v0 = _mm_mul_pd(m2v0, p2v0);
        m2v2 = _mm_mul_pd(m2v2, p2v2);

        m2v0 = _mm_add_pd(m2v0, m2v2);

        m1v0 = _mm_hadd_pd(m1v0, m2v0);

        m2v0 = _mm_unpacklo_pd(temp, m1v0);
        temp = _mm_unpackhi_pd(temp, m1v0);
        *pPartials++ = _mm_mul_pd(m2v0, temp);

        m1v0 = _mm_load_pd(&matrices1[w + 8]);
        m1v2 = _mm_load_pd(&matrices1[w + 10]);

        m2v0 = _mm_load_pd(&matrices2[w + 8]);
        m2v2 = _mm_load_pd(&matrices2[w + 10]);

        m1v0 = _mm_mul_pd(m1v0, p1v0);
        m1v2 = _mm_mul_pd(m1v2, p1v2);

        m1v0 = _mm_add_pd(m1v0, m1v2);

        m2v0 = _mm_mul_pd(m2v0, p2v0);
        m2v2 = _mm_mul_pd(m2v2, p2v2);

        m2v0 = _mm_add_pd(m2v0, m2v2);

        temp = _mm_hadd_pd(m1v0, m2v0);

        m1v0 = _mm_load_pd(&matrices1[w + 12]);
        m1v2 = _mm_load_pd(&matrices1[w + 14]);

        m2v0 = _mm_load_pd(&matrices2[w + 12]);
        m2v2 = _mm_load_pd(&matrices2[w + 14]);

        m1v0 = _mm_mul_pd(m1v0, p1v0);
        m1v2 = _mm_mul_pd(m1v2, p1v2);

        m1v0 = _mm_add_pd(m1v0, m1v2);

        m2v0 = _mm_mul_pd(m2v0, p2v0);
        m2v2 = _mm_mul_pd(m2v2, p2v2);

        m2v0 = _mm_add_pd(m2v0, m2v2);

        m1v0 = _mm_hadd_pd(m1v0, m2v0);

        m2v0 = _mm_unpacklo_pd(temp, m1v0);
        temp = _mm_unpackhi_pd(temp, m1v0);
        *pPartials++ = _mm_mul_pd(m2v0, temp);

        v += 4;
    }
}

static void _partials_undefined_4_SSE_cat(const SingleTreeLikelihood *tlk,
                                          const double *partials1,
                                          const double *matrices1, double *partials3) {
    int v = 0;
    __m128d *pPartials = (__m128d *)partials3;
    __m128d m1v0, m1v2, p1v0, p1v2, temp;

    for (int k = 0; k < tlk->pattern_count; k++) {
        int w = _cat_offset(tlk, k);

        p1v0 = _mm_load_pd(&partials1[v]);
        p1v2 = _mm_load_pd(&partials1[v + 2]);

        m1v0 = _mm_load_pd(&matrices1[w]);
        m1v2 = _mm_load_pd(&matrices1[w + 2]);

        m1v0 = _mm_mul_pd(m1v0, p1v0);
        m1v2 = _mm_mul_pd(m1v2, p1v2);
        temp = _mm_add_pd(m1v0, m1v2);

        m1v0 = _mm_load_pd(&matrices1[w + 4]);
        m1v2 = _mm_load_pd(&matrices1[w + 6]);

        m1v0 = _mm_mul_pd(m1v0, p1v0);
        m1v2 = _mm_mul_pd(m1v2, p1v2);
        m1v0 = _mm_add_pd(m1v0, m1v2);

        m1v2 = _mm_unpacklo_pd(temp, m1v0);
        temp = _mm_unpackhi_pd(temp, m1v0);
        *pPartials++ = _mm_add_pd(m1v2, temp);

        m1v0 = _mm_load_pd(&matrices1[w + 8]);
        m1v2 = _mm_load_pd(&matrices1[w + 10]);

        m1v0 = _mm_mul_pd(m1v0, p1v0);
        m1v2 = _mm_mul_pd(m1v2, p1v2);
        temp = _mm_add_pd(m1v0, m1v2);

        m1v0 = _mm_load_pd(&matrices1[w + 12]);
        m1v2 = _mm_load_pd(&matrices1[w + 14]);
        m1v0 = _mm_mul_pd(m1v0, p1v0);
        m1v2 = _mm_mul_pd(m1v2, p1v2);

        m1v0 = _mm_add_pd(m1v0, m1v2);

        m1v2 = _mm_unpacklo_pd(temp, m1v0);
        temp = _mm_unpackhi_pd(temp, m1v0);
        *pPartials++ = _mm_add_pd(m1v2, temp);

        v += 4;
    }
}

void update_partials_4_SSE_cat(SingleTreeLikelihood *tlk, int partialsIndex,
                               int partialsIndex1, int matrixIndex1, int partialsIndex2,
                               int matrixIndex2) {
    if (partialsIndex2 < 0) {
        if (tlk->partials[0][partialsIndex1] != NULL) {
            _partials_undefined_4_SSE_cat(
                tlk,
                tlk->partials[tlk->current_partials_indexes[partialsIndex1]]
                             [partialsIndex1],
                tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]]
                             [matrixIndex1],
                tlk->partials[tlk->current_partials_indexes[partialsIndex]]
                             [partialsIndex]);
        } else {
            _partials_states_4_SSE_cat(
                tlk, tlk->mapping[partialsIndex1],
                tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]]
                             [matrixIndex1],
                tlk->partials[tlk->current_partials_indexes[partialsIndex]]
                             [partialsIndex]);
        }
    } else if (tlk->partials[0][partialsIndex1] != NULL) {
        if (tlk->partials[0][partialsIndex2] != NULL) {
            _partials_undefined_and_undefined_4_SSE_cat(
                tlk,
                tlk->partials[tlk->current_partials_indexes[partialsIndex1]]
                             [partialsIndex1],
                tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]]
                             [matrixIndex1],
                tlk->partials[tlk->current_partials_indexes[partialsIndex2]]
                             [partialsIndex2],
                tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]]
                             [matrixIndex2],
                tlk->partials[tlk->current_partials_indexes[partialsIndex]]
                             [partialsIndex]);
        } else {
            _partials_states_and_undefined_4_SSE_cat(
                tlk, tlk->mapping[partialsIndex2],
                tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]]
                             [matrixIndex2],
                tlk->partials[tlk->current_partials_indexes[partialsIndex1]]
                             [partialsIndex1],
                tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]]
                             [matrixIndex1],
                tlk->partials[tlk->current_partials_indexes[partialsIndex]]
                             [partialsIndex]);
        }
    } else {
        if (tlk->partials[0][partialsIndex2] != NULL) {
            _partials_states_and_undefined_4_SSE_cat(
                tlk, tlk->mapping[partialsIndex1],
                tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]]
                             [matrixIndex1],
                tlk->partials[tlk->current_partials_indexes[partialsIndex2]]
                             [partialsIndex2],
                tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]]
                             [matrixIndex2],
                tlk->partials[tlk->current_partials_indexes[partialsIndex]]
                             [partialsIndex]);
        } else {
            _partials_states_and_states_4_SSE_cat(
                tlk, tlk->mapping[partialsIndex1],
                tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]]
                             [matrixIndex1],
                tlk->mapping[partialsIndex2],
                tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]]
                             [matrixIndex2],
                tlk->partials[tlk->current_partials_indexes[partialsIndex]]
                             [partialsIndex]);
        }
    }

    if (tlk->scale) {
        SingleTreeLikelihood_scalePartials(tlk, partialsIndex, partialsIndex1,
                                           partialsIndex2);
    }
}

#endif
