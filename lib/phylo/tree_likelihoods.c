/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: tree_likelihoods.c,v 1.14 2008-11-12 02:07:59 acs Exp $ */

#include <tree_likelihoods.h>
#include <subst_mods.h>
#include <markov_matrix.h>
#include <subst_mods.h>
#include <dgamma.h>
#include <sufficient_stats.h>

/* Computation of likelihoods for columns of a given multiple
   alignment, according to a given tree model.  */

/* NOTE: routines currently assume that the names at the leaves of the
   tree (as defined in the TreeModel object) are numbers from 1 to
   nseqs.  */

/* FIXME: inside and outside computations should be in log space --
   probably only an issue when the number of leaves is large */
void printTree(TreeNode* node);
int tuple_index_missing_data(char *tuple, int *inv_alph, int *is_missing,
                             int alph_size);



/* Compute the likelihood of a tree model with respect to an
   alignment.  Optionally retain column-by-column likelihoods,
   optionally compute posterior probabilities.  If 'post' is NULL, no
   posterior probabilities (or related quantities) will be computed.
   If 'post' is non-NULL each of its attributes must either be NULL or
   previously allocated to the required size. */
double tl_compute_log_likelihood(TreeModel *mod, MSA *msa, 
                                 double *col_scores, double *tuple_scores,
				 int cat, TreePosteriors *post) {

  int i, j;
  double retval = 0;
  int nstates = mod->rate_matrix->size;
  int alph_size = strlen(mod->rate_matrix->states); 
  int npasses = (mod->order > 0 && mod->use_conditionals == 1 ? 2 : 1); 
  int pass, col_offset, k, nodeidx, rcat, /* colidx, */ tupleidx, defined;
  TreeNode *n;
  double total_prob, marg_tot;
  List *traversal;
  double **inside_joint = NULL, **inside_marginal = NULL, 
    **outside_joint = NULL, **outside_marginal = NULL, 
    ****subst_probs = NULL;
  double *curr_tuple_scores=NULL;
  double rcat_prob[mod->nratecats];
  double tmp[nstates];

  checkInterrupt();

  /* allocate memory */
  inside_joint = (double**)smalloc(nstates * sizeof(double*));
  for (j = 0; j < nstates; j++) 
    inside_joint[j] = (double*)smalloc((mod->tree->nnodes+1) * 
                                       sizeof(double)); 
  outside_joint = (double**)smalloc(nstates * sizeof(double*));
  for (j = 0; j < nstates; j++) 
    outside_joint[j] = (double*)smalloc((mod->tree->nnodes+1) * 
                                        sizeof(double)); 
  /* only needed if post != NULL? */
  if (mod->order > 0) {
    inside_marginal = (double**)smalloc(nstates * sizeof(double*));
    for (j = 0; j < nstates; j++) 
      inside_marginal[j] = (double*)smalloc((mod->tree->nnodes+1) * 
                                            sizeof(double));
  }
  if (mod->order > 0 && post != NULL) {
    outside_marginal = (double**)smalloc(nstates * sizeof(double*));
    for (j = 0; j < nstates; j++) 
      outside_marginal[j] = (double*)smalloc((mod->tree->nnodes+1) * 
                                             sizeof(double));
  }
  if (post != NULL) {
    subst_probs = (double****)smalloc(mod->nratecats * sizeof(double***));
    for (rcat = 0; rcat < mod->nratecats; rcat++) {
      subst_probs[rcat] = (double***)smalloc(nstates * sizeof(double**));
      for (j = 0; j < nstates; j++) {
        subst_probs[rcat][j] = (double**)smalloc(nstates * sizeof(double*));
        for (k = 0; k < nstates; k++)
          subst_probs[rcat][j][k] = (double*)smalloc(mod->tree->nnodes * sizeof(double));          
      }
    }
  }

  /* create IUPAC mapping if needed */
  if (mod->iupac_inv_map == NULL)
    mod->iupac_inv_map = build_iupac_inv_map(mod->rate_matrix->inv_states, 
                                             alph_size);

  
  if (cat > msa->ncats)
    die("ERROR tl_compute_log_likelihood: cat (%i) > msa->ncats (%i)\n", cat, msa->ncats);

  if (!(cat < 0 || col_scores == NULL || msa->categories != NULL))
    die("ERROR tl_compute_log_likelihood: cat=%i, col_scores==NULL=%i, msa->categories==NULL=%i\n", cat, col_scores==NULL, msa->categories==NULL);
  /* if using categories and col-by-col
     scoring, then must have col-by-col
     categories */

  /* obtain sufficient statistics, if necessary */
  if (msa->ss != NULL){ 
    if (msa->ss->tuple_size <= mod->order)
      die("ERROR tl_compute_log_likelihood: tuple_size (%i) must be greater than mod->order (%i)\n",
	  msa->ss->tuple_size, mod->order);
  }
  else 
    ss_from_msas(msa, mod->order+1, col_scores == NULL ? 0 : 1, 
                 NULL, NULL, NULL, -1, subst_mod_is_codon_model(mod->subst_mod));

  /* set up leaf to sequence mapping, if necessary */
  if (mod->msa_seq_idx == NULL)
    tm_build_seq_idx(mod, msa);

  /* set up prob matrices, if any are undefined */
  for (i = 0, defined = TRUE; defined && i < mod->tree->nnodes; i++) {
    if (((TreeNode*)lst_get_ptr(mod->tree->nodes, i))->parent == NULL) 
      continue;  		/* skip root */
    for (j = 0; j < mod->nratecats; j++)
      if (mod->P[i][j] == NULL) defined = FALSE;
  }
/*   printf("mod %d\n", mod); */
  if (!defined) {
    tm_set_subst_matrices(mod);
  }
  if (col_scores != NULL && tuple_scores == NULL)
    curr_tuple_scores = (double*)smalloc(msa->ss->ntuples * sizeof(double));
  else if (tuple_scores != NULL) 
    curr_tuple_scores = tuple_scores;
  if (curr_tuple_scores != NULL) 
    for (tupleidx = 0; tupleidx < msa->ss->ntuples; tupleidx++)
      curr_tuple_scores[tupleidx] = 0;

  if (post != NULL && post->expected_nsubst_tot != NULL) {
    for (rcat = 0; rcat < mod->nratecats; rcat++)
      for (i = 0; i < nstates; i++) 
        for (j = 0; j < nstates; j++) 
          for (k = 0; k < mod->tree->nnodes; k++)
            post->expected_nsubst_tot[rcat][i][j][k] = 0;
  }
  if (post != NULL && post->rcat_expected_nsites != NULL)
    for (rcat = 0; rcat < mod->nratecats; rcat++)
      post->rcat_expected_nsites[rcat] = 0;

  
  /*Iterate over every column in MSA, computer likelihood L(i) for ith column. Add all
   likelihoods for overall tree.*/
  for (tupleidx = 0; tupleidx < msa->ss->ntuples; tupleidx++) {
    int skip_fels = FALSE;

    if ((cat >= 0 && msa->ss->cat_counts[cat][tupleidx] == 0) || 
        (cat < 0 && msa->ss->counts[tupleidx] == 0))
      continue;
    checkInterruptN(tupleidx, 1000);

    total_prob = 0;
    marg_tot = NULL_LOG_LIKELIHOOD;

    /* check for gaps and whether column is informative, if necessary */
    if (!mod->allow_gaps)
      for (j = 0; !skip_fels && j < msa->nseqs; j++) 
        if (ss_get_char_tuple(msa, tupleidx, j, 0) == GAP_CHAR) 
          skip_fels = TRUE;
    if (!skip_fels && mod->inform_reqd) {
      int ninform = 0;
      for (j = 0; j < msa->nseqs; j++) {
        if (msa->is_informative != NULL && !msa->is_informative[j])
          continue;
        else if (!msa->is_missing[(int)ss_get_char_tuple(msa, tupleidx, j, 0)])
          ninform++;
      }
      if (ninform < 2) skip_fels = TRUE;
    }
          
    if (!skip_fels) {
      for (pass = 0; pass < npasses; pass++) {
        double **pL = (pass == 0 ? inside_joint : inside_marginal);
        double **pLbar = (pass == 0 ? outside_joint : outside_marginal);
        /*         TreePosteriors *postpass = (pass == 0 ? post : postmarg); */

        if (pass > 0)
          marg_tot = 0;         /* will need to compute */

        for (rcat = 0; rcat < mod->nratecats; rcat++) {
          traversal = tr_postorder(mod->tree);      
          for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
            int partial_match[mod->order+1][alph_size];
            n = lst_get_ptr(traversal, nodeidx);      
            if (n->lchild == NULL) { 
              /* leaf: base case of recursion */
              int thisseq;

	      if (n->name == NULL)
		die("ERROR tl_compute_log_likelihood: n->name is NULL\n");
              thisseq = mod->msa_seq_idx[n->id];
	      if (thisseq < 0) 
		die("ERROR tl_compute_log_likelihood: expected a leaf node\n");

              /* first figure out whether there is a match for each
                 character in each position; we'll call this the record of
                 "partial_matches". */
              for (col_offset = -1*mod->order; col_offset <= 0; col_offset++) {
                int observed_state = -1;
                int *iupac_prob = NULL;

                if (pass == 0 || col_offset < 0) {
                  char thischar = (int)ss_get_char_tuple(msa, tupleidx, 
                                                         thisseq, col_offset);
                  observed_state = mod->rate_matrix->inv_states[(int)thischar];
                  if (observed_state < 0)
                    iupac_prob = mod->iupac_inv_map[(int)thischar];
                }
 
                /* otherwise, we're on a second pass and looking the
                   current base, so we want to use the "missing
                   information" principle */

                if (iupac_prob != NULL) {
                  for (i = 0; i < alph_size; i++) 
                    partial_match[mod->order+col_offset][i] = iupac_prob[i];
                }
                else {
                  for (i = 0; i < alph_size; i++) {
                    if (observed_state < 0 || i == observed_state) 
                      partial_match[mod->order+col_offset][i] = 1;
                    else
                      partial_match[mod->order+col_offset][i] = 0; 
                  }
                }
              }

              /* now find the intersection of the partial matches */
              for (i = 0; i < nstates; i++) {
                if (mod->order == 0)  /* handle 0th order model as special
                                         case, for efficiency.  In this case
                                         the partial match *is* the total
                                         match */
                  pL[i][n->id] = partial_match[0][i];
                else {
                  int total_match = 1;
                  /* figure out the "projection" of state i in the dimension
                     of each position, and see whether there is a
                     corresponding partial match. */
                  /* NOTE: mod->order is approx equal to log nstates
                     (prob no more than 2) */
                  for (col_offset = -1*mod->order; col_offset <= 0 && total_match; 
                       col_offset++) {
                    int projection = (i / int_pow(alph_size, -1 * col_offset)) % 
                      alph_size;

                    if (!partial_match[mod->order+col_offset][projection]) 
                      total_match = 0; /* must have partial matches in all
                                          dimensions for a total match */
                  }
                  pL[i][n->id] = total_match;
                }
              }
            }
            else {                    
              /* general recursive case */
              MarkovMatrix *lsubst_mat = mod->P[n->lchild->id][rcat];
              MarkovMatrix *rsubst_mat = mod->P[n->rchild->id][rcat];
              for (i = 0; i < nstates; i++) {
                double totl = 0, totr = 0;
                for (j = 0; j < nstates; j++) 
                  totl += pL[j][n->lchild->id] *
                    mm_get(lsubst_mat, i, j);
        
                for (k = 0; k < nstates; k++) 
                  totr += pL[k][n->rchild->id] *
                    mm_get(rsubst_mat, i, k);

                pL[i][n->id] = totl * totr;
              }
            }
          }

          if (post != NULL && pass == 0) {
            MarkovMatrix *subst_mat;
            double this_total, denom;

            /* do outside calculation */
            traversal = tr_preorder(mod->tree);
            for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
              n = lst_get_ptr(traversal, nodeidx);
              if (n->parent == NULL) { /* base case */ 
                for (i = 0; i < nstates; i++) 
                  pLbar[i][n->id] = vec_get(mod->backgd_freqs, i);
              }
              else {            /* recursive case */
                TreeNode *sibling = (n == n->parent->lchild ? 
                                     n->parent->rchild : n->parent->lchild);
                MarkovMatrix *par_subst_mat = mod->P[n->id][rcat];
                MarkovMatrix *sib_subst_mat = mod->P[sibling->id][rcat];

                /* breaking this computation into two parts as follows
                   reduces its complexity by a factor of nstates */

                for (j = 0; j < nstates; j++) { /* parent state */
                  tmp[j] = 0;
                  for (k = 0; k < nstates; k++) { /* sibling state */
                    tmp[j] += pLbar[j][n->parent->id] *
                      pL[k][sibling->id] * mm_get(sib_subst_mat, j, k);
                  }
                }

                for (i = 0; i < nstates; i++) { /* child state */
                  pLbar[i][n->id] = 0;
                  for (j = 0; j < nstates; j++) { /* parent state */
                    pLbar[i][n->id] += 
                      tmp[j] * mm_get(par_subst_mat, j, i);
                  }
                }
              }
              

              /* compute total probability based on current node, to
                 avoid numerical errors */
              this_total = 0;
              for (i = 0; i < nstates; i++)
                this_total += pL[i][n->id] * pLbar[i][n->id];

              if (post->expected_nsubst != NULL && n->parent != NULL)
                post->expected_nsubst[rcat][n->id][tupleidx] = 1;

              subst_mat = mod->P[n->id][rcat];
              for (i = 0; i < nstates; i++) {
                /* compute posterior prob of base (tuple) i at node n */
                if (post->base_probs != NULL) {
                  post->base_probs[rcat][i][n->id][tupleidx] = 
                    safediv(pL[i][n->id] * pLbar[i][n->id], this_total);
                }

                if (n->parent == NULL) continue;

                /* (intermediate computation used for subst probs) */
                denom = 0;
                for (k = 0; k < nstates; k++)
                  denom += pL[k][n->id] * mm_get(subst_mat, i, k);

                for (j = 0; j < nstates; j++) { 
                  /* compute posterior prob of a subst of base j at
                     node n for base i at node n->parent */
                  subst_probs[rcat][i][j][n->id] = 
                    safediv(pL[i][n->parent->id] * pLbar[i][n->parent->id], 
                            this_total) *
                    pL[j][n->id] * mm_get(subst_mat, i, j);
                  subst_probs[rcat][i][j][n->id] = 
                    safediv(subst_probs[rcat][i][j][n->id], denom);
                    
                  if (post->subst_probs != NULL)
                    post->subst_probs[rcat][i][j][n->id][tupleidx] = 
                      subst_probs[rcat][i][j][n->id];
                  
                  if (post->expected_nsubst != NULL && j == i)
                    post->expected_nsubst[rcat][n->id][tupleidx] -= 
                      subst_probs[rcat][i][j][n->id];

                }
              }
            }
          }

          if (pass == 0) {
            rcat_prob[rcat] = 0;
            for (i = 0; i < nstates; i++) {
              rcat_prob[rcat] += vec_get(mod->backgd_freqs, i) * 
                inside_joint[i][mod->tree->id] * mod->freqK[rcat];
            }
            total_prob += rcat_prob[rcat];
          }
          else { 
            for (i = 0; i < nstates; i++) 
              marg_tot += vec_get(mod->backgd_freqs, i) * 
                inside_marginal[i][mod->tree->id] * mod->freqK[rcat]; 
          }
        } /* for rcat */
      } /* for pass */
    } /* if skip_fels */

      /* compute posterior prob of each rate cat and related quantities */
    if (post != NULL) {
      if (skip_fels) die("ERROR: tl_compute_log_likelihood: skip_fels should be 0 but is %i\n", skip_fels);
      for (rcat = 0; rcat < mod->nratecats; rcat++) {
        double rcat_post_prob = safediv(rcat_prob[rcat], total_prob);
        if (post->rcat_probs != NULL) 
          post->rcat_probs[rcat][tupleidx] = rcat_post_prob;
        if (post->rcat_expected_nsites != NULL)
          post->rcat_expected_nsites[rcat] += rcat_post_prob * 
            (cat >= 0 ? msa->ss->cat_counts[cat][tupleidx] : 
             msa->ss->counts[tupleidx]);
        if (post->expected_nsubst_tot != NULL) {
          for (nodeidx = 0; nodeidx < mod->tree->nnodes; nodeidx++) {
            n = lst_get_ptr(mod->tree->nodes, nodeidx);
            if (n->parent == NULL) continue;
            for (i = 0; i < nstates; i++) 
              for (j = 0; j < nstates; j++) 
                post->expected_nsubst_tot[rcat][i][j][n->id] += 
                  subst_probs[rcat][i][j][n->id] * 
                  (cat >= 0 ? msa->ss->cat_counts[cat][tupleidx] : 
                   msa->ss->counts[tupleidx]) *
                  rcat_post_prob;
          }
        }
	if (post->expected_nsubst_col != NULL) {
	  for (nodeidx = 0; nodeidx < mod->tree->nnodes; nodeidx++) {
            n = lst_get_ptr(mod->tree->nodes, nodeidx);
            if (n->parent == NULL) continue;
            for (i = 0; i < nstates; i++) 
              for (j = 0; j < nstates; j++) 
                post->expected_nsubst_col[rcat][n->id][tupleidx][i][j] = 
                  subst_probs[rcat][i][j][n->id] * rcat_post_prob;
          }
        }
      }
    }

    if (mod->order > 0 && mod->use_conditionals == 1 && !skip_fels) 
      total_prob /= marg_tot; 

    /*    if (total_prob > 1.0) {
      if (total_prob - 1.0 < 1.0e-6) total_prob = 1.0;
      else die("got total_prob=%.10g\n", total_prob);
      }*/
    total_prob = log2(total_prob);

    if (curr_tuple_scores != NULL && 
        (cat < 0 || msa->ss->cat_counts[cat][tupleidx] > 0))
      curr_tuple_scores[tupleidx] = total_prob;
    /* NOTE: curr_tuple_scores contains the
       (log) probabilities *unweighted* by tuple counts */

    total_prob *= (cat >= 0 ? msa->ss->cat_counts[cat][tupleidx] : 
                   msa->ss->counts[tupleidx]); /* log space */

    retval += total_prob;     /* log space */
        
  } /* for tupleidx */

  for (j = 0; j < nstates; j++) {
    sfree(inside_joint[j]);
    sfree(outside_joint[j]);
    if (mod->order > 0) sfree(inside_marginal[j]);
    if (mod->order > 0 && post != NULL) sfree(outside_marginal[j]);
  }
  sfree(inside_joint);
  sfree(outside_joint);
  if (mod->order > 0) sfree(inside_marginal);
  if (mod->order > 0 && post != NULL) sfree(outside_marginal);
  if (col_scores != NULL) {
    if (cat >= 0) 
      for (i = 0; i < msa->length; i++)
        col_scores[i] = msa->categories[i] == cat ?
          curr_tuple_scores[msa->ss->tuple_idx[i]] :
          NEGINFTY;
    else
      for (i = 0; i < msa->length; i++)
        col_scores[i] = curr_tuple_scores[msa->ss->tuple_idx[i]];
    if (tuple_scores == NULL) sfree(curr_tuple_scores);
  }
  if (post != NULL) {
    for (rcat = 0; rcat < mod->nratecats; rcat++) {
      for (j = 0; j < nstates; j++) {
        for (k = 0; k < nstates; k++)
          sfree(subst_probs[rcat][j][k]);
        sfree(subst_probs[rcat][j]);
      }
      sfree (subst_probs[rcat]);
    }
    sfree(subst_probs);
  }
  return(retval);
}
/* =====================================================================================*/
/* Compute the likelihood of a tree model with respect to an
   alignment. Similar to above except uses an extended Felsestein
 * Tree Prunning algorithm as detailed by: Probabilistic Phylogenetic Inference 
 * with Insertions and Deletions (Rivas E, Eddy SR (2008)).
 * Optionally retain column-by-column likelihoods and/or posterior probabilities.  
   @param[in] mod Tree Model to compute likelihood for
   @param[in] msa Multiple Alignment containing data related to tree model
   @param[out] col_scores (Optional) Log likelihood score per column
   @param[out] tuple_scores (Optional) Log likelihood score per tuple
   @param[in] cat Whether to use categories
   @param[out] post (Optional) Computed posterior probabilities; If NULL, no
   posterior probabilities (or related quantities) will be computed.
   If non-NULL each of its attributes must either be NULL or
   previously allocated to the required size. 
   @result Log likelihood of entire tree model.
 *  */
double gapAwareLikelihood(TreeModel *mod, MSA *msa,double *col_scores, double *tuple_scores,
				 int cat, TreePosteriors *post) {
  int i, j,k;
  double retval = 0;
  int nstates = mod->rate_matrix->size;
  int alph_size = strlen(mod->rate_matrix->states); 
  int tupleidx;
  int rcat;

  double **inside_joint = NULL, **inside_marginal = NULL, 
    **outside_joint = NULL, **outside_marginal = NULL, 
    ****subst_probs = NULL;
  double *curr_tuple_scores=NULL;

  checkInterrupt();

  /* allocate memory */
  inside_joint = (double**)smalloc(nstates * sizeof(double*));
  for (j = 0; j < nstates; j++) 
    inside_joint[j] = (double*)smalloc((mod->tree->nnodes+1) * 
                                       sizeof(double)); 
  outside_joint = (double**)smalloc(nstates * sizeof(double*));
  for (j = 0; j < nstates; j++) 
    outside_joint[j] = (double*)smalloc((mod->tree->nnodes+1) * 
                                        sizeof(double)); 
  /* only needed if post != NULL? */
  if (mod->order > 0) {
    inside_marginal = (double**)smalloc(nstates * sizeof(double*));
    for (j = 0; j < nstates; j++) 
      inside_marginal[j] = (double*)smalloc((mod->tree->nnodes+1) * 
                                            sizeof(double));
  }
  if (mod->order > 0 && post != NULL) {
    outside_marginal = (double**)smalloc(nstates * sizeof(double*));
    for (j = 0; j < nstates; j++) 
      outside_marginal[j] = (double*)smalloc((mod->tree->nnodes+1) * 
                                             sizeof(double));
  }
  if (post != NULL) {
    subst_probs = (double****)smalloc(mod->nratecats * sizeof(double***));
    for (rcat = 0; rcat < mod->nratecats; rcat++) {
      subst_probs[rcat] = (double***)smalloc(nstates * sizeof(double**));
      for (j = 0; j < nstates; j++) {
        subst_probs[rcat][j] = (double**)smalloc(nstates * sizeof(double*));
        for (k = 0; k < nstates; k++)
          subst_probs[rcat][j][k] = (double*)smalloc(mod->tree->nnodes * sizeof(double));          
      }
    }
  }

  /* create IUPAC mapping if needed */
  if (mod->iupac_inv_map == NULL)
    mod->iupac_inv_map = build_iupac_inv_map(mod->rate_matrix->inv_states, alph_size);
  
  if (cat > msa->ncats)
    die("ERROR tl_compute_log_likelihood: cat (%i) > msa->ncats (%i)\n", cat, msa->ncats);

  if (!(cat < 0 || col_scores == NULL || msa->categories != NULL))
    die("ERROR tl_compute_log_likelihood: cat=%i, col_scores==NULL=%i, msa->categories==NULL=%i\n", cat, col_scores==NULL, msa->categories==NULL);
  /* if using categories and col-by-col scoring, then must have col-by-col categories */

  /* obtain sufficient statistics, if necessary */
  if (msa->ss != NULL){ 
    if (msa->ss->tuple_size <= mod->order)
      die("ERROR tl_compute_log_likelihood: tuple_size (%i) must be greater than mod->order (%i)\n",
	  msa->ss->tuple_size, mod->order);
  }
  else 
    ss_from_msas(msa, mod->order+1, col_scores == NULL ? 0 : 1, 
                 NULL, NULL, NULL, -1, subst_mod_is_codon_model(mod->subst_mod));

  /* set up leaf to sequence mapping, if necessary */
  if (mod->msa_seq_idx == NULL)
    tm_build_seq_idx(mod, msa);

  if (col_scores != NULL && tuple_scores == NULL)
    curr_tuple_scores = (double*)smalloc(msa->ss->ntuples * sizeof(double));
  else if (tuple_scores != NULL) 
    curr_tuple_scores = tuple_scores;
  if (curr_tuple_scores != NULL) 
    for (tupleidx = 0; tupleidx < msa->ss->ntuples; tupleidx++)
      curr_tuple_scores[tupleidx] = 0;

  if (post != NULL && post->expected_nsubst_tot != NULL) {
    for (rcat = 0; rcat < mod->nratecats; rcat++)
      for (i = 0; i < nstates; i++) 
        for (j = 0; j < nstates; j++) 
          for (k = 0; k < mod->tree->nnodes; k++)
            post->expected_nsubst_tot[rcat][i][j][k] = 0;
  }
  if (post != NULL && post->rcat_expected_nsites != NULL)
    for (rcat = 0; rcat < mod->nratecats; rcat++)
      post->rcat_expected_nsites[rcat] = 0;

  /**Do work here!*/
  retval = computeTotalTreeLikelihood(mod,msa,cat,post,inside_joint);

  for (j = 0; j < nstates; j++) {
    sfree(inside_joint[j]);
    sfree(outside_joint[j]);
    if (mod->order > 0) sfree(inside_marginal[j]);
    if (mod->order > 0 && post != NULL) sfree(outside_marginal[j]);
  }
  sfree(inside_joint);
  sfree(outside_joint);
  if (mod->order > 0) sfree(inside_marginal);
  if (mod->order > 0 && post != NULL) sfree(outside_marginal);
  if (col_scores != NULL) {
    if (cat >= 0) 
      for (i = 0; i < msa->length; i++)
        col_scores[i] = msa->categories[i] == cat ?
          curr_tuple_scores[msa->ss->tuple_idx[i]] :
          NEGINFTY;
    else
      for (i = 0; i < msa->length; i++)
        col_scores[i] = curr_tuple_scores[msa->ss->tuple_idx[i]];
    if (tuple_scores == NULL) sfree(curr_tuple_scores);
  }
  if (post != NULL) {
    for (rcat = 0; rcat < mod->nratecats; rcat++) {
      for (j = 0; j < nstates; j++) {
        for (k = 0; k < nstates; k++)
          sfree(subst_probs[rcat][j][k]);
        sfree(subst_probs[rcat][j]);
      }
      sfree (subst_probs[rcat]);
    }
    sfree(subst_probs);
  }
  return(retval);
}
/* =====================================================================================*/
/*Given a TreeModel, an MSA, the cat number, Tree Posteriors. Will return the total
 * likelihood for a tree. The rest of parameters are allocated memory. See gapAwareLikelihood
 * for allocations. This is function does the bulk of the work.
 * Note this is only guaranteed to work and tested on models of order zero with no column
 * offset and no rate categories set. */
  
double computeTotalTreeLikelihood(TreeModel* mod,MSA* msa,int cat,TreePosteriors* post,
        double **inside_joint){
  int tupleidx, i;
  int nstates = mod->rate_matrix->size;
  double total_prob;
  List* traversal;
  int nodeidx;
  double allGapProb;
  TreeNode *n;
  int alph_size = strlen(mod->rate_matrix->states);
  double retval = 0;
  double p = mod->geometricParameter;
      int rootNodeId = mod->tree->id;
  
  /*Diagonal Matrix will be the same for all sites as it doesn't change until we go
   *    back to optimization script, create and use for all sites.*/
  DiagonalMatrix* dMatrix = getDiagonalMatrix(mod->rate_matrix);
  
  /*Iterate over every column in MSA, compute likelihood L(i) for ith column. Add all
   *likelihoods for overall tree. Plus one since we need to take into account an "all 
   * gap" column probability. */
  for (tupleidx = 0; tupleidx < msa->ss->ntuples + 1; tupleidx ++) {    
    checkInterruptN(tupleidx, 1000);
    
    total_prob = 0;
    
    double **pL = inside_joint;
    
    /*Get traversal order so we iterate over nodes instead of recursing.*/
    traversal = tr_postorder(mod->tree);
    /* Iterate over traversal hitting all nodes in a post order matter.*/
    for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
      n = lst_get_ptr(traversal, nodeidx);

      /* Leaf: base case of recursion */
      if (n->lchild == NULL) {

        int sequenceNumber = mod->msa_seq_idx[n->id];
        /*Integer version of character in our alignment.*/
        int observedState;
        char thischar;

        if(tupleidx != msa->ss->ntuples)
          thischar = (int)ss_get_char_tuple(msa, tupleidx, sequenceNumber, 0);
        else
          thischar = '-'; /*All gap column case for extended algorithm.*/

        observedState = mod->rate_matrix->inv_states[(int)thischar];

        /*Iterate over all bases setting probability based on base case*/
        for (i = 0; i < alph_size; i++)
          pL[i][n->id] = probabilityOfLeaf(thischar,observedState,i);

      }
      else{/* General recursive case. Calculate probabilities at inner node for all bases.*/

        for (i = 0; i < nstates-1; i++){
          /*Case where this is not the extra (all gaps) column.*/
          if(tupleidx != msa->ss->ntuples)
            /*pL[k][n] :: Probability of base K given, node n.*/
            pL[i][n->id] = probForNodeResidue(i,pL,mod,n,msa,tupleidx,dMatrix);
          else /* Extra (all gap) column case*/
            pL[i][n->id] = probForNodeResidueAllGap(i,pL,mod,n,dMatrix);
        }
        /*Handle gap according to different formula.*/
        if(tupleidx != msa->ss->ntuples)
          pL[4][n->id] = probForNodeGap(pL,mod,n,msa,tupleidx,dMatrix);
        else
          pL[4][n->id] = probForNodeGapAllGap(pL,mod,n,dMatrix);
      }
    }
    
    total_prob = totalProbOfSite(pL, mod->backgd_freqs->data, rootNodeId, p);

    /*Debugging print info:*//*
    int paramIndex = mod->ratematrix_idx;
    double* paramArray = &(mod->all_params->data[paramIndex]);
    printf("Current Rate Matrix values:\n");
    printMatrix(mod->rate_matrix->matrix,5);
    printf("\n\n");

    printf("Current Rates:");
    printf("Lambda: %f\n", paramArray[0]);
    printf("Mu: %f\n", paramArray[1]);
    printf("Alpha: %f\n", paramArray[2]);
    printf("Betta: %f\n", paramArray[3]);
    printf("Geometric Distribution: %f\n",mod->geometricParameter);
    printf("\n\n");
    
    printf("Current background frequencies:\n");
    double* freq = mod->backgd_freqs->data;
    printf("A : %f\n", freq[0]);
    printf("C : %f\n", freq[1]);
    printf("T : %f\n", freq[2]);
    printf("G : %f\n", freq[3]);
    printf("- : %f\n", freq[4]);
    
    printf("\n\n");
    printf("Matrix Information:\n");
    printf("U matrix Information:\n");
    printMatrix(dMatrix->uMatrix,4);
    printf("U Inverse matrix Information:\n");
    printMatrix(dMatrix->uInverseMatrix,4);
    
    printf("Diagonal matrix Information:\n");
    printMatrix(dMatrix->diagonalMatrix,4);
    
    printf("Eigen Values:\n");
    for(i = 0; i < dMatrix->eigenNumber; i ++)
      printf("Eigenvalue %d: %f\n", i, dMatrix->allEigenvalues[i]);
    
    printf("\n\n");
    int j;
    for(i = 0; i < lst_size(traversal); i++)
      for (j = 0; j < 5; j++)
        printf("Likelihood at pL[%d][%d] = %f\n",j,i,pL[j][i]);
    printf("\n\n");*/
    
    /*Multiply by the amount of times this column appeared in the alignment.*/
    if(tupleidx != msa->ss->ntuples)
      retval += log2(total_prob) * msa->ss->counts[tupleidx];
    else /*Case for all gaps column.*/
      allGapProb = total_prob;
  }
    /* Calculate total probability for site in alignment, second modification of
     * extended pruning algorithm.*/

  retval = getProbZeroL(mod, p, allGapProb, retval);

  freeDiagonalMatrix(dMatrix);
  return retval;
}
/* =====================================================================================*/
/**
 * Recusively print tree starting from passed in note.
 * @param node
 */
void printTree(TreeNode* node){
  if(node->lchild == NULL)
    return;
  
  printf("Children of %d are: (node %d) and (node %d)\n", node->id, node->lchild->id, node->rchild->id);
  printf("Left child branch length: %f\n", node->lchild->dparent);
  printf("Right child branch length: %f\n", node->rchild->dparent);
  printTree(node->lchild);
  printTree(node->rchild);
  
  return;
}
/* =====================================================================================*/
/** According to the paper the score must be multiplied by the extra column contribution.
 * Formula (29).
 * @param mod, Tree model for out alignment.
 * @param p, parameter of geometric substitution.
 * @param probAllColumns, probability of all individual columns multiplied
 * @param allGapProb, the un-log(ed) probability of the all gap column.
 * @return total likelihood for whole function.
 */
double getProbZeroL(TreeModel* mod, double p, double allGapProb, double probAllColumns){
  int paramIndex = mod->ratematrix_idx;
  double lambda = mod->all_params->data[paramIndex];
  double mu = mod->all_params->data[paramIndex+1];
  double totalProb = 0;

  double extraColumn = probExtraColumn(mod,mu,lambda,p);
  /*Probability of all-gaps column using equation (25) minus one. Since we are working
   * in log though, we instead negate the value.*/
  double marginalizeValue = log2(1 - allGapProb);
  extraColumn = log2(extraColumn);
  
  totalProb = probAllColumns - marginalizeValue + extraColumn;
  
  return totalProb;
}
/* =====================================================================================*/
/**
 * Given the MSA which already has had it's sufficient statistics computed it will return
 * the parameter of geometric distribution p.
 * @param msa, our MSA for the file.
 * @return parameter of geometric distribution p.
 */
double getGeometricDistribution(MSA* msa){
  double p;
  /*Unique sites in our sufficient statistics.*/
  int uniqueSites = msa->ss->ntuples;
  /*Array holding the amount of times each unique column appears.*/
  double* countArray = msa->ss->counts;
  /*Average length of MSA.*/
  double averageLength = 0;
  int numberOfSpecies = msa->nseqs;
  /*Our sequence.*/
  char** sequence = msa->ss->col_tuples;
  char currentC = 0;
  int sum = 0;
  int totalSum = 0;
  int i,j;

  /*Count the number of non gaps. Notice this array is not how you might expect it to be.
   * It is the transposed version of a typical MSA array. */
  for(i = 0; i < uniqueSites; i++){
    sum = 0;
    for(j = 0; j < numberOfSpecies; j++ ){
      currentC = sequence[i][j];
      if(currentC != '-')
        sum += countArray[i];
    }
    totalSum += sum;
  }
  
  averageLength = (double) totalSum / (double) numberOfSpecies;
  /*This is the way it's done in dnaML, it ignores the rate categories though.*/
  p = averageLength / (averageLength + 1.0);
  
  return p;
}
/* =====================================================================================*/
/**
 * Calculates the extra column probability as described in dnaML-erate paper.
 * This is formula the first part of formula (26), (27), (28).
 * @param mod, model of tree.
 * @param mu, rate of deletion.
 * @param lambda, rate of insertion.
 * @param p, parameter of geometric distribution.
 * @return probability of extra column.
 */
double probExtraColumn(TreeModel* mod,double mu,double lambda, double p){

 TreeNode* root = mod->tree;
 return (1.0 - p) * starProb(root,mu,lambda);
}
/* =====================================================================================*/
/*Given some node, as well as our mu and lambda values. It will recursively traverse the
 * tree finding the starProb as defined by Equation (27) on the paper.
 * @param node, node k to start recursion down on, probably root.
 * @param mu, rate of deletion.
 * @param lambda, rate of insertion.
 * @return starProb value.
 */
double starProb(TreeNode* node,double mu, double lambda){
  /*Reached base case, leaf.*/
  if(node->rchild == NULL)
    return 1;
  
  TreeNode* lChild = node->lchild;
  TreeNode* rChild = node->rchild;
  
  /*Verify this is the correct length (dparent).*/
  double leftProb = starProb(lChild, mu, lambda);
  double rightProb = starProb(rChild, mu, lambda);
  double leftXiPrime = 1.0 - xi(lChild->dparent,mu,lambda);
  double rightXiPrime = 1.0 - xi(rChild->dparent,mu,lambda);

  return leftProb * rightProb * leftXiPrime * rightXiPrime;
}
/* =====================================================================================*/
/**
 * Calculates the probability at site site u of the MSA. This is the actual probability
 * that is added as we calculate the likelihood of the tree. Equation (25).
 * @param pL, array containing the probability pL[k][n], nucleotide k, node k.
 * @param freqs, background frequencies array.
 * @param rootNodeId, id of root node.
 * @param p, parameter of geometric distribution.
 * @return Total probability.
 */
double totalProbOfSite(double** pL,double* freqs,int rootNodeId, double p){
  int gapChar = 4;
  int nucleotides = 4;
  double totalProb = 0.0;
  double sum = 0.0;
  int i;
  
  for(i = 0; i < nucleotides; i++)
    sum += pL[i][rootNodeId] * freqs[i];

  totalProb = pL[gapChar][rootNodeId] + p * sum;
  
  return totalProb;
}
/* =====================================================================================*/
/**
 * Computes the likelihood of node curentNode (k) given the likelihoods up to node k and
 * a current base. Formula (20) on the paper.
 * @param i, given residue to compute for.
 * @param pL, probabilities calculated so far.
 * @param mod, TreeModel containing all information about our model.
 * @param k, node k we're iterating over.
 * @param msa, Multiple sequence alignment of our tree. Needed for deltaChar()
 * @param currSite, u column we are looking at.
 * @return likelihood as computed by formula 20.
 */
double probForNodeResidue(int i,double** pL,TreeModel* mod, TreeNode* k,MSA* msa,
        int currSite,DiagonalMatrix* dMatrix){
  int q;
  int gapChar = 4;
  int residues = 4;

  double* freqs= mod->backgd_freqs->data;
  int paramIndex = mod->ratematrix_idx;
  double* params= &(mod->all_params->data[paramIndex]);

  int lChild = k->lchild->id;
  double lChildDist = k->lchild->dparent;

  int rChild = k->rchild->id;
  double rChildDist = k->rchild->dparent;

  double xRight,yLeft,xLeft,yRight;
  /*Total summations for both sides.*/
  double totalRight = 0;
  double totalLeft = 0;

  /*Left hand summation over all bases.*/
  for(q = 0; q < residues; q++){
    /*Left hand side*/
    xLeft = pL[q][lChild] * singleEventCondProb(dMatrix,q,i,lChildDist,freqs,params);
    yLeft = deltaGap(k->lchild,msa,currSite) *
            singleEventCondProb(dMatrix,gapChar,i,lChildDist,freqs,params);
    /*Right hand side*/
    xRight =  pL[q][rChild] * singleEventCondProb(dMatrix,q,i,rChildDist,freqs,params);
    yRight = deltaGap(k->rchild,msa,currSite) *
            singleEventCondProb(dMatrix,gapChar,i,rChildDist,freqs,params);
    
    totalLeft += xLeft + yLeft;
    totalRight +=  xRight + yRight;
  }

  return totalLeft * totalRight;
}
/* =====================================================================================*/
/**
 * Computes the likelihood of node curentNode (k) given the likelihoods up to node k for
 * gap character. Formula (21) on the paper.
 * @param pL, probabilities calculated so far.
 * @param mod, TreeModel containing all information about our model.
 * @param k, node k we're iterating over.
 * @param msa, Multiple sequence alignment of our tree. Needed for deltaChar()
 * @param currSite, u column we are looking at.
 * @return likelihood as computed by formula 21.
 */
double probForNodeGap(double** pL,TreeModel* mod, TreeNode* k,MSA* msa,int currSite,
        DiagonalMatrix* dMatrix){
  int q;
  int gapChar = 4;
  int residues = 4;

  double* freqs= mod->backgd_freqs->data;
  int paramIndex = mod->ratematrix_idx;
  double* params = &(mod->all_params->data[paramIndex]);

  int lChild = k->lchild->id;
  double lChildDist = k->lchild->dparent;

  int rChild = k->rchild->id;
  double rChildDist = k->rchild->dparent;

  double firstLeft, firstRight, secondLeft, secondRight;
  double totalRight = 0;
  double totalLeft = 0;

  /*Left hand summation over all bases.*/
  for(q = 0; q < residues; q++){
    /*Left hand side*/
    firstLeft =  pL[q][lChild] * singleEventCondProb(dMatrix,q,gapChar,lChildDist,freqs,params);
    secondLeft = pL[gapChar][lChild];
    /*Right hand side*/
    firstRight = pL[q][rChild] * singleEventCondProb(dMatrix,q,gapChar,rChildDist,freqs,params);
    secondRight = pL[gapChar][rChild];
    
    totalLeft += firstLeft + secondLeft;
    totalRight += firstRight + secondRight;
            
  }

  return totalLeft * deltaGap(k->rchild,msa,currSite) +
          totalRight * deltaGap(k->lchild,msa,currSite);
}
/* =====================================================================================*/
/**
 * Same as probForNodeResidue except used for extra column containing all gaps therefore
 * we don not need to check deltaGap() to see if all children are gaps.
 * Computes the likelihood of node curentNode (k) given the likelihoods up to node k and
 * a current base. Formula (20) on the paper.
 * @param i, given residue to compute for.
 * @param pL, probabilities calculated so far.
 * @param mod, TreeModel containing all information about our model.
 * @param k, node k we're iterating over.
 * @param msa, Multiple sequence alignment of our tree. Needed for deltaChar()
 * @return likelihood as computed by formula 20.
 */
double probForNodeResidueAllGap(int i,double** pL,TreeModel* mod, TreeNode* k,
        DiagonalMatrix* dMatrix){
  int q;
  int gapChar = 4;
  int residues = 4;

  double* freqs= mod->backgd_freqs->data;
  int paramIndex = mod->ratematrix_idx;
  double* params= &(mod->all_params->data[paramIndex]);

  int lChild = k->lchild->id;
  double lChildDist = k->lchild->dparent;

  int rChild = k->rchild->id;
  double rChildDist = k->rchild->dparent;

  double xRight,yLeft,xLeft,yRight;
   /*Total summations for both sides.*/
  double totalRight = 0;
  double totalLeft = 0;

  /*Left hand summation over all bases.*/
  for(q = 0; q < residues; q++){
    /*Left hand side*/
    xLeft = pL[q][lChild] * singleEventCondProb(dMatrix,q,i,lChildDist,freqs,params);
    yLeft = singleEventCondProb(dMatrix,gapChar,i,lChildDist,freqs,params);
    /*Right hand side*/
    xRight =  pL[q][rChild] * singleEventCondProb(dMatrix,q,i,rChildDist,freqs,params);
    yRight = singleEventCondProb(dMatrix,gapChar,i,rChildDist,freqs,params);
    
    totalLeft += xLeft + yLeft;
    totalRight +=  xRight + yRight;
  }

  return totalLeft * totalRight;
}
/* =====================================================================================*/
/** Same as probForNodeGap except used for extra column containing all gaps therefore
 * we don not need to check deltaGap() to see if all children are gaps.
 * @param pL, probabilities calculated so far.
 * @param mod, TreeModel containing all information about our model.
 * @param k, node k we're iterating over.
 * @param msa, Multiple sequence alignment of our tree. Needed for deltaChar()
 * @return likelihood as computed by formula 21.
 */
double probForNodeGapAllGap(double** pL,TreeModel* mod, TreeNode* k, DiagonalMatrix* dMatrix){
  int q;
  int gapChar = 4;
  int residues = 4;

  double* freqs= mod->backgd_freqs->data;
  int paramIndex = mod->ratematrix_idx;
  double* params = &(mod->all_params->data[paramIndex]);

  int lChild = k->lchild->id;
  double lChildDist = k->lchild->dparent;

  int rChild = k->rchild->id;
  double rChildDist = k->rchild->dparent;

  double firstLeft, firstRight, secondLeft, secondRight;
  double totalRight = 0;
  double totalLeft = 0;

  /*Left hand summation over all bases.*/
  for(q = 0; q < residues; q++){
    /*Left hand side*/
    firstLeft =  pL[q][lChild] * singleEventCondProb(dMatrix,q,gapChar,lChildDist,freqs,params);
    secondLeft = pL[gapChar][lChild];
    /*Right hand side*/
    firstRight = pL[q][rChild] * singleEventCondProb(dMatrix,q,gapChar,rChildDist,freqs,params);
    secondRight = pL[gapChar][rChild];
    
    totalLeft += firstLeft + secondLeft;
    totalRight += firstRight + secondRight;
  }

  return totalLeft * totalRight;
}
/* =====================================================================================*/
/**
 * The single event conditional probability. Calculates: P(j|i,t) as like formula
 * (22),(23),(24) on dnaML-erate paper. Depending on whether either assumedBase or
 * currentBase is a gap. It uses the gamma-t (6) , xi-t (7) and P_t-epsilon (9)
 * functions to calculate the return value (functions defined below).
 * @param mod, model representing our phylogenetic tree.
 * @param assumedBase, base to assume (j)
 * @param currentBase, given base for probability (i)
 * @param branchLength, given branch length from j to i (t)
 * @param params, array of parameters containing: [mu,lambda,alpha,beta]
 * @return calculated probability.
 */
double singleEventCondProb(DiagonalMatrix* dMatrix,int j,int i, double branchLength,
        double* freqs,double* params){
  int gapCharacter = 4;
  double mu = params[0];
  double lambda = params[1];
  double xiValue = xi(branchLength,mu,lambda);
    
  /*Formula (23)*/
  if(j == gapCharacter){
    return (1-xiValue) * gammaML(branchLength,mu,lambda);
  /*Formula (24)*/
  }else if(i == gapCharacter){
    return xiValue * freqs[j];
  /*Formula (22)*/
  }else{
    return (1-xiValue) * epsilonProbability(dMatrix,j,i,branchLength,freqs,params);
  }
}
/* =====================================================================================*/
/**
 * Calculates the function xi based on branch length, lambda and mu (formula 7)
 * If both lambda and mu are zero then define the function as zero.
 * @param branchLength, given branch length from j to i (t)
 * @param mu, deletion rate.
 * @param lambda, insertion rate.
 * @return value.
 */
double xi(double branchLength,double mu,double lambda){
  if(mu <= 0.0001 && lambda <= 0.0001)
    return 0;
  
  double value;
  double muAndLambda = mu + lambda;
  
  value = lambda / muAndLambda * ( 1 - exp(- muAndLambda * branchLength) );
  
  return value;
}
/* =====================================================================================*/
/**
 * Calculates the function gamma based on branch length, lambda and mu (formula 6)
 * If both lambda and mu are zero then define the function as zero.
 * @param branchLength, given branch length from j to i (t)
 * @param mu, deletion rate.
 * @param lambda, insertion rate.
 * @return value.
 */
double gammaML(double branchLength,double mu,double lambda){
  if(mu <= 0.0001 && lambda <= 0.0001)
    return 0;
  double value;
  double muAndLambda = mu + lambda;

  value = mu / muAndLambda * (1 - exp(-muAndLambda * branchLength));

  return value;
}
/* =====================================================================================*/
/**
 * Calculates formula (9) on dnaML-erate paper. See apendix of dnaML-erate for extra
 * information and derivations of these formulas. Notice if mu and lambda are both zero
 * then model is defined by equation (8).
 * @param mod, model representing our phylogenetic tree.
 * @param i, assumed alphabet letter for probability.
 * @param j,given nucleotide.
 * @param t, branch length.
 * @param freqs, array frequencies.
 * @param params, array of parameters: [mu,lambda,alpha,beta]
 * @return computed formula (9) for M_t(i,j)
 */
double epsilonProbability(DiagonalMatrix* dMatrix,int i,int j,double t,double* freqs,double* params){
  double mu = params[0];
  double lambda = params[1];
  double firstPart = 0;
  double lastPart = 0;
  double eigenSum = 0;
  
  if(lambda != 0 || mu != 0)
    firstPart = freqs[j] * lambda / (lambda + mu);
  else
    firstPart = freqs[j];

  if(lambda != 0 || mu != 0)
    lastPart = freqs[j] * mu / (lambda + mu) * exp(-(lambda + mu) * t);
  else
    lastPart = 0;

  eigenSum = computeOhMatrixSummation(i,j,dMatrix,mu,t);

  return firstPart + eigenSum + lastPart;
}
/* =====================================================================================*/
/**
 * Given the parameters will compute the summation over all eigenvalues for our R matrix,
 * sum = 0
 * for eigen in eigenValues:
 *   sum += O(a,i,j)*exp(-(eigen+mu)*t)
 * Where Matrix O is defined in the apendix of the paper.
 * @param rMatrix, rateMatrix for our model.
 * @param mu, rate of deletions.
 * @param t, branch length.
 * @return sum.
 */
double computeOhMatrixSummation(int i, int j, DiagonalMatrix* dMatrix, double mu, double t){
  double currentEigen;
  double sum = 0;
  int k, l;
  int position;
  double oResults;
  
  /*Get from Matrix.*/
  int eigenNumber = dMatrix->eigenNumber;
  double* allEigenvalues = dMatrix->allEigenvalues;
  Matrix* diagonalMatrix = dMatrix->diagonalMatrix;
  Matrix* uMatrix = dMatrix->uMatrix;
  Matrix* uInverseMatrix = dMatrix->uInverseMatrix;
  
  /*If we have multiplicity higher than 1 we must account for this when finding the 
   * position of our eigenvalues in the diagonal matrix.*/
  int multiplicity[4] = {0, 0, 0, 0};
  
  for( k = 0; k < eigenNumber; k ++){
    currentEigen = allEigenvalues[k];
    for( l = 0; l <= k; l ++)
      if (abs(currentEigen - allEigenvalues[l]) < 0.00001)
        multiplicity[k]++;
  }

  /*Iterate through eigenvalues adding up values.*/
  for(k = 0; k < eigenNumber; k ++){
    currentEigen = allEigenvalues[k];
    position = findPositionOfEigenvalue(diagonalMatrix, currentEigen,multiplicity[k]);
    /* O_a(i,j) = U(i,pos)*U^(-1)(pos,j)*/
    oResults = mat_get(uMatrix, i, position) * mat_get(uInverseMatrix,position,j);
    /*Notice subtraction due to distributing the negative sign. Eigenvalues are supposed
     * to be turned positive by negating (they should all come out negative, except zero).
     */
    sum += oResults * exp((currentEigen - mu) * t);
  }

  return sum;
}
/* =====================================================================================*/
/**
 * Helper function to print matrices.
 * @param m, the actual matrix.
 * @param size, size of matrix.
 */
void printMatrix(Matrix* m, int size){
  int k,l;

  for(k=0;k < size; k++){
    for(l = 0;l < size; l++)
      printf("%f ",m->data[k][l]);
    printf("\n");
  }
  printf("\n");
  return;
}
/* =====================================================================================*/
/**
 * Get position of our matrix in our diagnolized matrix (k,k)
 * @param matrix
 * @param eigenval
 * @param mult, the multiplicity of this eigenvalue.
 * @return position on success, -1 on failure.
 */
int findPositionOfEigenvalue(Matrix* matrix,double eigenval, int mult){
  int l;
  int position = -1;
  int size = 4;

  for(l = 0; l < size;l++)
    
    /*Difference smaller than delta to account for round off error.*/
    if(fabs(mat_get(matrix,l,l) - eigenval) < 0.001){
      position = l;
      
      if (mult == 1)
        break;
      else
        mult--;
    }
  /*Position not found die!*/
  if(position == -1){
    printf("Eigen vector: %f not found in diagonal Matrix, sorry...", eigenval);
    exit(1);
  }
  return position;
}
/* =====================================================================================*/
/**
 * Recursive function. Given a node on a tree will return either true if node was
 * a gap or false if node was not a gap. Used by probForNodeResidue and probForNodeGap.
 * @param k, node being looked at.
 * @param msa, we need the MSA to see what the characters look like at tree.
 * @return predicate stating whether all leaves from k down were gaps or not.
 */
int deltaGap(TreeNode* k,MSA* msa,int currSite){

  /*Base case we hit the leaf. Notice assumes properly made tree*/
  if(k->lchild == NULL && k->rchild == NULL){
    char* name = k->name;
    /*Look up the characters at these nodes based on the name.*/
    char base = getCharacterForSpecie(name,msa,currSite);
    
    return isGap(base);
  }
  
  /*Else recurse down and return the results of our recursion.*/
  TreeNode* leftChild = k->lchild;
  TreeNode* rightChild = k->rchild;
  
  return deltaGap(leftChild,msa,currSite) * deltaGap(rightChild,msa,currSite);
}
/* =====================================================================================*/
/**
 * Returns whether the passed char was a gap ;)
 * @param c, character
 * @return 1 for true, 0 for false.
 */
int isGap(char c){
  if(c == '-')
    return 1;
  else
    return 0;
}
/* =====================================================================================*/
/**
 * Given a species name, a multiple sequence alignment and the index of the column where
 * you want the character from. Returns the character belonging to that species at that
 * position.
 * @param name, name of the specie.
 * @param msa, multiple sequence alignment.
 * @param index, site in the MSA to get char from.
 * @return char at that position.
 */
char getCharacterForSpecie(char* name,MSA* msa,int index){
  char myChar;
  char** nameTable = msa->names;
  int size = msa->nseqs;
  /*Index is sequence for my species.*/
  int specieIndex = -1;
  int i;
  /*Iterate through names until we find the index */
  for(i = 0; i < size; i++)
    if(strcmp(nameTable[i],name) == 0){
      specieIndex = i;
      break;
    }
  /*Index now holds the location of our species in the SS.*/
  myChar = msa->ss->col_tuples[index][specieIndex];

  return myChar;
}
/* =====================================================================================*/
/* P_u(L_k,i) = { 1, if leak k has residue i at position u;
 *              { 0, otherwise;
 * For gaps: P_u(L_k,'-') = 0  // Equation 18 on paper.
 */
int probabilityOfLeaf(char currentChar,int observedState,int iResidue){
  /*Gap case:*/
  if (currentChar == '-')
    return 0;
  /*Case where residue matches what we see.*/
  else if (observedState < 0 || iResidue == observedState) /*When would this ever be less than zero??*/
    return 1;
  else
    return 0;
}
/* =====================================================================================*/
/*EOO (End of Omar)*/
/* this is retained for possible use in the future; not using weight
   matrices for much anymore */
void tl_compute_log_likelihood_weight_matrix(TreeModel *mod, MSA *msa, 
                                             double *col_scores, int cat) {
  int i, seq, idx, alph_size = strlen(msa->alphabet);
  double retval = 0;
  char tuple[mod->order + 2];
  Vector *margfreqs = 
    get_marginal_eq_freqs(mod->rate_matrix->states, mod->order+1,
                          mod->backgd_freqs);
  int col_by_col = (col_scores != NULL || msa->ss == NULL);
                                /* evaluate the alignment
                                   column-by-column if either
                                   col-by-col scores are required or
                                   the sufficient stats are not
                                   available */
                                /* NOTE: !col_by_col -> msa->ss != NULL */

  checkInterrupt();
  tuple[mod->order+1] = '\0';

  if (mod->tree != NULL)
    die("ERROR tl_compute_log_likelihood_weight_matrix: mod->tree should be NULL\n");
  if (msa->ss == NULL && msa->seqs == NULL)
    die("ERROR tl_compute_log_likelihood_weight_matrix: mod->ss and mod->seqs are both NULL\n");

  if (cat >= 0) {
    if (col_by_col) {
      if (msa->categories == NULL)
	die("ERROR tl_compute_log_likelihood_weight_matrix: msa->categories is NULL\n");
    }
    /* if using categories and col-by-col
       scoring, then must have col-by-col
       categories */
    else if (msa->ss->cat_counts == NULL)
      die("ERROR tl_compute_log_likelihood_weight_matrix: msa->ss->cat_counts is NULL\n");
    /* if using categories and unordered
       sufficient statistics, must have
       category-by-category counts */
  }

  if (col_by_col)
    if (msa->seqs == NULL && msa->ss->tuple_idx == NULL)
      die("ERROR tl_compute_log_likelihood requires ordered alignment\n");
  /* if using col-by-col scoring, must
     have ordered representation */
    
  retval = 0;
  for (idx = 0; 
       idx < (col_by_col ? msa->length : msa->ss->ntuples);
       idx++) {
    int thisstate, col, tupleidx;
    double col_val = 0, prob = 0;

    /* NOTE: when evaluating col-by-col, idx is a column, but otherwise
       idx is a tuple index.  Let's be clear about this ... */
    col = (col_by_col ? idx : -1); 
    if (msa->ss == NULL) tupleidx = -1; 
    else if (col_by_col) tupleidx = msa->ss->tuple_idx[col];
    else tupleidx = idx;      /* NOTE: tupleidx will be defined
                                 whenever msa->ss != NULL */

    if (cat < 0 || (col_by_col && msa->categories[col] == cat) ||
        (!col_by_col && msa->ss->cat_counts[cat][tupleidx] > 0)) {

      for (seq = 0; seq < msa->nseqs; seq++) {

        for (i = -mod->order; i <= 0; i++) {
          if (msa->ss != NULL)
            tuple[mod->order+i] = ss_get_char_tuple(msa, tupleidx, seq, i); 
          else if (col + i >= 0)
            tuple[mod->order+i] = msa->seqs[seq][col+i];
          else 
            tuple[mod->order+i] = msa->missing[0];
        }

        if (!mod->allow_gaps && 
            msa->inv_alphabet[(int)tuple[mod->order]] < 0 &&
            !msa->is_missing[(int)tuple[mod->order]]) {
            
          col_val = NEGINFTY; /* we want to apply the strict penalty
                                 iff there is an unrecognized
                                 character in *this* (the rightmost)
                                 col; missing data is a special case -- don't
                                 penalize even if not in alphabet */

                                /* FIXME: seems too complicated --
                                   just check for gap? */
          break; 
        }
        else if (mod->allow_but_penalize_gaps &&
                 msa->inv_alphabet[(int)tuple[mod->order]] < 0) { 
          /* temporary */
          double tmp_prob;
          prob = 1;
          for (i = 0; i < alph_size; i++) {
            tuple[mod->order] = msa->alphabet[i];
            tmp_prob = vec_get(margfreqs, tuple_index_missing_data(tuple, msa->inv_alphabet, msa->is_missing, alph_size));
            if (tmp_prob < prob) prob = tmp_prob;
          }
          if (prob == 0) prob = 0.01;
        }
        else {
          thisstate = tuple_index_missing_data(tuple, msa->inv_alphabet, 
                                               msa->is_missing, alph_size);
          prob = vec_get(margfreqs, thisstate);
        }

        if (prob == 0) { col_val = NEGINFTY; break; }

        col_val += log2(prob);

        if (mod->use_conditionals && mod->order > 0) {
          tuple[mod->order] = msa->missing[0];
          col_val -= log2(vec_get(margfreqs, tuple_index_missing_data(tuple, msa->inv_alphabet, msa->is_missing, alph_size)));
        }
      }
    }
    if (!col_by_col)   /* tuple-by-tuple scoring */
      col_val *= (cat >= 0 ? msa->ss->cat_counts[cat][tupleidx] : 
                  msa->ss->counts[tupleidx]);
    retval += col_val;
    if (col_scores != NULL) col_scores[col] = col_val;
  }
  if (retval < NEGINFTY) retval = NEGINFTY; 
  /* must be true if any of the columns
     considered had prob NEGINFTY */
  vec_free(margfreqs);
}


TreePosteriors *tl_new_tree_posteriors(TreeModel *mod, MSA *msa, int do_bases, 
                                       int do_substs, int do_expected_nsubst, 
                                       int do_expected_nsubst_tot,
				       int do_expected_nsubst_col,
                                       int do_rate_cats, int do_rate_cats_exp) {
  int i, j, k, r, ntuples, nnodes, nstates;
  TreePosteriors *tp = (TreePosteriors*)smalloc(sizeof(TreePosteriors));

  if (mod->tree ==  NULL)
    die("ERROR tl_new_tree_posteriors: mod->tree is NULL\n");
  if (msa->ss == NULL)
    die("ERROR tl_new_tree_posteriors: msa->ss is NULL\n");

  ntuples = msa->ss->ntuples;
  nnodes = mod->tree->nnodes;
  nstates = mod->rate_matrix->size;

  if (do_bases) {
    tp->base_probs = (double****)smalloc(mod->nratecats * sizeof(double***));
    for (r = 0; r < mod->nratecats; r++) {
      tp->base_probs[r] = (double***)smalloc(nstates * sizeof(double**));
      for (i = 0; i < nstates; i++) {
        tp->base_probs[r][i] = (double**)smalloc(nnodes * sizeof(double*));
        for (j = 0; j < nnodes; j++) {
          tp->base_probs[r][i][j] = (double*)smalloc(ntuples * sizeof(double));
        }
      }
    }
  }
  else tp->base_probs = NULL;

  if (do_substs) {
    tp->subst_probs = (double*****)smalloc(mod->nratecats * sizeof(double****));
    for (r = 0; r < mod->nratecats; r++) {
      checkInterrupt();
      tp->subst_probs[r] = (double****)smalloc(nstates * sizeof(double***));
      for (i = 0; i < nstates; i++) {
        tp->subst_probs[r][i] = (double***)smalloc(nstates * sizeof(double**));
        for (j = 0; j < nstates; j++) {
          tp->subst_probs[r][i][j] = (double**)smalloc(nnodes * sizeof(double*));
          for (k = 0; k < nnodes; k++) 
            if (k != mod->tree->id) /* don't need one for the root */
              tp->subst_probs[r][i][j][k] = (double*)smalloc(ntuples * 
                                                            sizeof(double));
            else 
              tp->subst_probs[r][i][j][k] = NULL;
        }
      }
    }
  }
  else tp->subst_probs = NULL;

  if (do_expected_nsubst) {
    tp->expected_nsubst = (double***)smalloc(mod->nratecats * sizeof(double**));
    for (r = 0; r < mod->nratecats; r++) {
      tp->expected_nsubst[r] = (double**)smalloc(nnodes * sizeof(double*));
      for (i = 0; i < nnodes; i++) {
        if (i != mod->tree->id)
          tp->expected_nsubst[r][i] = (double*)smalloc(ntuples * sizeof(double));
        else
          tp->expected_nsubst[r][i] = NULL;
      }
    }
  }
  else tp->expected_nsubst = NULL;

  if (do_expected_nsubst_tot) {
    tp->expected_nsubst_tot = (double****)smalloc(mod->nratecats * sizeof(double***));
    for (r = 0; r < mod->nratecats; r++) {
      tp->expected_nsubst_tot[r] = (double***)smalloc(nstates * sizeof(double**));
      for (i = 0; i < nstates; i++) {
        tp->expected_nsubst_tot[r][i] = (double**)smalloc(nstates * 
                                                         sizeof(double*));
        for (j = 0; j < nstates; j++) 
          tp->expected_nsubst_tot[r][i][j] = (double*)smalloc(nnodes * 
                                                             sizeof(double));
      }
    }
  }
  else tp->expected_nsubst_tot = NULL;

  if (do_expected_nsubst_col) {
    tp->expected_nsubst_col = (double*****)smalloc(mod->nratecats * sizeof(double****));
    for (r=0; r < mod->nratecats; r++) {
      tp->expected_nsubst_col[r] = (double****)smalloc(nnodes * sizeof(double***));
      for (i=0; i < nnodes; i++) {
	tp->expected_nsubst_col[r][i] = (double***)smalloc(ntuples * sizeof(double**));
	for (j=0; j < ntuples; j++) {
	  tp->expected_nsubst_col[r][i][j] = (double**)smalloc(nstates * sizeof(double*));
	  for (k=0; k < nstates; k++) 
	    tp->expected_nsubst_col[r][i][j][k] = (double*)smalloc(nstates * sizeof(double));
	}
      }
    }
  }
  else tp->expected_nsubst_col = NULL;

  if (do_rate_cats) {
    tp->rcat_probs = (double**)smalloc(mod->nratecats * sizeof(double*));
    for (i = 0; i < mod->nratecats; i++)
      tp->rcat_probs[i] = (double*)smalloc(ntuples * sizeof(double));
  }
  else tp->rcat_probs = NULL;

  if (do_rate_cats_exp) 
    tp->rcat_expected_nsites = (double*)smalloc(mod->nratecats * sizeof(double));
  else tp->rcat_expected_nsites = NULL;

  return tp;
}

void tl_free_tree_posteriors(TreeModel *mod, MSA *msa, TreePosteriors *tp) {
  int i, j, k, r, ntuples, nnodes, nstates;

  if (mod->tree == NULL)
    die("ERROR tl_free_tree_posteriors: mod->tree is NULL\n");
  if (msa->ss == NULL)
    die("ERROR tl_free_tree_posteriors: msa->ss is NULL\n");
  ntuples = msa->ss->ntuples;
  nnodes = mod->tree->nnodes;
  nstates = mod->rate_matrix->size;

  if (tp->base_probs != NULL) {
    for (r = 0; r < mod->nratecats; r++) {
      for (i = 0; i < nstates; i++) {
        for (j = 0; j < nnodes; j++) 
          if (tp->base_probs[r][i][j] != NULL)
            sfree(tp->base_probs[r][i][j]);
        sfree(tp->base_probs[r][i]);
      }
      sfree(tp->base_probs[r]);
    }
    sfree(tp->base_probs);
  }
  if (tp->subst_probs != NULL) {
    for (r = 0; r < mod->nratecats; r++) {
      for (i = 0; i < nstates; i++) {
        for (j = 0; j < nstates; j++) {
          for (k = 0; k < nnodes; k++) 
            if (k != mod->tree->id) 
              sfree(tp->subst_probs[r][i][j][k]);
          sfree(tp->subst_probs[r][i][j]);
        }
        sfree(tp->subst_probs[r][i]);
      }
      sfree(tp->subst_probs[r]);
    }
    sfree(tp->subst_probs);
  }
  if (tp->expected_nsubst != NULL) {
    for (r = 0; r < mod->nratecats; r++) {
      for (i = 0; i < nnodes; i++) 
        if (i != mod->tree->id)
          sfree(tp->expected_nsubst[r][i]);
      sfree(tp->expected_nsubst[r]);
    }
    sfree(tp->expected_nsubst);
  }
  if (tp->expected_nsubst_tot != NULL) {
    for (r = 0; r < mod->nratecats; r++) {
      for (i = 0; i < nstates; i++) {
        for (j = 0; j < nstates; j++) 
          sfree(tp->expected_nsubst_tot[r][i][j]);
        sfree(tp->expected_nsubst_tot[r][i]);
      }
      sfree(tp->expected_nsubst_tot[r]);
    }
    sfree(tp->expected_nsubst_tot);
  }
  if (tp->expected_nsubst_col != NULL) {
    for (r = 0; r < mod->nratecats; r++) {
      for (i=0; i < nnodes; i++) {
	for (j=0; j < ntuples; j++) {
	  for (k=0; k < nstates; k++) 
	    sfree(tp->expected_nsubst_col[r][i][j][k]);
	  sfree(tp->expected_nsubst_col[r][i][j]);
	}
	sfree(tp->expected_nsubst_col[r][i]);
      }
      sfree(tp->expected_nsubst_col[r]);
    }
    sfree(tp->expected_nsubst_col);
  }
  if (tp->rcat_probs != NULL) {
    for (i = 0; i < mod->nratecats; i++)
      sfree(tp->rcat_probs[i]);
    sfree(tp->rcat_probs);           
  }
  if (tp->rcat_expected_nsites != NULL) {
    sfree(tp->rcat_expected_nsites);           
  }

  sfree(tp);
}

/* compute the expected (posterior) complete log likelihood of a tree
   model based on a TreePosteriors object.  Equilibrium frequencies
   are not considered. */
double tl_compute_partial_ll_suff_stats(TreeModel *mod, TreePosteriors *post) {
  double retval = 0;
  int i, j, k, cat;
  TreeNode *n;
  int nstates = mod->rate_matrix->size;

  for (cat = 0; cat < mod->nratecats; cat++) {
    for (i = 0; i < mod->tree->nnodes; i++) {
      MarkovMatrix *subst_mat;
      if (i == mod->tree->id) continue; /* skip root */
      n = lst_get_ptr(mod->tree->nodes, i);
      subst_mat = mod->P[n->id][cat];
      for (j = 0; j < nstates; j++) { /* from tuple */
        for (k = 0; k < nstates; k++) { /* to tuple */
          retval += (post->expected_nsubst_tot[cat][j][k][i] *
                     log2(mm_get(subst_mat, j, k)));
        }
      }
    }
  }
  return retval;
}

/* The functions below are used for computing likelihoods with
   weight-matrix models.  They should possibly be moved to misc.c.  */

/* given an alphabet, a tuple size, and a vector of equilibrium
   frequences, create a new vector of marginal equilibrium
   frequencies describing the space of "meta-tuples", which contain
   actual characters *or* missing data characters.  Each meta-tuple is
   given an equilibrium frequency equal to the sum of the frequencies
   of all "matching" ordinary tuples.  Missing data characters are
   assumed to be gap characters or Ns. */
Vector *get_marginal_eq_freqs (char *alphabet, int tuple_size, 
			       Vector *eq_freqs) {
  int alph_size = strlen(alphabet);
  int ntuples = int_pow(alph_size, tuple_size);
  int i;
  Vector *retval = vec_new(int_pow(alph_size+1, tuple_size));
  vec_zero(retval);

  /* loop through the ordinary (non-meta) tuples */
  for (i = 0; i < ntuples; i++) {
    int digits[tuple_size];
    int j, k, remainder;
    
    /* first decompose the tuple into its "digits" */
    remainder = i;
    for (j = 0; j < tuple_size; j++) { /* from least sig. to most */
      digits[j] = remainder % alph_size;
      remainder /= alph_size;
    }

    /* now consider every pattern of missing-data characters that can
       be overlaid on it.  The equilibrium frequency of the tuple
       contributes to the marginal frequency corresponding to every
       such pattern.  There are 2^tuple_size of them to consider. */
    for (k = 0; k < (1 << tuple_size); k++) {
      int newtuple = 0, base = 1;
      for (j = 0; j < tuple_size; j++) {
        if (k & (1 << j)) 
          newtuple += alph_size * base;
        else 
          newtuple += digits[j] * base;
        base *= (alph_size + 1);
      }
      vec_set(retval, newtuple, 
                     vec_get(retval, newtuple) + 
                     vec_get(eq_freqs, i));
    }
  }
  return retval;
}

/* given a tuple consisting of actual characters and/or missing data,
   return the corresponding state number in the "meta-tuple" space.
   Returns -1 for unallowed tuples */
int tuple_index_missing_data(char *tuple, int *inv_alph, int *is_missing,
                             int alph_size) {
  int retval = 0, i;
  int tuple_size = strlen(tuple);
  for (i = 0; i < tuple_size; i++) {
    int charidx = inv_alph[(int)tuple[tuple_size-i-1]];
    if (charidx < 0) {
      if (tuple[tuple_size-i-1] == GAP_CHAR || 
          is_missing[(int)tuple[tuple_size-i-1]])
        charidx = alph_size;
      else return -1;
    }
    retval += charidx * int_pow(alph_size+1, i);
                                /* i == 0 => least sig. dig; 
                                   i == tuple_size-1 => most sig. dig */
  }
  return retval;
}

