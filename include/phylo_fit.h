/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: phylo_fit.h,v 1.17 2008-11-12 02:07:59 acs Exp $ */

#ifndef PHYLO_FIT_H
#define PHYLO_FIT_H


/* default minimum number of informative sites (see -I) */
#define DEFAULT_NSITES_THRESHOLD 50

/* default starting alpha for dgamma */
#define DEFAULT_ALPHA 1

#include <tree_model.h>
#include <msa.h>
#include <category_map.h>
#include <stringsplus.h>
#include <lists.h>
#include <gff.h>
#include <list_of_lists.h>


struct phyloFit_struct {
  MSA *msa;
  List *msas;
  List *msa_file_names;
  char *output_fname_root, 
    *reverse_group_tag,
    *root_seqname, *subtree_name, *error_fname,
    *see_for_help, *parsimony_cost_fname,
    *msa_fname;  //note: msa_fname will not be read by run_phyloFit.  
                 //The msa should already be loaded.  msa_fname can be
                 //optionally set for producing informative messages
  int subst_mod, quiet, nratecats, use_em, 
    window_size, window_shift, use_conditionals,
    precision, likelihood_only, do_bases, 
    do_expected_nsubst, do_expected_nsubst_tot, 
    do_expected_nsubst_col,
    random_init, estimate_backgd, estimate_scale_only,
    do_column_probs, nonoverlapping, gaps_as_bases,
    no_freqs, no_rates, assume_clock, 
    init_parsimony, parsimony_only, no_branchlens,
    label_categories, symfreq, init_backgd_from_data,
    use_selection, max_em_its;
  unsigned int nsites_threshold;
  TreeNode *tree;
  List *trees;
  CategoryMap *cm;
  String *nooptstr;
  List *cats_to_do_str,  *window_coords, 
    *ignore_branches, *alt_mod_str,
    *bound_arg, *rate_consts, 
    *label_str, *label_type;
  double alpha, selection;
  GFF_Set *gff;
  TreeModel *input_mod;
  List *input_mods;
  List *mod_file_names;
  FILE *logf;
  int dnaMlNormalize;              /**<Flag leting us know user has selected tree to be computed
                               * how dnaML does. That is, all branchlengths are divided by
                               * the fracChange and tree is midpoint rooted. */
  int reroot;                      /**<reroot tree if user choses for it to do so. */
  //results go here if not-null
  ListOfLists *results;
};


struct phyloFit_struct* phyloFit_struct_new();
int run_phyloFit(struct phyloFit_struct*);
int run_phyloFit_multi(struct phyloFit_struct *pf);
List* setupmsa(struct phyloFit_struct *pf, MSA* msa);
List* setupmod(struct phyloFit_struct *pf, TreeModel *input_mod, MSA *msa, TreeNode *tree, List* cats_to_do, List* newparams,int toSkip,List* toSkipList);
void printmodMulti(struct phyloFit_struct *pf, TreeModel* mod, MSA* msa, List* cats_to_do,int index);
void print_post_prob_stats(TreeModel *mod, MSA *msa, char *output_fname_root, 
                           int do_bases, int do_expected_nsubst, 
                           int do_expected_nsubst_tot, 
			   int do_expected_nsubst_col,
			   int do_every_site,
			   int cat, int quiet,
			   ListOfLists *results);


/**
* Printing functin printing additional information computed from running phyloFit
 * with the extended Prunning Algorithm.
 * @param fileName, name of file to write to.
 * @param tm, model containing information we shall write.
 */
void printExtendedInfo(char* fileName, TreeModel* tm);

/**
 * Give the matrix rate parameters and the background frequencies including gaps return
 * the fractional rate of change as done by dnaMl.
 * @param lambda
 * @param mu
 * @param alpha
 * @param betta
 * @param freqA
 * @param freqC
 * @param freqG
 * @param freqT
 * @param freqGap
 * @return 
 */
double computeFracChange(double lambda, double mu, double alpha, double betta,
        double freqA, double freqC, double freqG, double freqT, double freqGap);

#define BRANCH_TYPE 0
#define SUBTREE_TYPE 1

#endif
