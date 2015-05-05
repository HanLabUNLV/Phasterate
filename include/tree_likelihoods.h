/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file tree_likelihoods.h
    Computation of likelihoods for columns of a given multiple
    alignment, according to a given tree model.
    @ingroup phylo
 */

#ifndef TREE_LIK_H
#define TREE_LIK_H

#include <tree_model.h>
#include <msa.h>
#include <math.h>
#include <misc.h>

/** Structure for information related to posterior probability of tree
   model wrt an alignment.  
     Each array is indexed as appropriate for
   rate categories, bases in the alignment (indexed according to a
   model's inv_states; with higher order models, actually tuples of
   bases), nodes or edges in the tree (indexed by node ids; a node is
   associated with the edge that connects it to its parent), and
   column tuples in a "sufficient statistics" representation of the
   alignment (all quantities will be the same for all instances of a
   column tuple). */
struct tp_struct {
  double ****base_probs;        /**< Posterior probability of each base
                                   given a node, a column
                                   tuple, and a rate category.  
					- First index is rate category
					- Second is base
					- Third is node
					- Fourth is column tuple
				 */
  double *****subst_probs;      /**< Posterior probability of a
                                   substitution of each base for each
                                   other, given a branch, column
                                   tuple, and rate category.  
					- First index is rate category
					- Second is original base
					- Third is replacement base
					- Fourth is branch
					- Fifth is column tuple 
				*/
  double ***expected_nsubst;    /**< Expected number of substitutions for each
                                   branch x column tuple, given a rate
                                   category (conditioned on rate
                                   category in case posterior
                                   probabilities of rate categories
                                   depend on an HMM or similar).
                                   	- First index is rate category
					- Second is branch 
					- Third is column tuple 
				*/ 
  double ****expected_nsubst_tot; 
                                /**< Total expected number of
                                   substitutions of each type along
                                   each branch for each rate category,
                                   summed over all column tuples
                                   (considering the number of
                                   instances of each tuple).  These
                                   are the sufficient statistics for
                                   computing the likelihood of a tree
                                   model.  Note that they are based on
                                   *joint* probabilities with rate
                                   categories, rather than being
                                   conditioned on rate categories (the
                                   posterior probability of each rate
                                   at each site in incorporated).
                                   	- First index is rate category
					- Second is original base
					- Third is replacement base
					- Fourth is branch 
 */
  double *****expected_nsubst_col;
                                /**< Expected number of substitutions of each
                                   type along each branch for each rate 
				   category, for each tuple column.
				   	- First index is rate category
					- Second is branch 
					- Third is tuple
					- Fourth is original base
					- Fifth is replacement base 
				*/
  double **rcat_probs;          /**< Posterior probability of each rate
                                   category for each column tuple.
                                    	- First index is rate category
					- Second is column tuple 
				*/
  double *rcat_expected_nsites; /**< Expected number of sites in each
                                   rate category */
};

typedef struct tp_struct TreePosteriors;
                                /* see incomplete type in tree_model.h */

#define NULL_LOG_LIKELIHOOD 1   /** Safe value for null when dealing with
                                   log likelihoods (should always be <= 0) FIXME? */

/* does not appear to be implemented */
void tl_dump_matrices(TreeModel *mod, double **inside_vals, 
                      double **outside_vals, double **posterior_probs);

/** Compute the likelihood of a tree model with respect to an
   alignment; Optionally retain column-by-column likelihoods and/or posterior probabilities.  
   @param[in] mod Tree Model to compute likelihood for
   @param[in] msa Multiple Alignment containing data related to tree model
   @param[out] col_scores (Optional) Log likelihood score per column
   @param[out] tuple_scores (Optional) Log likelihood score per tuple
   @param[in] cat Whether to use categories
   @param[out] post (Optional) Computed posterior probabilities; If NULL, no
   posterior probabilities (or related quantities) will be computed.
   If non-NULL each of its attributes must either be NULL or
   previously allocated to the required size. 
   @result Log likelihood of entire tree model specified
*/
double tl_compute_log_likelihood(TreeModel *mod, MSA *msa, 
                                 double *col_scores, 
				 double *tuple_scores, 
				 int cat,
                                 TreePosteriors *post);

/** Create a new TreePosteriors object.
    @param mod Tree Model of which the posterior probabilities are calculated
    @param msa Multiple Alignment
    @param do_bases Whether to allocate space for base probabilities
    @param do_subst Whether to allocate space for substitution probabilities
    @param do_expected_nsubst Whether to allocate space for expected number of substitutions matrix
    @param do_expected_nsubst_tot Whether to allocate space for total expected number of substitutions
    @param do_expected_nsubst_col Whether to allocate space for expected number of substitutions per column
    @param do_rate_cats Whether to allocate space for rate categories
    @param do_rate_cats_exp Whether to allocate space for expected rate categories
    @result Newly allocated TreePosteriors object
*/
TreePosteriors *tl_new_tree_posteriors(TreeModel *mod, MSA *msa, int do_bases, 
                                       int do_substs, int do_expected_nsubst, 
                                       int do_expected_nsubst_tot,
				       int do_expected_nsubst_col,
                                       int do_rate_cats, int do_rate_cats_exp);

/** Free TreePosteriors object
   @param mod Tree model of which posterior are calculated
   @param msa Multiple Alignment
   @param tp TreePosteriors object to free
 */
void tl_free_tree_posteriors(TreeModel *mod, MSA *msa, TreePosteriors *tp);

/** Compute the expected (posterior) complete log likelihood of a tree
   model based on a TreePosteriors object.  
   @param[in] mod Tree Model
   @param[in] post Pre-calculated posterior probabilities
   @note Equilibrium frequencies are not considered
   @result Log Likelihood of tree
*/
double tl_compute_partial_ll_suff_stats(TreeModel *mod, TreePosteriors *post);

/* Could not find implementation */
double tl_compute_ll_suff_stats(TreeModel *mod, MSA *msa, TreePosteriors *post);

/** Given an alphabet, a tuple size, and a vector of equilibrium
   frequencies, create a new vector of marginal equilibrium
   frequencies describing the space of "meta-tuples", which contain
   actual characters *or* missing data characters.  
   Each meta-tuple is
   given an equilibrium frequency equal to the sum of the frequencies
   of all "matching" ordinary tuples.  
    Missing data characters are
   assumed to be gap characters or Ns. 
   @param alphabet List of possible characters
   @param tuple_size Size of tuples
   @param eq_freqs Equilibrium frequencies
   @param New vector of marginal equilibrium frequencies
*/
Vector *get_marginal_eq_freqs (char *alphabet, int tuple_size, 
                                   Vector *eq_freqs);

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
   @result Log likelihood of entire tree model specified
 *  */
double gapAwareLikelihood(TreeModel *mod, MSA *msa,double *col_scores, double *tuple_scores,
				 int cat, TreePosteriors *post);

/*Given a TreeModel, an MSA, the cat number, Tree Posteriors. Will return the total
 * likelihood for a tree. The rest of parameters are allocated memory. See gapAwareLikelihood
 * for allocations. This is function does the bulk of the work.
 * Note this is only guaranteed to work and tested on models of order zero with no column
 * offset and no rate categories set. */
  
double computeTotalTreeLikelihood(TreeModel* mod,MSA* msa,int cat,TreePosteriors* post,
        double **inside_joint);

/** According to the paper the score must be multiplied by the extra column contribution.
 * 
 * @param mod, Tree model for out alignment.
 * @param p, parameter of geometric substitution.
 * @param probAllColumns, probability of all individual columns multiplied.
 * @return total likelihood for whole function.
 */
double getProbZeroL(TreeModel* mod,double p, double allgapProb, double probAllColumns);

/**
 * Given the MSA which already has had it's sufficient statistics computed it will return
 * the parameter of geometric distribution p.
 * @param msa, our MSA for the file.
 * @return parameter of geometric distribution p.
 */
double getGeometricDistribution(MSA* msa);

/**
 * Calculates the extra column probability as described in dnaML-erate paper.
 * This is formula the first part of formula (26).
 * @param mod, model of tree.
 * @param mu, rate of deletion.
 * @param lambda, rate of insertion.
 * @param p, parameter of geometric distribution.
 * @return probability of extra column.
 */
double probExtraColumn(TreeModel* mod,double mu,double lambda, double p);

/*Given some node, as well as our mu and lambda values. It will recursively traverse the
 * tree finding the starProb as defined by Equation (27) on the paper.
 * @param node, node k to start recursion down on, probably root.
 * @param mu, rate of deletion.
 * @param lambda, rate of insertion.
 * @param mod, the tree model representing our data.
 * @return starProb value.
 */
double starProb(TreeNode* node,double mu, double lambda, TreeModel* mod);

/**
 * Calculates the probability at site site u of the MSA. This is the actual probability
 * that is added as we calculate the likelihood of the tree.
 * @param pL, array containing the probability pL[k][n], nucleotide k, node k.
 * @param freqs, background frequencies array.
 * @param rootNodeId, id of root node.
 * @param p, parameter of geometric distribution.
 * @return Total probability.
 */
double totalProbOfSite(double** pL,double* freqs,int rootNodeId, double p);

/**
 * Computes the likelihood of node curentNode (k) given the likelihoods up to node k and
 * a current base. Formula (20) on the paper.
 * @param i, given residue to compute for.
 * @param pL, probabilities calculated so far.
 * @param lMatrix, conditional probability matrix.
 * @param rMatrix, conditional probability matrix.
 * @param lChild, leftChild id for pL.
 * @param rChild, rightChild id for pL.
 * @param msa, Multiple sequence alignment of our tree. Needed for deltaChar()
 * @param k, tree node used for calculating delta gap.
 * @param currSite, u column we are looking at.
 * @return likelihood as computed by formula 20.
 */
double probForNodeResidue(int i,double** pL,double** lMatrix, double** rMatrix,
        int lChild, int rChild, MSA* msa, TreeNode* k, int currSite);

/**
 * Computes the likelihood of node curentNode (k) given the likelihoods up to node k for
 * gap character. Formula (21) on the paper.
 * @param pL, probabilities calculated so far.
 * @param lMatrix, conditional probability matrix.
 * @param rMatrix, conditional probability matrix.
 * @param lChild, leftChild id for pL.
 * @param rChild, rightChild id for pL.
 * @param msa, Multiple sequence alignment of our tree. Needed for deltaChar()
 * @param currSite, u column we are looking at.
 *  * @param currSite, u column we are looking at.
 * @return likelihood as computed by formula 21.
 */
double probForNodeGap(double** pL, double** lMatrix, double** rMatrix, int lChild,
        int rChild, MSA* msa, TreeNode* k, int currSite);

/**
 * Same as probForNodeResidue except used for extra column containing all gaps therefore
 * we don not need to check deltaGap() to see if all children are gaps.
 * Computes the likelihood of node curentNode (k) given the likelihoods up to node k and
 * a current base. Formula (20) on the paper.
 * @param i, given residue to compute for.
 * @param pL, probabilities calculated so far.
 * @param lMatrix, conditional probability matrix.
 * @param rMatrix, conditional probability matrix.
 * @param lChild, leftChild id for pL.
 * @param rChild, rightChild id for pL.
 * @return likelihood as computed by formula 20.
 */
double probForNodeResidueAllGap(int i,double** pL,double** lMatrix, double** rMatrix,
        int lChild, int rChild);

/** Same as probForNodeGap except used for extra column containing all gaps therefore
 * we don not need to check deltaGap() to see if all children are gaps.
 * @param pL, probabilities calculated so far.
 * @param lMatrix, conditional probability matrix.
 * @param rMatrix, conditional probability matrix.
 * @param lChild, leftChild id for pL.
 * @param rChild, rightChild id for pL.
 * @return likelihood as computed by formula 21.
 */
double probForNodeGapAllGap(double** pL,double** lMatrix, double** rMatrix,
        int lChild, int rChild);

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
double singleEventCondProb(int j,int i, double branchLength,
        double* freqs,double* params);

/**
 * Calculates the function xi based on branch length, lambda and mu (formula 7)
 * @param branchLength, given branch length from j to i (t)
 * @param mu, deletion rate.
 * @param lambda, insertion rate.
 * @return value.
 */
double xi(double branchLength,double mu,double lambda);

/**
 * Calculates the function gamma based on branch length, lambda and mu (formula 6)
 * @param branchLength, given branch length from j to i (t)
 * @param mu, deletion rate.
 * @param lambda, insertion rate.
 * @return value.
 */
double gammaML(double branchLength,double mu,double lambda);
  
 /**
 * Calculates formula (9) on dnaML-erate paper.
 * param mod, model representing our phylogenetic tree.
 * @param i, assumed alphabet letter for probability.
 * @param j,given nucleotide.
 * @param t, branch length.
 * @param freqs, array frequencies.
 * @param params, array of parameters: [mu,lambda,alpha,beta]
 * @return computed formula (9) for M_t(i,j)
 */
double epsilonProbability(int j,int i,double t,double* freqs,double* params);

/**
 * Recursive function. Given a node on a tree will return either true if node was
 * a gap or false if node was not a gap. Used by probForNodeResidue and probForNodeGap.
 * @param k, node being looked at.
 * @param msa, we need the MSA to see what the characters look like at tree.
 * @return predicate stating whether all leaves from k down were gaps or not.
 */


int deltaGap(TreeNode* k,MSA* msa,int currSite);

/**
 * Returns whether the passed char was a gap ;)
 * @param c, character
 * @return 1 for true, 0 for false.
 */
int isGap(char c);

/**
 * Given a species name, a multiple sequence alinment and the index of the column where
 * you want the character from. Returns the character belonging to that species at that
 * position.
 * @param name, name of the specie.
 * @param msa, multiple sequence alignment.
 * @param index, site in the MSA to get char from.
 * @return char at that position.
 */
char getCharacterForSpecie(char* name,MSA* msa,int index);

/* Pu(L_k,i) = { 1, if leak k has residue i at position u;
 *              { 0, otherwise;
 * For gaps: Pu(L_k,'-') = 0  // Equation 18 on paper.
 */
int probabilityOfLeaf(char currentChar,int observedState,int iResidue);

 

#endif
