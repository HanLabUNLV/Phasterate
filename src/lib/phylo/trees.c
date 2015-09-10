/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "stacks.h"
#include "trees.h"
#include "misc.h"
#include "stringsplus.h"


static int idcounter = 0;
/* NOTE: when tree is parsed from Newick file, node ids are assigned
   sequentially in a preorder traversal.  Some useful properties
   result.  For example, if two nodes u and v are such that v->id >
   u->id, then the first ancestor a of v s.t. a->id < u->id is the LCA
   of u and v (see tr_lca) */

/* coords for postscript printing */
/* top-left x */
#define TL_X 15                 

/* top-left y */
#define TL_Y 10                 
                                
/* bottom-right x */
#define BR_X 500                
                                
/* bottom-right y */
#define BR_Y 700                
                                
/*Local Functions declarations. Should not be seen through the header file.*/
/**
 * Recurse through the tree setting all lists for common ancestors. This function
 * allocates memory that needs freeing! 
 * Like always we assume the tree is properly formed.
 * @param currentNode: node of tree to set @childrenSet for.
 * @return void: Sets the @childrenSet field for all nodes.
 */
void populateCommonAncestor(TreeNode* currentNode);

/**
 * Helper function for @findCommonAncestor. Does the work recursively. Finds and returns
 * the common ancestor for a tree that has ran @populateCommonAncestor.
 * @param currentNode: Current node we are recursing over. User should probably pass root.
 * @param nodeOneIndex: Child one to find common ancestor to with respect to child two.
 * @param nodeTwoIndex: Child two to find common ancestor to with respect to child two.
 * @return index of Common ancestor.
 */
int findCommonAncestorHelper(TreeNode* currentNode, int nodeOneIndex, int nodeTwoIndex);

/**
 * Given a node and a common ancestor of that node. It will iterate up adding values to
 * find the distance from our node to the ancestor.
 * @param node
 * @param commonAncestorId
 * @return 
 */
double addUpToCommonAncestor(TreeNode* node, int commonAncestorId);

/**
 * Given an allocated node. Populate by copying the important values of the original root.
 * @param root: Original Root of tree.
 * @param newRoot: new root to copy values to.
 * return: void.
 */
void initNewRoot(TreeNode* root, TreeNode* newRoot);

/**
 * Given the tree information it will iterate from the leaf until the common ancestor
 * attempting to insert the new root. If successful the new root is returned and the
 * lengths of the two children to the new root are permanently altered. We cannot know
 * ahead of time which nodes are changed (that's what this function computes!). Anyway
 * if you are here then that's what you wanted. If it fails to insert then we return NULL
 * and nothing is changed. Else we also populate newLeftDistance and newRightDistance with
 * the correct lengths.
 * @param root: The root of the entire tree.
 * @param leaf: Leaf to start attempting insertion from.
 * @param maxDistance: Distance from leaf to other leaf. This should be the total distance
 * as it will be halved.
 * @param commonAncestor: Common ancestor to attempt insertion up to.
 * @param newLeftDistance: will filled  with new length if successful in inserting root.
 * @param newRightDistance: will filled  with new length if successful in inserting root.
 * @return: new root if successful otherwise NULL.
 */
TreeNode* attemptInsert(TreeNode* root, TreeNode* leaf, double maxDistance,
        int commonAncestor, double* newLeftDistance, double* newRightDistance);

/**
 * Helper function to @midPointRooting() should not be used otherwise! This function
 * fixes the messed up length of some children and nodes as explained in
 * @midPointRooting(). It does so by recusively descending through the tree making sure
 * the lengths and parents are set properly.
 * @param node: Tree to recurse through fixing lengths.
 * @return: void, operates on passed tree.
 */
void fixParentsAndLengths(TreeNode* node);

/** Parse a tree from a file in Newick (New Hampshire) format */
TreeNode *tr_new_from_file(FILE *f) { 
  String *s = str_new(STR_VERY_LONG_LEN);
  TreeNode *retval;

  str_slurp(s, f);

  str_double_trim(s);
  if (s->chars[0] != '(')
    die("ERROR: This doesn't look like a tree (Newick format): \"%s\".\n", 
        s->chars);

  if (s->chars[s->length-1] == ';') 
    s->chars[--s->length] = '\0';

  retval = tr_new_from_string(s->chars);
  str_free(s);
  return retval;
}

/** Parse a Newick-formatted tree from a character string */
TreeNode *tr_new_from_string(const char *treestr) { 
  TreeNode *root, *node, *newnode;
  int i, in_distance = FALSE, in_label=FALSE, len = strlen(treestr), nopen_parens = 0,
    nclose_parens = 0, already_allowed = FALSE;
  char c;
  String *diststr = str_new(STR_SHORT_LEN), *labelstr = str_new(STR_SHORT_LEN);

  tr_reset_id();
  root = tr_new_node(); root->nnodes = 1;
  node = root;
  for (i = 0; i < len; i++) {
    c = treestr[i];

    if (in_distance) {
      if (c != '(' && c != ',' && c != ')' && c != ':' && c != '#' && c != '!') {
        str_append_char(diststr, c);
        continue;
      }
      else {
	str_double_trim(diststr);
        if (str_as_dbl(diststr, &node->dparent) != 0)
          die("ERROR: Can't parse distance in tree (\"%s\").\n", 
              diststr->chars);
	in_distance = FALSE;
      }
    }

    if (in_label) {
      if (c != '(' && c != ',' && c != ')' && c != ':' && c != '#' && c != '!') {
        str_append_char(labelstr, c);
        continue;
      }
      else {
	str_double_trim(labelstr);
	node->label = copy_charstr(labelstr->chars);
	in_label = FALSE;
      }
    }

    if (c == '(') {
      tr_add_child(node, newnode = tr_new_node());
      node = newnode;
      root->nnodes++;
      nopen_parens++;
    }
    else if (c == ',') {
      if (node->parent == NULL)
        die("ERROR: invalid binary tree (check parens).\n");
      if (node->parent->lchild != NULL && node->parent->rchild != NULL){
        if (node->parent == root && !already_allowed)
          already_allowed = TRUE;
        else
          die("ERROR (tree parser): invalid binary tree (too many children)\n");
                                /* we'll prohibit multinary
                                   branchings, except that we'll allow
                                   a single trinary branch immediately
                                   below the root (common with
                                   reversible models) */
      }

      tr_add_child(node->parent, newnode = tr_new_node());
      node = newnode;
      root->nnodes++;
    }
    else if (c == ')') {
      node = node->parent;
      nclose_parens++;
    }
    else if (c == ':') {
      str_clear(diststr);
      in_distance = TRUE;
    }
    else if (c == '#') {
      str_clear(labelstr);
      in_label = TRUE;
    }
    else if (c == '!') {
      node->hold_constant = 1;
    }
    else {			/* has to be part of name */
      if (!isspace(c) || node->name[0] != '\0')	/* avoid leading spaces */
        strncat(node->name, &c, 1);
    }
  }

  if (nopen_parens != nclose_parens)
    die("ERROR: mismatching parens in tree.\n");

  tr_set_nnodes(root);
  str_free(diststr);
  str_free(labelstr);
  return root;
}

/* traverse tree to set nnodes at each node (postorder); also set
   height at each node, and create "nodes" list */
void tr_set_nnodes(TreeNode *tree) {
  Stack *stack;
  TreeNode *node;
  int i, j;

  tree->nodes = lst_new_ptr(tree->nnodes);
  for (i = 0; i < tree->nnodes; i++) lst_push_ptr(tree->nodes, NULL);
  stack = stk_new_ptr(tree->nnodes);
  stk_push_ptr(stack, tree);
  while ((node = stk_pop_ptr(stack)) != NULL) {
    if (! ((node->lchild == NULL && node->rchild == NULL) ||
           (node->lchild != NULL && node->rchild != NULL)))
      die("Invalid tree");

    if (node->lchild == NULL) {
      node->nnodes = 1;
      node->height = 0;
      if (node->id >= tree->nnodes) 
        for (j = tree->nnodes; j <= node->id; j++) 
          lst_push_ptr(tree->nodes, NULL); /* this hack necessary
                                              because original estimate
                                              of size of list may be
                                              wrong */
      lst_set_ptr(tree->nodes, node->id, node);
    }
    else if (node->lchild->nnodes != -1 && node->rchild->nnodes != -1) {
      node->nnodes = node->lchild->nnodes + node->rchild->nnodes + 1;
      node->height = max(node->lchild->height, node->rchild->height) + 1;

      for (j = tree->nnodes; j <= node->id; j++) 
        lst_push_ptr(tree->nodes, NULL); /* this hack necessary
                                            because original estimate
                                            of size of list may be
                                            wrong */
      lst_set_ptr(tree->nodes, node->id, node);
    }
    else {			/* internal node whose children have
                       not yet been visited */
      stk_push_ptr(stack, node);
      stk_push_ptr(stack, node->lchild);
      stk_push_ptr(stack, node->rchild);
    }
  }
  stk_free(stack);
}

/* Create and initialize a new tree node */
TreeNode *tr_new_node() {
  TreeNode *n = (TreeNode*)smalloc(sizeof(TreeNode));
  n->parent = n->lchild = n->rchild = NULL;
  n->data = NULL;
  n->id = idcounter++;
  n->name[0] = '\0';
  n->dparent = 0;
  n->nnodes = -1;
  n->height = 0;
  n->label = NULL;
  n->nodes = n->preorder = n->inorder = n->postorder = NULL;
  n->hold_constant = 0;
  n->childList = NULL;
  return(n);
}

/* Add specified child to specified parent, creating all requisite
   links.  If the parent already has two children, add a new node
   (to simulate an nary tree with a binary tree) */
void tr_add_child(TreeNode *parent, TreeNode *child) {
  if (parent->lchild == NULL) {
    parent->lchild = child;
  }
  else if (parent->rchild == NULL) {
    parent->rchild = child;
  }
  else {
    /* add intermediate node to accommodate extra child */
    TreeNode *tmp = tr_new_node();
    tmp->lchild = parent->lchild;
    tmp->rchild = parent->rchild;
    tmp->parent = parent;
    tmp->lchild->parent = tmp;
    tmp->rchild->parent = tmp;
    parent->lchild = tmp;
    parent->rchild = child;
  }
  child->parent = parent;
}


/*  returns a newly allocated char * with newick representation of
    tree.
 */
char *tr_to_string(TreeNode *root, int show_branch_lengths) {
  char *str;
  str = smalloc((STR_MED_LEN+20)*(root->nnodes+1)*sizeof(char));
  str[0] = '\0';
  tr_to_string_recur(str, root, show_branch_lengths);
  str = srealloc(str, (strlen(str)+1)*sizeof(char));
  return str;
}


/* recursive subroutine used by tr_to_string */
void tr_to_string_recur(char *str, TreeNode *n, int show_branch_lengths) {
  char temp[100];
  if (!((n->lchild == NULL && n->rchild == NULL) ||
	(n->lchild != NULL && n->rchild != NULL)))
    die("ERROR: tr_to_string_recur, either both children should be NULL or neither\n");
  if (n->lchild != NULL) {
    strcat(str, "(");
    tr_to_string_recur(str, n->lchild, show_branch_lengths);
    strcat(str, ",");
    tr_to_string_recur(str, n->rchild, show_branch_lengths);
    strcat(str, ")");
    strcat(str, n->name);
  }
  else {
    strcat(str, n->name);
  }
  if (show_branch_lengths && n->parent != NULL) {
    sprintf(temp, ":%g", n->dparent);
    strcat(str, temp);
  }
  if (n->label != NULL) {
    sprintf(temp, " # %s", n->label);
    strcat(str, temp);
  }
}


/** Print tree in New Hampshire format. */
void tr_print(FILE* f, TreeNode *root, int show_branch_lengths) {
  int len;

  /* It's simplest to do this recursively. */
  tr_print_recur(f, root, show_branch_lengths);
  len = strlen(root->name);
  if (len == 0 || root->name[len-1] != ';')
    fprintf(f, ";");
  fprintf(f, "\n");
}

/* Recursive subroutine used by print_tree */
void tr_print_recur(FILE* f, TreeNode *n, int show_branch_lengths) {

  if (!((n->lchild == NULL && n->rchild == NULL) || 
	(n->lchild != NULL && n->rchild != NULL)))
    die("ERROR tr_print_recur: either both children should be NULL or neither\n");

  if (n->lchild != NULL) {
    fprintf(f, "(");
    tr_print_recur(f, n->lchild, show_branch_lengths);
    fprintf(f, ",");
    tr_print_recur(f, n->rchild, show_branch_lengths);
    fprintf(f, ")");
    fprintf(f, "%s", n->name);
  }
  else {
    fprintf(f, "%s", n->name);
  }

  if (show_branch_lengths && n->parent != NULL)
    fprintf(f, ":%g", n->dparent);
  if (n->hold_constant)
    fprintf(f, "!");
  if (n->label != NULL)
    fprintf(f, " # %s", n->label);
}

/** Free memory for tree */
void tr_free(TreeNode *tree) {
  Stack *stack;
  TreeNode *n;
  stack = stk_new_ptr(tree->nnodes);
  stk_push_ptr(stack, tree);
  while ((n = stk_pop_ptr(stack)) != NULL) {
    if (n->lchild != NULL) stk_push_ptr(stack, n->lchild);
    if (n->rchild != NULL) stk_push_ptr(stack, n->rchild);
    if (n->nodes != NULL) lst_free(n->nodes);
    if (n->preorder != NULL) lst_free(n->preorder);
    if (n->inorder != NULL) lst_free(n->inorder);
    if (n->postorder != NULL) lst_free(n->postorder);
    if (n->label != NULL) sfree(n->label);
    sfree(n);
  }
  stk_free(stack);
}

void tr_reset_id() {
  idcounter = 0;
}

/** Copy tree */
void tr_cpy(TreeNode *dest, TreeNode *src) {
  Stack *stack, *nodes, *cpystack;
  TreeNode *n, *ncpy, *lcpy, *rcpy;

  /* first just flatten dest into a list; we won't try to match
   * corresponding nodes, point is just to reuse memory */
  stack = stk_new_ptr(src->nnodes);
  nodes = stk_new_ptr(src->nnodes);
  cpystack = stk_new_ptr(src->nnodes);
  stk_push_ptr(stack, dest); 
  while ((n = stk_pop_ptr(stack)) != NULL) {
    if (n != dest)              /* leave the root aside */
      stk_push_ptr(nodes, n);

    if (n->lchild != NULL)
      stk_push_ptr(stack, n->lchild);
    if (n->rchild != NULL)
      stk_push_ptr(stack, n->rchild);
  }

  /* now copy node by node */
  tr_node_cpy(dest, src);       /* copy root */
  stk_push_ptr(stack, src); 
  stk_push_ptr(cpystack, dest);
  while ((n = stk_pop_ptr(stack)) != NULL) {
    ncpy = stk_pop_ptr(cpystack); 

    if (n->lchild != NULL) {
      lcpy = stk_pop_ptr(nodes);
      tr_node_cpy(lcpy, n->lchild);
      tr_add_child(ncpy, lcpy);
      stk_push_ptr(stack, n->lchild);
      stk_push_ptr(cpystack, lcpy);
    }
    if (n->rchild != NULL) {
      rcpy = stk_pop_ptr(nodes);
      tr_node_cpy(rcpy, n->rchild);
      tr_add_child(ncpy, rcpy);
      stk_push_ptr(stack, n->rchild);
      stk_push_ptr(cpystack, rcpy);
    }
  }

  stk_free(stack);
  stk_free(cpystack);
  stk_free(nodes);
}

/** Create a new tree that's a copy of another one */
TreeNode *tr_create_copy(TreeNode *src) {
  Stack *stack, *cpystack;
  TreeNode *n, *ncpy, *lcpy, *rcpy, *dest;
  
  tr_reset_id();

  stack = stk_new_ptr(src->nnodes);
  cpystack = stk_new_ptr(src->nnodes);
  stk_push_ptr(stack, src); 
  dest = tr_new_node();
  tr_node_cpy(dest, src);
  stk_push_ptr(cpystack, dest);
  while ((n = stk_pop_ptr(stack)) != NULL) {
    ncpy = stk_pop_ptr(cpystack);
    if (n->lchild != NULL) {
      lcpy = tr_new_node();
      tr_node_cpy(lcpy, n->lchild);
      tr_add_child(ncpy, lcpy);
      stk_push_ptr(stack, n->lchild);
      stk_push_ptr(cpystack, lcpy);
    }
    if (n->rchild != NULL) {
      rcpy = tr_new_node();
      tr_node_cpy(rcpy, n->rchild);
      tr_add_child(ncpy, rcpy);
      stk_push_ptr(stack, n->rchild);
      stk_push_ptr(cpystack, rcpy);
    }
  }  
  stk_free(stack);
  stk_free(cpystack);

  /* leave inorder, preorder, postorder NULL at dest; will be created
     on demand (can't copy from src because pointers will be wrong) */

  dest->nnodes = src->nnodes;
  tr_set_nnodes(dest);
  return dest;
}


/* Copy contents of tree node (ignore pointers) */
void tr_node_cpy(TreeNode *dest, TreeNode *src) {
  dest->parent = dest->lchild = dest->rchild = NULL;
  dest->id = src->id;
  strcpy(dest->name, src->name); 
  dest->dparent = src->dparent;
  if (src->label != NULL) dest->label = copy_charstr(src->label);
  dest->hold_constant = src->hold_constant;
  /* don't copy data, nnodes, height, preorder, inorder, postorder */
}


/** Print tree in Newick format.  This version imposes an
   ordering on the leaves (useful when comparing trees that have been
   rearranged).  At every internal node, we store the name of the leaf
   beneath it that comes first alphanumerically.  When recursively
   printing the tree, at each internal node, we call its children in
   the order of these names.  */
void tr_print_ordered(FILE* f, TreeNode *root, int show_branch_lengths) {
  int *left_right, *mark;
  char **names;
  TreeNode *n;
  Stack *stack;
  int i;

  left_right = (int*)smalloc(root->nnodes * sizeof(int));
  mark = (int*)smalloc(root->nnodes * sizeof(int));
  names = (char**)smalloc(root->nnodes * sizeof(char*));
  for (i = 0; i < root->nnodes; i++) { 
    left_right[i] = 0; 
    mark[i] = 0; 
    names[i] = NULL;
  }

  stack = stk_new_ptr(root->nnodes);
  stk_push_ptr(stack, root);
  while ((n = stk_pop_ptr(stack)) != NULL) {

    if (!((n->lchild == NULL && n->rchild == NULL) ||
	  (n->lchild != NULL && n->rchild != NULL)))
      die("ERROR tr_print_ordered: either both children should be NULL or neither\n");
	
    if (n->lchild == NULL) {
      names[n->id] = n->name;
      mark[n->id] = 1;
    }
    else if (mark[n->lchild->id] == 1 && mark[n->rchild->id] == 1) {
      if (names[n->rchild->id] == NULL)
        left_right[n->id] = 1;
      else if (names[n->lchild->id] == NULL) 
        left_right[n->id] = 0;
      /* now we know neither is NULL */
      else if (strcmp(names[n->lchild->id], names[n->rchild->id]) <= 0)
        left_right[n->id] = 1;
      else
        left_right[n->id] = 0;

      names[n->id] = left_right[n->id] == 1 ? names[n->lchild->id] :
        names[n->rchild->id];
      
      mark[n->id] = 1;
    }
    else {                      /* push back on stack for post-order
                                 * traversal */
      stk_push_ptr(stack, n);
      stk_push_ptr(stack, n->rchild);
      stk_push_ptr(stack, n->lchild);
    }
  }

  tr_print_ordered_recur(f, root, left_right, show_branch_lengths);
  fprintf(f, ";\n");
  
  stk_free(stack);
  sfree(left_right);
  sfree(mark);
  sfree(names);
}

/* Recursive subroutine for tr_print_ordered */
void tr_print_ordered_recur(FILE* f, TreeNode *n, int *left_right,
                            int show_branch_lengths) {

  if (!((n->lchild == NULL && n->rchild == NULL) || 
	(n->lchild != NULL && n->rchild != NULL)))
    die("ERROR tr_print_ordered_recur: either both children should be NULL or neither\n");

  if (n->lchild != NULL) {
    fprintf(f, "(");
    if (left_right[n->id]) {
      tr_print_ordered_recur(f, n->lchild, left_right, show_branch_lengths);
      fprintf(f, ",");
      tr_print_ordered_recur(f, n->rchild, left_right, show_branch_lengths);
    }
    else {
      tr_print_ordered_recur(f, n->rchild, left_right, show_branch_lengths);
      fprintf(f, ",");
      tr_print_ordered_recur(f, n->lchild, left_right, show_branch_lengths);
    }
    fprintf(f, ")");
  }
  else {
    fprintf(f, "%s", n->name);
  }

  if (show_branch_lengths)
    fprintf(f, ":%g", n->dparent);

  if (n->hold_constant)
    fprintf(f, "!");

  if (n->label != NULL)
    fprintf(f, " # %s", n->label);
}

/** Obtain a list representing a preorder traversal of the tree.
    Allows for simplified loops (no stacks!) and improved
    efficiency */
List *tr_preorder(TreeNode *tr) {
  if (tr->preorder == NULL) {   /* produce on demand */
    Stack *stack;
    TreeNode *n;

    tr->preorder = lst_new_ptr(tr->nnodes);
    stack = stk_new_ptr(tr->nnodes);
    stk_push_ptr(stack, tr);
    while ((n = stk_pop_ptr(stack)) != NULL) {
      if (!((n->lchild == NULL && n->rchild == NULL) ||
	    (n->lchild != NULL && n->rchild != NULL)))
	die("ERROR tr_preorder: either both children should be NULL or neither\n");
      if (n->lchild != NULL) {
        stk_push_ptr(stack, n->rchild);
        stk_push_ptr(stack, n->lchild);
      }
      lst_push_ptr(tr->preorder, n);
    }
    stk_free(stack);
  }

  return tr->preorder;
}

/** Obtain a list representing an in-order traversal of the tree.
    Allows for simplified loops (no stacks!) and improved
    efficiency */
List *tr_inorder(TreeNode *tr) {
  if (tr->inorder == NULL) {    /* produce on demand */
    int i;
    int *mark;
    Stack *stack;
    TreeNode *n;

    tr->inorder = lst_new_ptr(tr->nnodes);
    mark = (int*)smalloc(tr->nnodes * sizeof(int));
    for (i = 0; i < tr->nnodes; i++) mark[i] = 0;
    stack = stk_new_ptr(tr->nnodes);
    stk_push_ptr(stack, tr);
    while ((n = stk_pop_ptr(stack)) != NULL) {
      if (!((n->lchild == NULL && n->rchild == NULL) ||
	    (n->lchild != NULL && n->rchild != NULL)))
	die("ERROR: tr_inorder: either both children should be NULL or neither\n");
      if (n->lchild != NULL && mark[n->lchild->id] == 0) {
        stk_push_ptr(stack, n->rchild);
        stk_push_ptr(stack, n);
        stk_push_ptr(stack, n->lchild);
      }
      else {
        lst_push_ptr(tr->inorder, n);
        mark[n->id] = 1;
      }
    }
    stk_free(stack);
    sfree(mark);
  }

  return tr->inorder;
}

/** Obtain a list representing a postorder traversal of the tree.
    Allows for simplified loops (no stacks!) and improved
    efficiency */
List *tr_postorder(TreeNode *tr) {
  if (tr->postorder == NULL) {  /* produce on demand */
    int i;
    int *mark;
    Stack *stack;
    TreeNode *n;

    tr->postorder = lst_new_ptr(tr->nnodes);
    mark = (int*)smalloc(tr->nnodes * sizeof(int));
    for (i = 0; i < tr->nnodes; i++) mark[i] = 0;
    stack = stk_new_ptr(tr->nnodes);
    stk_push_ptr(stack, tr);
    while ((n = stk_pop_ptr(stack)) != NULL) {
      if (!((n->lchild == NULL && n->rchild == NULL) ||
	    (n->lchild != NULL && n->rchild != NULL)))
	die("ERROR tr_postorder: either both children should be NULL or neither\n");
      if (n->lchild != NULL && mark[n->lchild->id] == 0) {
        stk_push_ptr(stack, n); /* order? */
        stk_push_ptr(stack, n->rchild);
        stk_push_ptr(stack, n->lchild);
      }
      else {
        lst_push_ptr(tr->postorder, n);
        mark[n->id] = 1;
      }
    }
    stk_free(stack);
    sfree(mark);
  }
  
  return tr->postorder;
}

/** Provide x-y coordinates for layout. */
void tr_layout_xy(TreeNode *tree, 
                  int x0,       /* Upper left x bound */
                  int y0,       /* Upper left y bound  */
                  int x1,       /* Lower right x bound */
                  int y1,       /* Lower right y bound  */
                  int *x,       /* On return, will contain
                                   x-coordinates for nodes, in order
                                   of tree->nodes.  Must be
                                   preallocated. */
                  int *y,       /* On return, will contain
                                   y-coordinates for nodes, in order
                                   of tree->nodes.  Must be
                                   preallocated. */
                  int use_branch_lens, 
                                /* If TRUE, tree will be laid out
                                   such that edges are proportional to
                                   branch lengths (dparent
                                   attributes) */
                  int horizontal
                                /* If TRUE, tree will be laid out
                                   with root on left and leaves on
                                   right; otherwise, root will be at
                                   top and leaves at bottom */
                  ) {

  int delt_x, delt_y, i; 
  List *traversal; 
  TreeNode *n;
  double scale = 0;

  delt_x = x1 - x0; 
  delt_y = y0 - y1; 

  /* if scaling according to branch lens, we need the total height of
     the tree (in terms of branch lengths, not levels) */
  if (use_branch_lens) {
    double *total_height = (double*)smalloc(tree->nnodes * sizeof(double));
    traversal = tr_postorder(tree);
    for (i = 0; i < tree->nnodes; i++) {
      n = lst_get_ptr(traversal, i);
      if (n->lchild == NULL)
        total_height[n->id] = 0;
      else
        total_height[n->id] = 
          max(total_height[n->lchild->id] + n->lchild->dparent, 
              total_height[n->rchild->id] + n->rchild->dparent);
    }
    scale = (horizontal == 1 ? abs(delt_x)/total_height[tree->id] :
             abs(delt_y)/total_height[tree->id]);
    sfree(total_height);
  }

  /* set x coords (or y's if horizontal) by spacing evenly in an
     inorder traversal; if not scaling according to branch lens, can
     simultaneously set opposite coordinates, using the "height" of
     each node in the tree */
  traversal = tr_inorder(tree); 
  for (i = 0; i < tree->nnodes; i++) {  
    n = lst_get_ptr(traversal, i); 
    if (horizontal == 0) {
      x[n->id] = x0 + i*delt_x/tree->nnodes;     
      if (!use_branch_lens) y[n->id] = y1 + n->height * delt_y/tree->height;
    }
    else {
      y[n->id] = y1 + i*delt_y/tree->nnodes;     
      if (!use_branch_lens) x[n->id] = x1 - n->height * delt_x/tree->height;
    }
  } 

  /* if scaling according to branch lens, set y coords (or x's if
     horizontal) incrementally in a preorder traversal */
  if (use_branch_lens) {
    traversal = tr_preorder(tree);
    for (i = 0; i < tree->nnodes; i++) {  
      n = lst_get_ptr(traversal, i); 
      if (horizontal == 0) 
        y[n->id] = n->parent == NULL ? y0 :
          y[n->parent->id] - n->dparent*scale;
      else
        x[n->id] = n->parent == NULL ? x0 :
          x[n->parent->id] + n->dparent*scale;
    }
  }
} 

/** Print a (very basic!) postscript rendering of a tree. */
void tr_print_ps(FILE *F,       /* Destination file */
                 TreeNode *tree, 
                                /**< Tree root */
                 int show_branch_lens,
                                /* Whether to print branch lengths
                                   by edges */
                 int square_branches,
                                /* If TRUE, branches will be
                                   right-angled, otherwise will be
                                   diagonal */
                 int use_branch_lens, 
                                /* If TRUE, tree will be laid out
                                   such that edges are proportional to
                                   branch lengths (dparent
                                   attributes) */
                 int horizontal_layout
                                /* If TRUE, tree will be laid out
                                   with root on left and leaves on
                                   right; otherwise, root will be at
                                   top and leaves at bottom */
                 ) {
  int i, xoffset, yoffset;
  int *x, *y;
  List *traversal;

  x = (int*)smalloc(tree->nnodes * sizeof(int));
  y = (int*)smalloc(tree->nnodes * sizeof(int));
  tr_layout_xy(tree, TL_X, TL_Y, BR_X, BR_Y, x, y, 
               use_branch_lens, horizontal_layout);

  /* print header */
  fprintf(F, "%%!\n\
1 setlinecap 1 setlinejoin 1 setlinewidth 0 setgray\n\
/basefont /Times-Roman findfont 12 scalefont def\n\
50 50 translate\n\
basefont setfont\n");

  traversal = tr_postorder(tree);
  for (i = 0; i < lst_size(traversal); i++) {
    TreeNode *n = lst_get_ptr(traversal, i);
    if (n->lchild == NULL) {    /* leaf */

      /* offsets for labels; should parameterize this a bit better,
         perhaps accepting font size as a function parameter */
      if (horizontal_layout) {
        xoffset = 10;
        yoffset = -6;
      }
      else {
        xoffset = -3 * (n->name != NULL ? strlen(n->name) : 0);
        yoffset = -18;
      }

      /* draw leaf label */
      fprintf(F, "%d %d moveto\n(%s) show\n", x[n->id]+xoffset, 
              y[n->id]+yoffset, 
              n->name != NULL ? n->name : "");
    }
    else {                      /* internal node */
      /* draw branches from node to each of its children; this is all
         a bit messy but it will have to do for now ... */
      if (square_branches) {
        if (horizontal_layout) {
          fprintf(F, "%d %d moveto\n%d %d lineto stroke\n", x[n->lchild->id], 
                  y[n->lchild->id], x[n->id], y[n->lchild->id]);
          fprintf(F, "%d %d moveto\n%d %d lineto stroke\n", x[n->id], 
                  y[n->lchild->id], x[n->id], y[n->id]);
          fprintf(F, "%d %d moveto\n%d %d lineto stroke\n", x[n->rchild->id], 
                  y[n->rchild->id], x[n->id], y[n->rchild->id]);
          fprintf(F, "%d %d moveto\n%d %d lineto stroke\n", x[n->id], 
                  y[n->rchild->id], x[n->id], y[n->id]);

          if (show_branch_lens) {
            fprintf(F, "%d %d moveto\n(%.4f) show\n", 
                    (x[n->lchild->id]-x[n->id])/2 + x[n->id], 
                    y[n->lchild->id] + 6, n->lchild->dparent);
            fprintf(F, "%d %d moveto\n(%.4f) show\n", 
                    (x[n->rchild->id]-x[n->id])/2 + x[n->id], 
                    y[n->rchild->id] + 6, n->rchild->dparent);
          }
        }
        else {
          fprintf(F, "%d %d moveto\n%d %d lineto stroke\n", x[n->lchild->id], 
                  y[n->lchild->id], x[n->lchild->id], y[n->id]);
          fprintf(F, "%d %d moveto\n%d %d lineto stroke\n", x[n->lchild->id], 
                  y[n->id], x[n->id], y[n->id]);
          fprintf(F, "%d %d moveto\n%d %d lineto stroke\n", x[n->rchild->id], 
                  y[n->rchild->id], x[n->rchild->id], y[n->id]);
          fprintf(F, "%d %d moveto\n%d %d lineto stroke\n", x[n->rchild->id], 
                  y[n->id], x[n->id], y[n->id]);

          if (show_branch_lens) {
            fprintf(F, "%d %d moveto\n(%.4f) show\n", 
                    x[n->lchild->id] + 6,
                    (y[n->id]-y[n->lchild->id])/2 + y[n->lchild->id],
                    n->lchild->dparent);
            fprintf(F, "%d %d moveto\n(%.4f) show\n", 
                    x[n->rchild->id] + 6,
                    (y[n->id]-y[n->rchild->id])/2 + y[n->rchild->id],
                    n->rchild->dparent);
          }
        }
      }
      else {                    /* diag branches */
        fprintf(F, "%d %d moveto\n%d %d lineto stroke\n", x[n->lchild->id], 
                y[n->lchild->id], x[n->id], y[n->id]);
        fprintf(F, "%d %d moveto\n%d %d lineto stroke\n", x[n->rchild->id], 
                y[n->rchild->id], x[n->id], y[n->id]);

        if (show_branch_lens) {
          if (horizontal_layout) {
            fprintf(F, "%d %d moveto\n(%.4f) show\n", 
                    (x[n->lchild->id]-x[n->id])/2 + x[n->id], 
                    (y[n->id]-y[n->lchild->id])/2 + y[n->lchild->id] + 3, 
                    n->lchild->dparent);
            fprintf(F, "%d %d moveto\n(%.4f) show\n", 
                    (x[n->rchild->id]-x[n->id])/2 + x[n->id], 
                    (y[n->rchild->id]-y[n->id])/2 + y[n->id] - 12, 
                    n->rchild->dparent);
          }
          else {
            fprintf(F, "%d %d moveto\n(%.4f) show\n", 
                    (x[n->lchild->id]-x[n->id])/2 + x[n->id] + 6,
                    (y[n->id]-y[n->lchild->id])/2 + y[n->lchild->id],
                    n->lchild->dparent);
            fprintf(F, "%d %d moveto\n(%.4f) show\n", 
                    (x[n->id]-x[n->rchild->id])/2 + x[n->rchild->id] + 6,
                    (y[n->id]-y[n->rchild->id])/2 + y[n->rchild->id],
                    n->rchild->dparent);
          }
        }
      }
    }
  }
  fprintf(F, "showpage\n");     /* complete PS file */

  sfree(x);
  sfree(y);
}

/** Compute and return sum of lengths at all edges */
double tr_total_len(TreeNode *t) {
  double retval = 0;
  int i;
  for (i = 0; i < t->nnodes; i++) {
    TreeNode *n = lst_get_ptr(t->nodes, i);
    if (n->parent != NULL) 
      retval += n->dparent;
  }
  return retval;
}

/** Compute and return sum of lengths of edges in subtree below
    given node */
double tr_total_len_subtree(TreeNode *sub_root) {
  TreeNode *n;
  Stack *stack = stk_new_ptr(sub_root->nnodes);
  double retval = 0;
  stk_push_ptr(stack, sub_root->lchild);
  stk_push_ptr(stack, sub_root->rchild);
  while ((n = stk_pop_ptr(stack)) != NULL) {
    if (n->lchild != NULL) {
      stk_push_ptr(stack, n->rchild);
      stk_push_ptr(stack, n->lchild);
    }
    retval += n->dparent;
  }
  stk_free(stack);
  return retval;
}

/** Compute and return maximum branch length in tree */
double tr_max_branchlen(TreeNode *t) {
  double retval = 0;
  int i;
  for (i = 0; i < t->nnodes; i++) {
    TreeNode *n = lst_get_ptr(t->nodes, i);
    if (n->parent == NULL) continue;
    if (n->dparent > retval)
      retval = n->dparent;
  }
  return retval;  
}


double tr_distance_to_root(TreeNode *node) {
  if (node->parent == NULL) return 0.0;
  return node->dparent + tr_distance_to_root(node->parent);
}

/** Return node having specified name or NULL if none found.  */
TreeNode *tr_get_node(TreeNode *t, const char *name) {
  int i;
  for (i = 0; i < t->nnodes; i++) {
    TreeNode *n = lst_get_ptr(t->nodes, i);
    if (n->name[0] != '\0' && !strcmp(n->name, name))
      return n;
  }
  return NULL;
}

/** Return node having specified id. */
TreeNode *tr_get_node_id(TreeNode *t, int id){
  int i;
  for (i = 0; i < t->nnodes; i++) {
    TreeNode *n = lst_get_ptr(t->nodes, i);
    if(n->id == id)
      return n;
  }
  return NULL;
}

/** Scale all branch lengths by constant factor. */
void tr_scale(TreeNode *t, double scale_const) {
  int i;
  for (i = 0; i < t->nnodes; i++) {
    TreeNode *n = lst_get_ptr(t->nodes, i);
    if (n->parent != NULL) 
      n->dparent *= scale_const;

  }
}

/** Scale all branch lengths by constant factor in subtree beneath
    given node. */
void tr_scale_subtree(TreeNode *t, TreeNode *sub, double scale_const,
		      int include_leading) {
  int i;
  List *inside = lst_new_ptr(t->nnodes);
  tr_partition_nodes(t, sub, inside, NULL);
  for (i = 0; i < lst_size(inside); i++) {
    TreeNode *n = lst_get_ptr(inside, i);
    if (n != sub || include_leading) n->dparent *= scale_const;
  }
  lst_free(inside);
}

/** Prune away all leaves whose names are in (or not in) the specified
    list.  Nodes will be removed and branches combined (branch lengths
    added) to restore as a proper binary tree.  */
void tr_prune(TreeNode **t,     /* Tree to prune (may be altered
                                   because root can change) */
              List *names,      /* List of names.  On return, will
                                   contain list of names of leaves
                                   that were pruned away.  */
              int all_but,       /* if FALSE, prune leaves *in*
                                   'names'; if TRUE, prune leaves *not
                                   in* 'names'  */
              int *id_map        /* if not NULL, should be allocated to 
				    the number of nodes in original tree.
				    On return, will be filled in with the
				    new id for each node */
              ) {

  TreeNode *n;
  int i, new_nnodes = (*t)->nnodes;
  int *is_leaf;
  List *traversal, *pruned_leaves = lst_new_ptr((*t)->nnodes / 2);

  /* first identify original leaves; will need to distinguish them
     from leaves that are created by pruning */
  is_leaf = smalloc((*t)->nnodes * sizeof(int));
  for (i = 0; i < (*t)->nnodes; i++) {
    n = lst_get_ptr((*t)->nodes, i);
    is_leaf[i] = (n->lchild == NULL && n->rchild == NULL);
  }

  /* get rid of nodes, preorder, and inorder lists (do now because
     root may change) */
  if ((*t)->nodes != NULL) { lst_free((*t)->nodes); (*t)->nodes = NULL; }
  if ((*t)->preorder != NULL) { lst_free((*t)->preorder); (*t)->preorder = NULL; }
  if ((*t)->inorder != NULL) { lst_free((*t)->inorder); (*t)->inorder = NULL; }

  /* remove nodes and combine branches in postorder traversal */
  traversal = tr_postorder(*t);
  for (i = 0; i < lst_size(traversal); i++) {
    TreeNode *n = lst_get_ptr(traversal, i);
    if (n->lchild == NULL && n->rchild == NULL){ /* missing both children */
      String *s;
      int prune;
      
      if (!is_leaf[n->id])      /* if not originally a leaf, must be pruned */
        prune = TRUE;
      else {
        s = str_new_charstr(n->name);
        prune = str_in_list(s, names);
        if (all_but) prune = !prune;

        if (prune) lst_push_ptr(pruned_leaves, s);
        else str_free(s);
      }

      if (prune) {
        if (n->parent == NULL) 
          *t = NULL;            /* entire tree has been pruned away! */
        else {
          if (n == n->parent->lchild) n->parent->lchild = NULL;
          else n->parent->rchild = NULL;
        }
        sfree(n);
        new_nnodes--;
      }
    }
    else if (n->lchild == NULL) { /* missing left child only */
      if (n->parent == NULL) {
	if (n != *t)
	  die("ERROR tr_prune n is not root!\n");
        n->rchild->parent = NULL; /* redefine root */
        *t = n->rchild;
        (*t)->dparent = 0;
      }
      else {                    /* mid-level node; remove and combine
                                   branch lengths */
        n->rchild->parent = n->parent;
        n->rchild->dparent += n->dparent;
        if (n == n->parent->lchild) n->parent->lchild = n->rchild;
        else n->parent->rchild = n->rchild;
      }
      sfree(n);
      new_nnodes--;
    }
    else if (n->rchild == NULL) { /* missing right child only */
      if (n->parent == NULL) {
	if (n != *t)
	  die("ERROR tr_prune n must be root\n");
        n->lchild->parent = NULL; /* redefine root */
        *t = n->lchild;         
        (*t)->dparent = 0;
      }
      else {                    /* mid-level node; remove and combine
                                   branch lengths */
        n->lchild->parent = n->parent;
        n->lchild->dparent += n->dparent;
        if (n == n->parent->lchild) n->parent->lchild = n->lchild;
        else n->parent->rchild = n->lchild;
      }
      sfree(n);
      new_nnodes--;
    }
  }

  /* finally, free postorder list */
  lst_free(traversal); 
  if (*t != NULL && (*t)->postorder != NULL) (*t)->postorder = NULL;

  /* reset ids, nodes, nnodes, heights */
  if (id_map != NULL) {
    for (i=0; i < lst_size(traversal); i++)
      id_map[i] = -1;
  }
  if (*t != NULL) {
    (*t)->nnodes = new_nnodes;
    traversal = tr_preorder(*t);
    for (i = 0; i < lst_size(traversal); i++) {
      n = lst_get_ptr(traversal, i);
      if (id_map != NULL)
	id_map[n->id] = i;
      n->id = i;
      if (n != *t) n->nnodes = -1;
    }
    tr_set_nnodes(*t);
  }

  lst_free_strings(names);
  lst_clear(names);
  for (i = 0; i < lst_size(pruned_leaves); i++) 
    lst_push_ptr(names, lst_get_ptr(pruned_leaves, i));

  lst_free(pruned_leaves);
  sfree(is_leaf);
}

void tr_prune_supertree(TreeNode **t, TreeNode *node) {
  List *prune_names=lst_new_ptr((*t)->nnodes);
  int i, *inSub = tr_in_subtree(*t, node);       
  String *s;
  TreeNode *n;

  for (i=0; i<(*t)->nnodes; i++) {
    n = lst_get_ptr((*t)->nodes, i);
    if (inSub[i]==0 && n->lchild == NULL) {
      s = str_new_charstr(n->name);
      lst_push_ptr(prune_names, s);
    }
  }
  tr_prune(t, prune_names, FALSE, NULL);
  lst_free_strings(prune_names);
  lst_free(prune_names);
  sfree(inSub);
}

void tr_prune_subtree(TreeNode **t, TreeNode *node) {
  List *prune_names=lst_new_ptr((*t)->nnodes);
  int i, *inSub = tr_in_subtree(*t, node);       
  String *s;
  TreeNode *n;

  for (i=0; i<(*t)->nnodes; i++) {
    n = lst_get_ptr((*t)->nodes, i);
    if (inSub[i]==0 && n->lchild == NULL) {
      s = str_new_charstr(n->name);
      lst_push_ptr(prune_names, s);
    }
  }
  tr_prune(t, prune_names, TRUE, NULL);
  lst_free_strings(prune_names);
  lst_free(prune_names);
}


/** Return the LCA of the given species.  Assumes ids are numbered in
    preorder (a node's parent always has a smaller id than it does and
    left descendants have smaller ids than right descendants). */
TreeNode *tr_lca(TreeNode *tree, List *names) {
  int i, min = tree->nnodes, max = -1, idx;
  String *tmpstr = str_new(STR_MED_LEN);
  TreeNode *n;
  int *found = smalloc(lst_size(names) * sizeof(int));

  for (i = 0; i < lst_size(names); i++) found[i] = FALSE;

  for (i = 0; i < tree->nnodes; i++) {
    n = lst_get_ptr(tree->nodes, i);
    if (n->lchild == NULL && n->rchild == NULL && n->name[0] != '\0') {
      str_cpy_charstr(tmpstr, n->name);
      if (str_in_list_idx(tmpstr, names, &idx)) {
        found[idx] = TRUE;
        if (n->id < min) min = n->id;
        if (n->id > max) max = n->id;
      }
    }
  }

  for (i = 0; i < lst_size(names); i++)
    if (!found[i])
      die("ERROR: species name not found in tr_lca ('%s')\n",
          ((String*)lst_get_ptr(names, i))->chars);

  /* now the LCA must be the first ancestor of the node with the max
     id that has an id smaller than the min id */
  for (n = lst_get_ptr(tree->nodes, max); n->id > min; n = n->parent);

  str_free(tmpstr);
  sfree(found);
  return n;
}

/** Given two trees, one of which is a (proper) subtree of the
    other, create a hybrid tree composed of the smaller tree and a
    scaled version of the larger tree.  First, a copy of the larger
    tree will be created and scaled such that the total branch length
    in the subtree in question is equal to the total branch length of
    the smaller tree.  Then, (a copy of) the smaller tree will be used
    in place of the overlapping subtree in the larger tree.  This
    function can be used to extrapolate from a small phylogeny for
    which accurate branch length estimation is possible (e.g., of
    eutherian mammals) to a larger phylogeny for which approximate
    branch length proportions are available, but absolute branch
    length estimates are not (e.g., of more distant vertebrates). */
TreeNode *tr_hybrid(TreeNode *sub, TreeNode *super) {
  TreeNode *retval, *lca, *sub_copy;
  double lfrac, sum;
  List *names = tr_leaf_names(sub);

  if (lst_size(names) <= 2)
    die("ERROR: subtree must have at least two leaves in tr_hybrid.\n");

  /* copy supertree then find LCA corresponding to subtree */
  retval = tr_create_copy(super);
  lca = tr_lca(retval, names);
  lst_free_strings(names);
  lst_free(names);

  sub_copy = tr_create_copy(sub);

  /* verify subtree is proper.  We know that all names in the subtree
     were found in the supertree (otherwise tr_lca will abort), so
     it's sufficient to verify that the numbers of nodes in the
     subtree and beneath the LCA are equal */
  if (lca->nnodes != sub->nnodes)
    die("ERROR: subtree must be proper in tr_hybrid (must contain all leaves beneath LCA in supertree).\n");

  if (lca == super) {           /* trivial case -- trees equal */
    tr_free(retval);
    return sub_copy;
  }

  /* scale supertree so that overlapping portions have equal total length */
  tr_scale(retval, tr_total_len(sub_copy) / tr_total_len_subtree(lca));

  /* now recombine */
  sub_copy->parent = lca->parent;
  sub_copy->dparent = lca->dparent;
  if (lca == lca->parent->lchild) lca->parent->lchild = sub_copy;
  else lca->parent->rchild = sub_copy;

  /* also ensure that the subtree is rooted with proportions equal to
     those of the original supertree (usually hard to get root right
     in subtree) */
  lfrac = lca->lchild->dparent / (lca->lchild->dparent + lca->rchild->dparent);
  sum = sub_copy->lchild->dparent + sub_copy->rchild->dparent;
  sub_copy->lchild->dparent = lfrac * sum;
  sub_copy->rchild->dparent = sum - sub_copy->lchild->dparent;

  tr_free(lca);                 /* this works recursively */

  return retval;
}

/* Scale a tree so that the total branch length of some subtree is as
   defined by a second tree.  Return value is scale factor.  The leaf
   names of the second tree ('sub') must be a subset of those in the
   first.  This function is similar to tr_extrapolate, but the branch
   length proportions of the larger tree are used without change and
   the smaller tree need not be a proper subtree of the larger
   tree. */
double tr_scale_by_subtree(TreeNode *tree, TreeNode *sub) {
  TreeNode *n;
  int i;
  double scale;
  List *leaf_names = tr_leaf_names(tree);
  List *sub_names = tr_leaf_names(sub);

  if (lst_size(sub_names) <= 2)
    die("ERROR: (tr_scale_for_subtree) subtree must have at least two leaves.\n");

  for (i = 0; i < lst_size(sub_names); i++)
    if (!str_in_list(lst_get_ptr(sub_names, i), leaf_names))
      die("ERROR: (tr_scale_for_subtree) leaf names in subtree must be subset of those in main tree.\n");

  /* Now scale */
  n = tr_create_copy(tree);
  tr_prune(&n, sub_names, TRUE, NULL);
  scale = tr_total_len(sub) / tr_total_len(n);
  tr_scale(tree, scale);
  tr_free(n);
  lst_free_strings(leaf_names); lst_free(leaf_names);
  lst_free_strings(sub_names); lst_free(sub_names);
  return scale;
}

/** Partition leaves of tree at (branch above) given node.  All
    descendant leaves of 'sub' will be added to 'inside' list and all
    non-descendants will be added to 'outside' list.  Lists must be
    pre-allocated. */
void tr_partition_leaves(TreeNode *tree, TreeNode *sub, List *inside, 
                         List *outside) {
  int i;
  TreeNode *n;
  int *mark = smalloc(tree->nnodes * sizeof(int));
  Stack *stack = stk_new_ptr(sub->nnodes);

  for (i = 0; i < tree->nnodes; i++) mark[i] = FALSE;

  lst_clear(inside);
  lst_clear(outside);
  stk_push_ptr(stack, sub);
  while ((n = stk_pop_ptr(stack)) != NULL) {
    if (n->lchild == NULL) {
      lst_push_ptr(inside, n);
      mark[n->id] = TRUE;
    }
    else {
      stk_push_ptr(stack, n->rchild);
      stk_push_ptr(stack, n->lchild);
    }
  }
  for (i = 0; i < tree->nnodes; i++) {
    n = lst_get_ptr(tree->nodes, i);
    if (n->lchild == NULL && !mark[n->id])
      lst_push_ptr(outside, n);
  }
  stk_free(stack);
  sfree(mark);
}

/** Similar to above, but partition all nodes; if either 'inside' or
    'outside' is NULL, it will be ignored */
void tr_partition_nodes(TreeNode *tree, TreeNode *sub, List *inside, 
                        List *outside) {
  int i;
  TreeNode *n;
  int *mark = smalloc(tree->nnodes * sizeof(int));
  Stack *stack = stk_new_ptr(sub->nnodes);

  for (i = 0; i < tree->nnodes; i++) mark[i] = FALSE;

  if (inside != NULL) lst_clear(inside);
  if (outside != NULL) lst_clear(outside);
  stk_push_ptr(stack, sub);
  while ((n = stk_pop_ptr(stack)) != NULL) {
    if (inside != NULL) lst_push_ptr(inside, n);
    mark[n->id] = TRUE;
    if (n->lchild != NULL) {
      stk_push_ptr(stack, n->rchild);
      stk_push_ptr(stack, n->lchild);
    }
  }
  if (outside != NULL) {
    for (i = 0; i < tree->nnodes; i++) {
      n = lst_get_ptr(tree->nodes, i);
      if (!mark[n->id])
        lst_push_ptr(outside, n);
    }
  }
  stk_free(stack);
  sfree(mark);
}


/* recursive version of tr_leaf_names that works for nodes that are
   not the root of the tree */
void tr_leaf_names_rec(TreeNode *tree, List *rv) {
  if (tree->lchild == NULL) 
    lst_push_ptr(rv, str_new_charstr(tree->name));
  else {
    tr_leaf_names_rec(tree->lchild, rv);
    tr_leaf_names_rec(tree->rchild, rv);
  }
}

/** Return a list of the leaf names in a given tree */
List *tr_leaf_names(TreeNode *tree) {
  List *retval = lst_new_ptr((tree->nnodes + 1) / 2);
  int i;
  if (tree->nodes == NULL) {
    // do this recursively if tree->nodes is NULL (usually true for non-root nodes)
    tr_leaf_names_rec(tree, retval);
  } else {
    for (i = 0; i < tree->nnodes; i++) {
      TreeNode *n = lst_get_ptr(tree->nodes, i);
      if (n->lchild == NULL && n->rchild == NULL)
	lst_push_ptr(retval, str_new_charstr(n->name));
    }
  }
  return retval;
}

/* Ensure all ancestral nodes have names.  If a node is unnamed, give
   it a name that is a concatenation of the name of a leaf from its
   left subtree and the name of a leaf from its right subtree.
   Leftmost decendants are selected, for lack of any better
   criterion.  */
void tr_name_ancestors(TreeNode *tree) {
  int i;
  List *traversal = tr_postorder(tree);
  char **repname = smalloc(tree->nnodes * sizeof(void*));
  for (i = 0; i < lst_size(traversal); i++) {
    TreeNode *n = lst_get_ptr(traversal, i);
    
    if ((n->lchild == NULL && n->rchild != NULL) ||
	(n->lchild != NULL && n->rchild == NULL))
      die("ERROR: malformed tree.\n");

    if (n->lchild == NULL) {
      if (n->name[0] == '\0') die("ERROR: unnamed leaf.\n");
      repname[n->id] = n->name;
    }  
    else {
      if (n->name[0] == '\0') {
	strcat(n->name, repname[n->lchild->id]);
	strcat(n->name, "-");
	strcat(n->name, repname[n->rchild->id]);
      }
      repname[n->id] = repname[n->lchild->id];
    }
  }
  sfree(repname);
}

/** Print verbose description of each node */
void tr_print_nodes(FILE *F, TreeNode *tree) {
  int i;
  List *l = tr_preorder(tree);
  for (i = 0; i < lst_size(l); i++) {
    TreeNode *n = lst_get_ptr(l, i);
    fprintf(F, "Node %d:\n", n->id);
    fprintf(F, "\tparent = %d\n", n->parent == NULL ? -1 : n->parent->id);
    fprintf(F, "\tlchild = %d\n", n->lchild == NULL ? -1 : n->lchild->id);
    fprintf(F, "\trchild = %d\n", n->rchild == NULL ? -1 : n->rchild->id);
    fprintf(F, "\tname = '%s'\n", n->name);
    fprintf(F, "\tdparent = %g", n->dparent);
    if (n->hold_constant) fprintf(F, " (constant)\n");
    else fprintf(F, "\n");
    if (n->label != NULL)
      fprintf(F, "\tlabel = %s\n", n->label);
    fprintf(F, "\n");
  }
}

/** Reroot tree.  Subtree originally beneath selected node will become
    right subtree of root, and remainder of tree will be left
    subtree. If include_branch == FALSE, the selected node will become
    the new root, and a zero-length branch to its right will connect
    it to its original subtree.  If instead include_branch == TRUE,
    then the branch above the selected node will also be included in the
    right subtree.  In this case, the selected node will become the
    right child of the new root and the branch in question will become
    the right branch beneath the new root.  The left branch beneath
    the new root will have length zero and will connect to the former
    parent of the selected node.  Warning: ids will not be altered, so
    they will no longer be consistent with a preorder traversal of the
    tree  */
void tr_reroot(TreeNode *tree, TreeNode *selected_node, int include_branch) {
  TreeNode *n, *p, *gp, *root_lchild, *root_rchild, *newroot;
  double d, d2;
  double root_d;
  int i;
  
  if (tree->nnodes < 3)
    die("Error: tr_reroot tree should have at least 3 nodes, has %i\n",
	tree->nnodes);

  if (tree == selected_node) {
    if (include_branch) 
      die("ERROR: strange call to tr_reroot -- rerooting at existing root with include_branch == TRUE.\n");
    else {
      tree->lchild->dparent += tree->rchild->dparent;
      tree->rchild->dparent = 0; 
      /* (make effect on branch lengths same as in true rerooting) */
      return;
    }
  }

  if (!include_branch &&
      (selected_node->lchild == NULL || selected_node->rchild == NULL))
    die("ERROR: cannot reroot at leaf unless including branch above it.\n");

  /* save left and right children of old root */
  root_lchild = tree->lchild;
  root_rchild = tree->rchild;
  root_d = root_lchild->dparent + root_rchild->dparent;

  /* move old root to position beneath selected node, connected by 0-length
     branch */
  tree->lchild = selected_node->lchild;
  tree->rchild = selected_node->rchild;
  tree->parent = selected_node;
  selected_node->rchild = tree;
  selected_node->lchild = NULL;
  tree->name[0] = '\0';		/* name may no longer make sense */
  /* branch lengths are okay */

  /* swap pointers above new root */
  n = selected_node;
  p = n->parent;
  d = n->dparent;               /* d is length of branch between n and p */
  while (p != tree) {
    /* save grandparent */
    gp = p->parent;
    d2 = p->dparent;            /* d2 is length of branch between p and gp */

    /* swap edge between n and parent */
    p->parent = n;    
    p->dparent = d;
    if (n->lchild == NULL) 
      n->lchild = p;
    else
      n->rchild = p;

    /* set extra pointer to NULL */
    if (p->lchild == n)
      p->lchild = NULL;
    else
      p->rchild = NULL;

    /* repeat */
    n = p;
    p = gp;    
    d = d2;
  }

  /* now n points to either root_rchild or root_lchild */
  if (!(n == root_rchild || n == root_lchild))
    die("ERROR tr_reroot n should be root->rchild or root->lchild\n");
  if (n == root_lchild) {
    if (n->lchild == NULL) 
      n->lchild = root_rchild;
    else
      n->rchild = root_rchild;

    root_rchild->parent = n;
    root_rchild->dparent = root_d;
  }
  else {                      /* n == root_rchild */
    if (n->lchild == NULL) 
      n->lchild = root_lchild;
    else
      n->rchild = root_lchild;

    root_lchild->parent = n;
    root_lchild->dparent = root_d;
  }
  /* now selected node is the root */

  /* if include_branch == TRUE, tweak the tree slightly (swap
     selected_node and its right child [separated by 0-length branch],
     adjust branch lengths beneath root) */
  if (include_branch) {
    TreeNode *root_lchild = selected_node->lchild;
    newroot = selected_node->rchild;

    selected_node->parent = newroot;
    selected_node->lchild = newroot->lchild;
    selected_node->rchild = newroot->rchild;

    newroot->rchild = selected_node;
    newroot->lchild = root_lchild;

    selected_node->dparent = root_lchild->dparent;
    root_lchild->dparent = 0;    
  }
  else                          /* otherwise new root is just selected node */
    newroot = selected_node;

  newroot->parent = NULL;
  newroot->dparent = 0;

  /* reset nnodes and height */
  for (i = 0; i < lst_size(tree->nodes); i++) {
    n = lst_get_ptr(tree->nodes, i);
    n->nnodes = (n == newroot ? lst_size(tree->nodes) : -1);
    n->height = 0;
    n->hold_constant = 0;
  }
  lst_free(tree->nodes);

  tr_set_nnodes(newroot);

  if (tree->preorder != NULL) {
    lst_free(tree->preorder);
    tree->preorder = NULL;
  }

  if (tree->postorder != NULL) {
    lst_free(tree->postorder);
    tree->postorder = NULL;
  }

  if (tree->inorder != NULL) {
    lst_free(tree->inorder);
    tree->inorder = NULL;
  }
}

/** Return an array indicating whether each node is in the designated
    subtree */
int* tr_in_subtree(TreeNode *t, TreeNode *sub) {
  int *in_subtree = smalloc(t->nnodes * sizeof(int));
  List *inside = lst_new_ptr(t->nnodes);
  int i;
  for (i = 0; i < t->nnodes; i++) in_subtree[i] = FALSE;
  tr_partition_nodes(t, sub, inside, NULL);
  for (i = 0; i < lst_size(inside); i++) {
    TreeNode *n = lst_get_ptr(inside, i);
    in_subtree[n->id] = TRUE;
  }
  lst_free(inside);
  return in_subtree;
}


void tr_label(TreeNode *t, const char *label) {
  if (t->label != NULL) sfree(t->label);
  t->label = copy_charstr(label);
}

void tr_label_node(TreeNode *tree, const char *nodename,
		   const char *label) { 
  TreeNode *node = tr_get_node(tree, nodename);
  if (node == NULL) die("ERROR: unknown node %s\n", nodename);
  tr_label(node, label);
}


void tr_label_subtree(TreeNode *tree, const char *subtreeNode, 
		      int include_leading_branch,
		      const char *label) {
  List *inside = lst_new_ptr(10);
  TreeNode *parent_node, *node;
  int i;
  parent_node = tr_get_node(tree, subtreeNode);
  if (parent_node == NULL) die("ERROR: unknown node %s\n", subtreeNode);
  tr_partition_nodes(tree, parent_node, inside, NULL);
  for (i=0; i < lst_size(inside); i++) {
    node = lst_get_ptr(inside, i);
    if (node != parent_node || include_leading_branch)
      tr_label(node, label);
  }
  lst_free(inside);
}


/** Sets rv to a list of nodes in tree with the given label */
void tr_get_labelled_nodes(TreeNode *tree, const char *label, List *rv) {
  List *traversal = tr_preorder(tree);
  TreeNode *node;
  int i;
  lst_clear(rv);
  for (i=0; i < lst_size(traversal); i++) {
    node = lst_get_ptr(traversal, i);
    if (node->label != NULL && strcmp(node->label, label)==0)
      lst_push_ptr(rv, node);
  }
}

/**
 * Given the optarg which should be a name of a file, it will open the .mod file
 * ignore everything except the tree and return a treeNode obeject.
 * 
 * @param string
 * @return parsedTree
 */
char* tr_only_from_file(const char* string){
    FILE* fin = phast_fopen(string,"r");
    /*Temporary holder for tags until tree is found*/
    char tag[STR_MED_LEN];
    String* tempStr = str_new(STR_LONG_LEN);
    int results;
           
    results = fscanf(fin,"%s",tag);
   
    /*Iterate until TREE is reached.*/
    while(strcmp(tag,"TREE:") && results != EOF){
        str_readline(tempStr,fin);
        results = fscanf(fin,"%s",tag);
    }
   
    if (!strcmp(tag,"TREE:")){
        str_readline(tempStr, fin);
        str_double_trim(tempStr);

        if (tempStr->chars[tempStr->length-1] == ';')
                tempStr->chars[--tempStr->length] = '\0';

     
    return tempStr->chars;
    
    }else
        return NULL;   
}
/* =====================================================================================*/
/**
 * Given a tree it will do a midpoint rooting. That is, consider the two branches on
 * each side, add their values and make the root the middle of these two distances based
 * on the longest branch. This function destroys the original root. The passed in node
 * should never be used again! Instead set your tree to the return value:
 * tree = midPointRooting(tree);
 * @param tree, root node of the tree to midpoint root. This will be permanently destroyed
 * and set NULL.
 * return, the new root to use for the tree.
 */
TreeNode* midpointRooting(TreeNode* tree){
  
  /*Find the most distant nodes and their distance.*/
  Pair myPair = findMostDistantNode(tree);
  TreeNode* newRoot;
  /*Parameters passed to @attempInsert() to know newLengths for children of newRoot.*/
  double newLeftDistance, newRightDistance;

  /* We need this total distance, the two species and their common ancestor.*/
  double maxDistance = myPair.distance;
  TreeNode* leaf1 = myPair.node1;
  TreeNode* leaf2 = myPair.node2;

  int commonAncestor = findCommonAncestor(tree, leaf1->id, leaf2->id);
  
  /* Let's start with leaf1 from here we shall check distances all the way up the tree
   * until the common ancestor to see if we insert the new root here.*/
  newRoot = attemptInsert(tree, leaf1, maxDistance, commonAncestor, &newLeftDistance,
          &newRightDistance);

  /* Check if we were successful! Otherwise we will succeed from the other side: */
  if(newRoot == NULL)
    newRoot = attemptInsert(tree, leaf2, maxDistance, commonAncestor, &newLeftDistance,
            &newRightDistance);

  /* From here we know for sure that we inserted a new node. Now newRoot->left holds the
   * parent, P and newRoot->right holds the child. */
  TreeNode* p = newRoot->lchild;
  TreeNode* c = newRoot->rchild;

  /*General case where the new node wasn't inserted between a child of the root and the
   root.*/
  if(p->id != tree->id){
    /*Go through all nodes up to the root flipping vertices.*/
    while(p->id != tree->id){
      if(p->lchild->id == c->id)
        p->lchild = p->parent;
      else
        p->rchild = p->parent;
      
      c = p;
      p = p->parent;
    }

    /* C holds the node on the side we came from that needs to link to the opposite side of
     * the root. Now we delete the old root and do this linking. */
    /*Left side case:*/
    if(c->lchild->id == tree->id){
      /* If we are here we connect this node to the left side of C, but are we on the
       right side of the root or left? Let's check! */
      TreeNode* correctChild = (tree->lchild->id == c->id) ? tree->rchild : tree->lchild;
      c->lchild = correctChild;
      correctChild->parent = c;
      /* Add distances to nodes.*/
      correctChild->dparent = c->dparent + correctChild->dparent;
    }
    else{ /*Right side case:*/
      /*If we are here we connect this node to the right side of C. */
      TreeNode* correctChild = (tree->lchild->id == c->id) ? tree->rchild : tree->lchild;
      c->rchild = correctChild;
      correctChild->parent = c;
      /* Add distances to nodes.*/
      correctChild->dparent = c->dparent + correctChild->dparent;
    }
    newRoot->lchild->parent = newRoot;
    newRoot->rchild->parent = newRoot;
    /* The tree is connected properly! But notice that every node knows who their parent is.
     * Yet when we inverted the root we didn't set these right, or the branch lengths. So we
     * do one last pass through the tree to fix this. */
    fixParentsAndLengths(newRoot);
    /*Set the lengths for the children of the new root. */
  newRoot->lchild->dparent = newLeftDistance;
  newRoot->rchild->dparent = newRightDistance;
  }
  /*Special Case where we inserted right below the root.*/
  else{
    /*Set the parent and lengths for the children of the new root. Notice this case needs
     * it done before we do work, the previous case needs it after... */
    newRoot->lchild->parent = newRoot;
    newRoot->rchild->parent = newRoot;
    newRoot->lchild->dparent = newLeftDistance;
    newRoot->rchild->dparent = newRightDistance;

    /*Figure out whether we are on the right side or the left side.*/
    TreeNode* correctChild = (p->lchild->id == c->id) ? tree->rchild : tree->lchild;
    /*We are on the left side. */
    newRoot->lchild = correctChild;
    /*Add up distances for final distance.*/
    correctChild->dparent = newLeftDistance + correctChild->dparent;
    correctChild->parent = newRoot;
  }

  /*Not sure the proper way to delete the root... but that's okay...*/
  free(tree);
  lst_set_ptr(newRoot->nodes, 0, newRoot);
  
  return newRoot;
}
/* =====================================================================================*/
/**
 * Helper function to @midPointRooting() should not be used otherwise! This function
 * fixes the messed up length of some children and nodes as explained in
 * @midPointRooting(). It does so by recusively descending through the tree making sure
 * the lengths and parents are set properly.
 * @param node: Tree to recurse through fixing lengths.
 * @return: void, operates on passed tree.
 */
void fixParentsAndLengths(TreeNode* node){
  /*Base Case: This is a leaf. Do nothing!*/
  if(node->lchild == NULL)
    return;

  /*Recurse on left and right sides!*/
  fixParentsAndLengths(node->rchild);
  fixParentsAndLengths(node->lchild);
  
  /* General Case. Attempt to fix node and recurse on children. For both sides, check if
   * our child knows we are their parent. */
  if(node->lchild->parent->id != node->id){
    node->lchild->parent = node;
    node->lchild->dparent = node->dparent;
  }
  /* Do same for the right side. */
  if(node->rchild->parent->id != node->id){
    node->rchild->parent = node;
    node->rchild->dparent = node->dparent;
  }
  
  return;
}
/* =====================================================================================*/
/**
 * Given the tree information it will iterate from the leaf until the common ancestor
 * attempting to insert the new root. If successful the new root is returned and the
 * lengths of the two children to the new root are permanently altered. We cannot know
 * ahead of time which nodes are changed (that's what this function computes!). Anyway
 * if you are here then that's what you wanted. If it fails to insert then we return NULL
 * and nothing is changed. Else we also populate newLeftDistance and newRightDistance with
 * the correct lengths.
 * @param root: The root of the entire tree.
 * @param leaf: Leaf to start attempting insertion from.
 * @param maxDistance: Distance from leaf to other leaf. This should be the total distance
 * as it will be halved.
 * @param commonAncestor: Common ancestor to attempt insertion up to.
 * @param newLeftDistance: will filled  with new length if successful in inserting root.
 * @param newRightDistance: will filled  with new length if successful in inserting root.
 * @return: new root if successful otherwise NULL.
 */
TreeNode* attemptInsert(TreeNode* root, TreeNode* leaf, double maxDistance,
        int commonAncestor, double* newLeftDistance, double* newRightDistance){
  TreeNode* newRoot = NULL;

  /* We need to insert the new root exactly half way. */
  double halfDistance = maxDistance / 2.0;
  /*Round off threshold.*/
  double epsilon = 0.00001;
  TreeNode* current;
  double distanceSoFar = 0;

  for(current = leaf;                /* Start at leaf2. */
      current->id != commonAncestor; /* Iterate until we hit the common ancestor. */
      current = current->parent){    /* Continue up the tree until the parent. */

    distanceSoFar += current->dparent;
    /*Two cases: 1: For sure the root goes here (distanceSoFar > halfDistance).
     * Case 2: The node falls exactly here (due to floating point arithmetic we use an
     * epsilon value  to check for this. 
     * Root goes between these two nodes!*/
    if( (distanceSoFar > halfDistance) ||
        (fabs(distanceSoFar - halfDistance) < epsilon) ){
      /* We know for sure this function call should create a new root! Create and init!*/
      newRoot = tr_new_node();
      initNewRoot(root, newRoot);

      /* We are in vertex (P, c) where P is the parent, and is the child node. We alway
       * let the left side of the root point to the parent and the right side to the
       * child. */
      newRoot->lchild = current->parent;
      newRoot->rchild = current;

       /* Compute distances and set them to nodes. Note: We got in here through either of
        * the above methods. Do the proper thing according to the case. */
      (*newLeftDistance) = (distanceSoFar > halfDistance) ?
        (distanceSoFar - halfDistance) :
        0.0001;
      
      (*newRightDistance) = halfDistance - distanceSoFar + current->dparent;
      /*We are done!*/
      break;
    }
  }

  return newRoot;
}
/* =====================================================================================*/
/**
 * Given a TreeNode, find the most distant pair of leaves in the whole tree.
 * We must find the distance between all leaves and pick the max. For n leafs
 * that is n*n different paths we have to compare. To calculate the distance
 * for (a,b) by summing up the distance from the a to the common ancestor
 * of a and b with the distance of b to the common ancestor. We then pick the
 * maximum. Note this is only meant to work on roots! Mainly cause of the use of the
 * preorder function.
 * @param root: Root of tree.
 * @return Pair containing the index of the two nodes furthest from each other.
 */
Pair findMostDistantNode(TreeNode* root){
  Pair maxPair;
  List* allNodeList = tr_preorder(root);
  
  int i, j;
    /*Create List with of nodes containing all children.*/
  List* childrenList = lst_new_ptr(EXPECTED_MAX_NODES);

  /*If the node is a child, add it's index to our childrenList.*/
  for(i = 0; i < lst_size(allNodeList); i++){
    TreeNode* n = lst_get_ptr(allNodeList, i);

    if(n->lchild == NULL)
      lst_push_ptr(childrenList, (void*)n);
  }
  /*Get distance and create distance matrix.*/
  int size = lst_size(childrenList);
  double childDistances[size][size];
  
  /*Initialize all nodes..*/
  for(i = 0; i < size; i++)
    for(j = 0; j < size; j++)
      childDistances[i][j] = -1;
  
  /*For all children find maximum distance to all other children. This is done by finding
   * the common ancestor of the two nodes and then adding up the lengths of the node up
   * to the common ancestor then adding these two values together. (This could be done
   * more efficiently, but meh. */
  for(i = 0; i < size; i++)
    for(j = 0; j < size; j++){
      /*Skip finding distance to oneself*/
      if(i == j) continue;
      /*Skip if we already found distance of [i]-[j] and we are now in [j]-[i]*/
      if(childDistances[i][j] != -1) continue;
      
      TreeNode* node1 = lst_get_ptr(childrenList, i);
      TreeNode* node2 = lst_get_ptr(childrenList, j);
      
      /*We find the common ancestor id. Now add the lengths from the node to the common
       ancestor.*/
      int commonAncestorId = findCommonAncestor(root, node1->id, node2->id);
      
      double sum1 = addUpToCommonAncestor(node1, commonAncestorId);
      double sum2 = addUpToCommonAncestor(node2, commonAncestorId);
      /*Set distances for nodes both ways.*/
      childDistances[i][j] = sum1 + sum2;
      childDistances[j][i] = sum1 + sum2;
    }
  
  /* Iterate until we find the max.*/
  double max = -1;
  for(i = 0; i < size; i++)
    for(j = 0; j < size; j++){
      /*Extract the proper ids and set equal to node.*/
      if(max < childDistances[i][j]){
        max = childDistances[i][j];

        maxPair.node1 = lst_get_ptr(childrenList, i);
        maxPair.node2 = lst_get_ptr(childrenList, j);
        maxPair.distance = max;
      }
    }
  
  /*Max pair now holds the max distance and correct nodes. :)*/
  return maxPair;
  }
/* =====================================================================================*/
/**
 * Given a node and a common ancestor of that node. It will iterate up adding values to
 * find the distance from our node to the ancestor.
 * @param node
 * @param commonAncestorId
 * @return 
 */
double addUpToCommonAncestor(TreeNode* node, int commonAncestorId){
  TreeNode* current;
  double sum = 0;
  
  /*Starting at current node, iterate up until common ancestor.*/
  for(current = node; current->id != commonAncestorId; current = current->parent){
    /*Error checking. Make sure we are in fact a child of our common ancestor.*/
    if(current == NULL){
      char* error =
      "Error: addUpToCommonAncestor(): node [%d] does not have node[%d]"
       "as a common ancestor!\n";
      die(error, node->id, commonAncestorId);
     }
      
    sum += current->dparent;
  }
  
  return sum;
}
/* =====================================================================================*/
/**
 * Given two leaf it will find the common ancestor between these two nodes and return it.
 * Every TreeNode has a List object called commonAncestor. This function those two passes
 * through the tree. First pass populates the list, second pass finds the set that
 * contains both nodes.
 * @param treeRoot: Root of our tree.
 * @param nodeIndex1: Child one to find common ancestor to with respect to child two.
 * @param nodeIndex2: Child two to find common ancestor to with respect to child two.
 * @return index of common ancestor, else -1 on failure.
 */
int findCommonAncestor(TreeNode* root, int nodeIndex1, int nodeIndex2){
  /*Error Checking...*/
  TreeNode* node1 = tr_get_node_id(root,nodeIndex1);
  TreeNode* node2 = tr_get_node_id(root,nodeIndex2);
  
  if(nodeIndex1 == nodeIndex2)
   die("ERROR: Invalid call to findCommonAncestor. Nodes [%d] are the same.\n",
           nodeIndex1);
  if(node1 == NULL || node2 == NULL){
    char* error =
    "ERROR: Invalid call to findCommonAncestor. Nodes [%d] [%d] not in tree.\n";
     die(error, nodeIndex1, nodeIndex2);
  }
  if(node1->lchild != NULL || node2->rchild != NULL){
    char* error =
    "ERROR: Invalid call to findCommonAncestor. Nodes [%d] or [%d] not leafs!\n";
    die(error, nodeIndex1, nodeIndex2);
  }
  
  /*First set all node's lists for their children. */
  populateCommonAncestor(root);
  
  int commonAncestor = findCommonAncestorHelper(root, nodeIndex1, nodeIndex2);
  return commonAncestor;
}
/* =====================================================================================*/
/**
 * Helper function for @findCommonAncestor. Does the work recursively. Finds and returns
 * the common ancestor for a tree that has ran @populateCommonAncestor.
 * @param currentNode: Current node we are recursing over. User should probably pass root.
 * @param nodeIndex1: Child one to find common ancestor to with respect to child two.
 * @param nodeIndex2: Child two to find common ancestor to with respect to child two.
 * @return index of Common ancestor.
 */
int findCommonAncestorHelper(TreeNode* node, int nodeIndex1, int nodeIndex2){
  /*Base case: We hit a leaf. Leaf cannot be common ancestor return -1.*/
  if(node->lchild == NULL)
    return -1;
  
  /*Recursive Case: Check if either children found common ancestor. If yes then pass the
   id up. Else try to see if we are the common ancestor! Note it must be done this way
   as if the child finds it first but we don't check for this first, we automatically
   already have both numbers in our childrenSet but we are not the common ancestor.*/
  int leftResult = findCommonAncestorHelper(node->lchild, nodeIndex1, nodeIndex2);
  int rightResult = findCommonAncestorHelper(node->rchild, nodeIndex1, nodeIndex2);
  
  /*Left or right node found common ancestor! Pass the news up the tree!*/
  if(leftResult != -1)
    return leftResult;
  if(rightResult != -1)
    return rightResult;
  
  /*Else we check if we are the common ancestor. Check if our list contains both nodes.*/
  if(lst_find_int(node->childList,nodeIndex1) != -1 &&
     lst_find_int(node->childList,nodeIndex2) != -1 ){
    /*We are the common ancestor. Let the world know...*/
    return node->id;
  }

  /*Else we are not... */
    return -1;
}
/* =====================================================================================*/
/**
 * Recurse through the tree setting all lists for common ancestors. This function
 * allocates memory that needs freeing!
 * Like always we assume the tree is properly formed.
 * @param node: node of tree to set @childrenSet for.
 * @return void: Sets the @childrenSet field for all nodes.
 */
void populateCommonAncestor(TreeNode* node){
  TreeNode* leftChild = node->lchild;
  TreeNode* rightChild = node->rchild;
  int i;
  
  /*Free old list and create a new one! */
  if(node->childList != NULL)
    lst_free(node->childList);
  /*Create a new list.*/
  node->childList = lst_new_int(EXPECTED_MAX_NODES);
      
  /*Base case: We hit a leaf set the list to just yourself.*/
  if(leftChild == NULL){
    lst_push_int(node->childList, node->id);
    return;
  }

  /*Recursive Case: Recurse on left and right child. This set it the union of left
   *and right children's sets. */
  populateCommonAncestor(leftChild);
  populateCommonAncestor(rightChild);
  
  /*Now set this set to the union of the two children sets.*/
  int leftSize = lst_size(leftChild->childList);
  int rightSize = lst_size(rightChild->childList);
  
  /*Iterate over children's lists to create the parent list as the union of the children's
   list.*/
  for(i = 0; i < leftSize; i++){
    int currentIndex = lst_get_int(leftChild->childList, i);
    lst_push_int(node->childList, currentIndex);
  }
  for(i = 0; i < rightSize; i++){
    int currentIndex = lst_get_int(rightChild->childList, i);
    lst_push_int(node->childList, currentIndex);
  }
  
  return;
}
/* =====================================================================================*/
/**
 * Given an allocated node. Populate by copying the important values of the original root.
 * @param root: Original Root of tree.
 * @param newRoot: new root to copy values to.
 * return: void.
 */
void initNewRoot(TreeNode* root, TreeNode* newRoot){
  newRoot->id = 0;
  newRoot->hold_constant = root->hold_constant;
  newRoot->dparent = root->dparent;
  newRoot->label = root->label;
  newRoot->nnodes = root->nnodes;
  newRoot->parent = root->parent;
  newRoot->lchild = root->lchild;
  newRoot->rchild = root->rchild;
  newRoot->nodes = root->nodes;
  
  return;
  }
/* =====================================================================================*/