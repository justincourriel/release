#ifndef RED_BLACK_TREE_H
#define RED_BLACK_TREE_H

#include <arb.h>
#include <doubly_linked_list.h>

typedef enum{red,black,dblack}color;

typedef enum{left,right,root}child_type;

typedef struct rbnode{

  color node_color;
  struct rbnode *lchild;
  struct rbnode *rchild;
  struct rbnode *parent;
  char *key;
  void *value;
}rbnode;

typedef struct rbtree{

  rbnode *root;
  void (*free_value)(void *);
  void (*write_value)(void *,void *,stream_type);
  void *(*read_value)(void *,stream_type);
  void *(*copy_value)(void *);
  long size;
}rbtree;

int isLeaf(rbnode *root);
int isRoot(rbnode *root);
rbnode *rotate_left(rbnode **root,rbnode *pivot);
rbnode *rotate_right(rbnode **root,rbnode *pivot);
rbnode *initializeRBNode(char *strval,void *val,color ncol,rbnode *parent);
rbtree initializeRBTree(void (*func)(void *),
                        void *(*funcRead)(void *,stream_type),
                        void (*funcWrite)(void *,void *,stream_type),
                        void *(*funcCopy)(void *));
int isEmpty(rbnode *root);
rbnode *grandParent(rbnode *rt);
rbnode *uncle(rbnode *rt);
rbnode *insertBinary(rbtree *rbt,rbnode *root,char *d,void *v);
char *printRBTreeToString(rbtree *rbt,rbnode *subtree,int depth);
char *printRBNodeToString(rbnode *node);
int testRBTree(rbnode *root);
int isWhatChildType(rbnode *rbn);
void insertRBNode(rbtree *rbt,rbnode *root,char *d,void *v);
void recolor(rbnode **root,rbnode *inserted,rbnode *gp,rbnode *unc);
void restructure(rbnode **root,rbnode *inserted,rbnode *gp,rbnode *unc);
rbnode *sibling(rbnode *n);
rbnode *findRBNode(rbnode *root, char *d);
rbnode *findRBNodeVal(void *val,rbnode *curr);
void freeRBNode(rbtree *rbt,rbnode *node,int v);
rbnode *deleteBinary(rbtree *rbt, rbnode *ins_pos, char *d1);
int numChildren(rbnode *parent);
long deleteRBNode(rbtree *rbt, rbnode *root, char *d1);
void recurDeletion(rbnode **root,rbnode *u,rbnode *s);
void transferKeyVal(rbtree *rbt,rbnode **source, rbnode **dest);
void freeRBTreeValues(rbtree *rbt,rbnode *curr);
void freeRBTreeNodes(rbtree *rbt,rbnode *curr);
void freeRBTree(rbtree *rbt,rbnode *st,int v);
rbnode *predecessor(rbnode *root);
long testEquivalentKeySets(rbtree *t1,rbtree *t2);
void testEquivalentKeySetsSub(rbnode *curr,rbtree *t2,long *test);
char *printRBNodeToStringVS(rbnode *node);
char *printRBTreeToStringVS(rbnode *root,int depth);
rbnode *fillRBDummyTree(int maxd,char *keys,color tc,int currd,
                         rbnode *parent);
int maxDepth(rbnode *root,int count);
void *keyVal(rbtree *rbt,char *key);
void layers(rbtree *rbt,rbnode *pos,dl_list *layersll,int currdep,
            long node_width);
char *NumSpaces(long n);
void PadString(char **ss,char pad,long nleft,long nright);
void fillRBTreeWithDummyNodes(rbtree *rbt,rbnode **pos,int depth,int maxd);
void fixSpacing(dl_list *layers,long dlev,long maxd,long node_width);
char *prettyPrintRBTreeToString(rbtree *rbt);
long largestKeySize(rbnode *current);
rbtree copyRBTree(rbtree *rbt);
void copyRBTreeSub(rbnode **node,rbnode *source,rbnode *par,
                   void * (*cpy_func)(void *));
void writeRBTreeToStream(rbtree *rbt,rbnode *current,void *fs,stream_type st);
void writeRBTreeToFile(rbtree *rbt,char *fname);
rbnode *readRBNode(rbtree *rbt, void *fs, stream_type st);
void readRBTreeFromFile(rbtree *rbt,char *fname);
void readRBTreeFromStream(rbtree *rbt, void *fs, stream_type st);
void mergeRBTrees(rbtree *t1, rbtree *t2, int ow, int addnew);
void mergeRBTreesSub(rbtree *t1, rbtree *t2, rbnode *start, int ow, int addnew);
child_type childType(rbnode *rbn);

#endif
