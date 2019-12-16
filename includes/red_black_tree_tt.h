#ifndef RED_BLACK_TREE_TT_H
#define RED_BLACK_TREE_TT_H

#include <doubly_linked_list_tt.h>

typedef enum{red_tt,black_tt,dblack_tt}color_tt;

typedef enum{left_tt,right_tt,root_tt}child_type_tt;

typedef struct rbnode_tt{

  color_tt node_color;
  struct rbnode_tt *lchild;
  struct rbnode_tt *rchild;
  struct rbnode_tt *parent;
  long key;
  void *value;
}rbnode_tt;

typedef struct rbtree_tt{

  rbnode_tt *root;
  void (*free_value)(void *);
  void (*write_value)(void *,void *,stream_type);
  void *(*read_value)(void *,stream_type);
  void *(*copy_value)(void *);
  long size;
}rbtree_tt;

int rbt_isLeaf(rbnode_tt *root_tt);
rbnode_tt *rbt_rotate_left(rbnode_tt **root_tt,rbnode_tt *pivot);
rbnode_tt *rbt_rotate_right(rbnode_tt **root_tt,rbnode_tt *pivot);
rbnode_tt *rbt_initializeNode(long strval, void *val, color_tt ncol, rbnode_tt *parent);
rbtree_tt rbt_initialize(void (*func)(void *),
                        void *(*funcRead)(void *,stream_type),
                        void (*funcWrite)(void *,void *,stream_type),
                        void *(*funcCopy)(void *));
int rbt_isEmpty(rbnode_tt *root_tt);
rbnode_tt *rbt_grandParent(rbnode_tt *rt);
rbnode_tt *rbt_uncle(rbnode_tt *rt);
rbnode_tt *rbt_insertBinary(rbtree_tt *rbt, rbnode_tt *root_tt, long d, void *v);
int rbt_isWhatChildType(rbnode_tt *rbn);
void rbt_insertNode(rbtree_tt *rbt, rbnode_tt *root_tt, long d, void *v);
void rbt_recolor(rbnode_tt **root_tt,rbnode_tt *inserted,rbnode_tt *gp,rbnode_tt *unc);
void rbt_restructure(rbnode_tt **root_tt,rbnode_tt *inserted,rbnode_tt *gp,rbnode_tt *unc);
rbnode_tt *rbt_sibling(rbnode_tt *n);
rbnode_tt *rbt_findNode(rbnode_tt *root_tt,long d);
void rbt_freeNode(rbtree_tt *rbt,rbnode_tt *node,int v);
rbnode_tt *rbt_deleteBinary(rbtree_tt *rbt, rbnode_tt *ins_pos,long data);
int rbt_numChildren(rbnode_tt *parent);
long rbt_deleteNode(rbtree_tt *rbt,rbnode_tt *root_tt,long d);
void rbt_recurDeletion(rbnode_tt **root_tt,rbnode_tt *u,rbnode_tt *s);
void rbt_transferKeyVal(rbtree_tt *rbt,rbnode_tt **source, rbnode_tt **dest);
void rbt_freeValues(rbtree_tt *rbt,rbnode_tt *curr);
void rbt_freeNodes(rbtree_tt *rbt,rbnode_tt *curr);
void rbt_free(rbtree_tt *rbt,rbnode_tt *st,int v);
rbnode_tt *rbt_predecessor(rbnode_tt *root_tt);
child_type_tt rbt_childType(rbnode_tt *rbn);
rbtree_tt rbt_copy(rbtree_tt *rbt);
void rbt_copySub(rbnode_tt **node,rbnode_tt *source,rbnode_tt *par,
                 void * (*cpy_func)(void *));

#endif
