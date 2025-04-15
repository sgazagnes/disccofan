void build_tree_step(Arguments *args, Node *tree, ulong *dims, int step);
void tree_flood_tee_par(Node *tree, bool *visited, PrioQueue *q, ulong *ranks, ulong* ranks_inv, ulong *dims,  ulong lwb, ulong upb, int connectivity);
void tree_flood_tee_att(Node *tree, AuxDataStore *store,  bool *visited,  PrioQueue *q, ulong *ranks, ulong* ranks_inv,  ulong *dims,  ulong *attr_off, ulong lwb, ulong upb, int connectivity);
long tree_par_sal(Node *tree,  Queue *q,  idx *levelroot, bool *reached, ulong *dims, ulong lwb, ulong upb, long level, int connectivity);
void merge_nodes_par(Node *tree,idx x, idx y);
void merge_nodes_att(Node *tree, AuxDataStore *store, idx x, idx y) ;
