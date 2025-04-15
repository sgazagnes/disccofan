void write_inertia_json(Node *tree, char *filename, LambdaVec *lvec, double *(*attribute)(void *), double lambda, ulong size);
void write_inertia_txt(Node *tree, char *filename, double *(*attribute)(void *), double lambda, ulong size, ulong *dims, int attrib);
/*Annexes */
void write_inertia_bin(Node *tree, char *filename, double *(*attribute)(void *), double lambda, ulong size, ulong *dims, ulong *attr_off, int attrib);
void compare_density(Node *tree, char *fname_in, char *fname_out,  double *(*attribute)(void *));
double *hessian(value *gvals, ulong *dims, ulong p);
void tree_seg_dir(Node *tree, value *out, bool *reached, ulong *rank, ulong lwb, ulong upb, double max_var, double (*attribute)(void *), double lambda);
void write_inertiaall_bin(Arguments *args,Node *tree, char *filename, double *(*attribute)(void *), double lambda, ulong size, ulong *dims, ulong *attr_off, int attrib, double *extrema);
void inertia_attributes_all(Arguments *args, Node *tree,  double *(*attribute)(void *), ulong *dims, ulong *attr_off, double *extrema);
void inertia_attributes_bins(Arguments *args, Node *tree,  double *(*attribute)(void *), ulong *dims, ulong *attr_off, double *extrema);
void inertia_attributes_hii(Arguments *args, Node *tree,  double *(*attribute)(void *), ulong *dims);
