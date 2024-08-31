#ifndef ARGUMENTS_H
#define ARGUMENTS_H

Arguments *parse_args(int argc, char** argv, Arguments *args);
void check_bytes(void);
void print_args(Arguments *args, ulong *dims);
void print_args_recall(Arguments *args, ulong size);

#endif
