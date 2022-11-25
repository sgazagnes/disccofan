#ifndef COMMUNICATION_H
#define COMMUNICATION_H

void send_boundary(Boundary *b,  int dest);
Boundary *receive_boundary(int src, ulong size_item);
void send_updated_boundary(Boundary *b, int dest);
Boundary *receive_updated_boundary(Boundary *b, int src);
void send_boundary_par(Boundary *b, int dest);
Boundary *receive_boundary_par(int src);
Boundary *receive_boundary_att(int src, ulong size_item);
void send_updated_boundary_par(Boundary *b, int dest);
Boundary *receive_updated_boundary_par(Boundary *b, int src);
Boundary *receive_updated_boundary_att(Boundary *b, int src);
void send_updated_boundary_att(Boundary *b, int dest);
#endif
