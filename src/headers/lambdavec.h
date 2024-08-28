#ifndef LAMBDAVEC_H
#define LAMBDAVEC_H

/******************************************************************************/
/*                             Lambda Vector handling                         */
/******************************************************************************/

LambdaVec *lambda_vector_create(uint size);
void lambda_vector_delete(LambdaVec *lvec);
LambdaVec *lambda_vector_read(char *name, double im_scale);
void print_lambda_vec(LambdaVec *lvec);
#endif
