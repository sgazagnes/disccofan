#ifndef LAMBDAVEC_H
#define LAMBDAVEC_H

/******************************************************************************/
/*                             Lambda Vector handling                         */
/******************************************************************************/

LambdaVec *lambda_vector_create(uint size);
void lambda_vector_delete(LambdaVec *lvec);
LambdaVec *lambda_vector_read(char *prefix, char *suffixe, double imScale);

#endif
