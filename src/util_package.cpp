#include "util_package.h"

using irtpp::weights;
using irtpp::quads;

//Please reference dat memory before this or set pointer to NULL.

irtpp::dataset* mat2dat(Rcpp::NumericMatrix mat)
{
  irtpp::dataset* dat;

  dat = new irtpp::dataset(0);

  for (int r = 0; r < mat.nrow(); r++)
  {
    std::vector<char> v(mat.ncol());

    for (int c = 0; c < mat.ncol(); c++)
    {
      v[c] = mat[c*mat.nrow()+r];
    }

    dat->size = mat.ncol();
    dat->push(v);
  }

  dat->size = mat.ncol();

  return dat;
}

Rcpp::NumericMatrix mat2rcpp(Matrix<double>* mat)
{
  Rcpp::NumericMatrix x(mat->nR(),mat->nC());

  for (int i = 0; i < x.nrow() ; i++)
  {
    for (int j = 0; j < x.ncol(); j++)
    {
      x[i+j*x.nrow()] = (*mat)(i,j);
    }
  }

  return x;
}

Rcpp::NumericVector getRVector(Matrix<double>* mat)
{
  Rcpp::NumericVector result(mat->nR()*mat->nC());

  for(int i = 0; i < mat->nR(); i++)
  {
    for(int j = 0; j < mat->nC(); j++)
    {
      result[i*mat->nC() + j] = (*mat)(i, j);
    }
  }

  return result;
}

//' uirtestimate
//' @param data -
//' @param model_ -
//' @return list
//' @export
// [[Rcpp::export]]
Rcpp::List uirtestimate(Rcpp::NumericMatrix data , int model_,double convergenceEpsilon)
{
  Matrix<double>*     args;
  Matrix<double>*     f;
  Matrix<double>*     r;
  Rcpp::NumericVector probability;
  double              loglikelihood;
  void**              status_list;
  irtpp::dataset*     d;
  irtpp::model*       m;

  loglikelihood = 0;
  d = mat2dat(data);

  if(model_ == 1)
  {
    irtpp::emestimation em(new irtpp::onepl(), d,convergenceEpsilon);
    status_list   = em.estimate();
    args          = em.coef();
    f             = em.getF();
    r             = em.getR();
    loglikelihood = em.LogLik();
    probability   = getRVector((Matrix<double>*) status_list[2]);
  }
  else if(model_ == 2)
  {
    irtpp::emestimation em(new irtpp::twopl(), d,convergenceEpsilon);
    status_list   = em.estimate();
    args          = em.coef();
    f             = em.getF();
    r             = em.getR();
    loglikelihood = em.LogLik();
    probability   = getRVector((Matrix<double>*) status_list[2]);
  }
  else
  {
    irtpp::emestimation em(new irtpp::threepl(), d,convergenceEpsilon);
    status_list   = em.estimate();
    args          = em.coef();
    f             = em.getF();
    r             = em.getR();
    loglikelihood = em.LogLik();
    probability   = getRVector((Matrix<double>*) status_list[2]);
  }

  double* returnpars;
  Rcpp::NumericVector pars(3*d->size);
  Rcpp::NumericVector ff(f->nR());
  Rcpp::NumericVector rr(r->nR() * r->nC());
  Rcpp::NumericVector theta(f->nR());
  Rcpp::NumericVector weightss(f->nR());

  returnpars = new double[3*d->size];

  for(int i = 0; i < args->nR(); i++)
  {
    for(int j = 0; j < args->nC(); j++)
    {
      returnpars[i*args->nC() + j] = (*args)(i, j);
    }
  }

  for(int i = 0 ; i < f->nR(); i ++ )
  {
    theta[i]    = quads(f->nR())[i];
    weightss[i] = weights(f->nR())[i];
  }

  for(int i=0; i< f->nR() ; i++) { ff[i] = (*f)(i, 0); }

  for(int i = 0; i < r->nR(); i++)
  {
    for(int j = 0; j < r->nC(); j++)
    {
      rr[i*r->nC() + j] = (*r)(i, j);
    }
  }

  //Now the quadnodes and weights
  if(model_ < 3) { for (int i = 2*d->size;i < 3*d->size;i++) { returnpars[i]=0; } }
  for (int i = 0; i < 3*d->size; i++) { pars[i] = returnpars[i]; }


  Rcpp::List result = Rcpp::List::create(
    Rcpp::_["z"]          = pars,
    Rcpp::_["iterations"] = *((int*)status_list[0]),
    Rcpp::_["LL"]         = loglikelihood,
    Rcpp::_["r"]          = rr,
    Rcpp::_["f"]          = ff,
    Rcpp::_["theta"]      = theta,
    Rcpp::_["weights"]    = weightss,
    Rcpp::_["prob_mat"]   = probability
  );

  delete (int*)status_list[0];
  delete (bool*)status_list[1];
  delete[] status_list;
  delete[] returnpars;

  return result;
}

//' abilityinterface
//'  @param zita_par -
//'  @param data -
//'  @param model_ -
//'  @param method -
//'  @param matrix_flag -
//'  @param prob_matrix -
//' @return list
//' @export
// [[Rcpp::export]]
Rcpp::List abilityinterface(Rcpp::NumericMatrix zita_par,
                            Rcpp::NumericMatrix data    ,
                            int                 model_  ,
                            int                 method  ,
                            bool                matrix_flag,
                            Rcpp::NumericVector prob_matrix)
{
  irtpp::model*   m;
  irtpp::dataset* d;
  Matrix<double>* z_temp;
  double**        result;

  d = mat2dat(data);
  // Now create the estimation
  irtpp::LatentTraitEstimation lte(d);

  if(model_ == 1)      { m = new irtpp::onepl(); }
  else if(model_ == 2) { m = new irtpp::twopl(); }
  else                 { m = new irtpp::threepl(); }

  z_temp = m->getZ(d->size);

  for (int i = 0; i < zita_par.ncol(); i++)
  {
    for (int j = 0; j < zita_par.nrow(); j++)
    {
      (*z_temp)(j, i) = zita_par(j, i);
    }
  }

  m->qnodes = 40;
  // Pass the model
  lte.setModel(m);

  // Ready to estimate
  if(matrix_flag)
  {
    m->probability = new Matrix<double>(m->qnodes, d->size);

    for(int i = 0; i < m->qnodes; i++)
    {
      for(int j = 0; j < d->size; j++)
      {
        (*(m->probability))(i, j) = prob_matrix[i * d->size + j];
      }
    }

    if(method == 0) { lte.estimateLatentTraitsEAP(); }
    else { lte.estimateLatentTraitsMAP(z_temp); }
  }
  else
  {
    if(method == 0) { lte.estimateLatentTraitsEAP(z_temp); }
    else { lte.estimateLatentTraitsMAP(z_temp); }
  }

  result = lte.lt->getListPatternTheta();

  // Return in list
  Rcpp::NumericVector pars1(lte.lt->pm->countItems() * lte.lt->pm->matrix.size());
  Rcpp::NumericVector pars_aux(1);
  Rcpp::NumericVector pars2(lte.lt->pm->matrix.size());

  for(unsigned int i = 0; i < lte.lt->pm->matrix.size(); i++)
  {
    for(int j = 0; j < d->size; j++)
    {
      pars1[i * d->size + j] = result[i][j];
    }

    pars2[i] = result[i][d->size];
  }

  lte.lt->deleteListPatternTheta(result);

  Rcpp::List z = Rcpp::List::create(Rcpp::_["patterns"] = pars1,
  Rcpp::_["trait"] = pars2,
  Rcpp::_["path"] = "No path");

  delete m;

  return z;
}

//' mapinterface
//'  @param zita_par -
//'  @param dat -
//'  @param e_model -
//'  @param matrix_flag -
//'  @param prob_matrix -
//' @return list
//' @export
//[[Rcpp::export]]
Rcpp::List mapinterface(Rcpp::NumericMatrix zita_par    ,
                        Rcpp::NumericMatrix dat         ,
                        int                 e_model     ,
                        bool                matrix_flag ,
                        Rcpp::NumericVector prob_matrix)
{
  Rcpp::List result = abilityinterface(zita_par   ,
                                       dat        ,
                                       e_model    ,
                                       1          ,
                                       matrix_flag,
                                       prob_matrix);

  return result;
}

//' eapinterface
//'  @param zita_par -
//'  @param dat -
//'  @param e_model -
//'  @param matrix_flag -
//'  @param prob_matrix -
//' @return list
//' @export
//[[Rcpp::export]]
Rcpp::List eapinterface(Rcpp::NumericMatrix zita_par   ,
                        Rcpp::NumericMatrix dat        ,
                        int                 e_model    ,
                        bool                matrix_flag,
                        Rcpp::NumericVector prob_matrix)
{
  Rcpp::List result = abilityinterface(zita_par   ,
                                       dat        ,
                                       e_model    ,
                                       0          ,
                                       matrix_flag,
                                       prob_matrix);

  return result;
}
