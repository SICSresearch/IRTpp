#include <Rcpp.h>
 using namespace Rcpp;
 using std::cout;
 
 #include <iostream>
 #include "MultiPoly/dicho-multi/estimation/estimation.h"

void convert_matrix ( Rcpp::IntegerMatrix mat, irtpp::matrix<char> &Y ) {
    Y = irtpp::matrix<char>(mat.nrow(), mat.ncol());
    for ( int i = 0; i < mat.nrow(); ++i )
        for ( int j = 0; j < mat.ncol(); ++j )
            Y(i, j) = '0'+mat(i,j);
}

// [[Rcpp::export]]
Rcpp::List multiTest_dico(Rcpp::IntegerMatrix RDataset){
  irtpp::matrix<char> Y;
    convert_matrix(RDataset, Y);
	irtpp::dichomulti::estimation e(2, Y, 2, 0.001);
    e.EMAlgortihm();
	
	Rcpp::NumericVector a(e.data.d);
	Rcpp::NumericVector d(1);
	Rcpp::NumericVector c(1);
	
	for ( int i = 0; i < e.data.p; ++i ) {
		int j = 0;
		if ( e.data.m.parameters > 1 )
			for ( ; j < e.data.d; ++j ) a[i] = e.data.zeta[i](j);
		d[0] = e.data.zeta[i](j);
		if ( e.data.m.parameters == 3 ) c = e.data.zeta[i](++j);
		else c[0] = 0;
	}
	
	Rcpp::List output(3);
	output[0] = a;
	output[1] = d;
	output[2] = c;
	return output;
}

/*
Rcpp::List multiTest_dico(Rcpp::IntegerMatrix RDataset){
	matrix<char> Y;
    convert_matrix(RDataset, 4, 10, Y);
	dichomulti::estimation e(2, Y, 2, 0.001);
    e.EMAlgortihm();
	
	Rcpp::NumericVector a(e.data.d);
	Rcpp::NumericVector d(1);
	Rcpp::NumericVector c(1);
	
	for ( int i = 0; i < e.data.p; ++i ) {
		int j = 0;
		if ( e.data.m.parameters > 1 )
			for ( ; j < e.data.d; ++j ) cout << e.data.zeta[i](j) << ' '; //alphas
		cout << e.data.zeta[i](j) << ' '; //d
		if ( e.data.m.parameters == 3 ) cout << e.data.zeta[i](++j) << ' '; //c
	}
	
	
	///
	
	Rcpp::NumericVector a(e.data.d);
	Rcpp::NumericVector d(e.data.d - e.data.zeta);
	Rcpp::NumericVector c(1);
	
	for ( int i = 0; i < e.data.p; ++i ) {
		int j = 0;
		for ( ; j < e.data.d; ++j ) cout << e.data.zeta[i](j) << ' '; //alphas
		for ( ; j < e.data.zeta[i].size(); ++j ) cout << e.data.zeta[i](j) << ' '; //d's
	}
}
*/
