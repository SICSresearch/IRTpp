////// SICS Cholesky algorithm matrix implementation for SICS compatible matrices.
#include <math.h>
#include <stdexcept>

void cholesky(Matrix<double> & A, Matrix<double> &L){
	// Start with first column in the matrix
	int col = A.nC();
	int row = A.nR();

	if(row != col){
		// Matrix must be square
		throw std::invalid_argument( "Matrix<double> cholesky :: Matrix must be square" );
	}
	if(A.nC() != L.nC()){
		// Matrix must be square
		throw std::invalid_argument( "Matrix<double> cholesky :: Matrix A and L must be of the same dimensions" );
	}
	if(A.nR() != L.nR()){
		// Matrix must be square
		throw std::invalid_argument( "Matrix<double> cholesky :: Matrix A and L must be of the same dimensions" );
	}
	//Matrix<double> L(col,row);
	L.reset();

	for (int i = 0; i < col; i++) {
		for (int j = 0; j < row; j++) {
			L(i,j) = 99;
			std::cout<<L(i,j)<<" ł ";
		}
	}std::cout<<std::endl;

	L.reset();

	std::cout<<" col row "<<col<<" "<<row<<std::endl;
	for (int i = 0; i < col; i++)
	{
		std::cout<<"i "<<i<<",";
		for (int k = i; i < row; k++)
		{
			std::cout<<"k "<<k<<std::endl;
			if( k == i)
			{
				int sum = 0;
				if(k>0)
				{
					for (int j = 0; j < (k-1); j++)
					 {
						sum += L(k,j)*L(k,j);
					}
				}
				L(k,k) = sqrt(A(k,k) - sum);
				std::cout<<"A("<<k<<","<<k<<") = "<<A(k,k)<<std::endl;
				std::cout<<"L("<<k<<","<<k<<") = "<<L(k,k)<<std::endl;
			}
			else
			{
				int sum = 0;
				if(i>0)
				{
					for (int j = 0; j < (i-1); j++)
					 {
						sum += L(i,j)*L(k,j);
					}
				}
				L(k,i) = ((A(k,i) - sum)/L(i,i));
				std::cout<<"A("<<k<<","<<i<<") = "<<A(k,i)<<std::endl;
				std::cout<<"L("<<k<<","<<i<<") = "<<L(k,i)<<std::endl;
			}

			std::cout<<"đ k i L"<<k<<" "<<i<<std::endl<<L<<A;
			std::cout<<"Una linea mas"<<std::endl;
		}
	}
	std::cout<<"Cholesky computation result : "<<std::endl;
	std::cout<<L<<std::endl;
}
