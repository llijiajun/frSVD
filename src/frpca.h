// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <cstdlib>
#include <cmath>
#include <random>
#include <time.h>
#include <Rcpparmadillo.h>
#include <iostream>


using namespace std;
using namespace arma;

namespace frpca
{
	void QR_Decomposition(mat& Q)
	{

		static const double EPS(1E-4);
		for(std::size_t i = 0; i < Q.n_cols; ++i)
		{
			for(std::size_t j = 0; j < i; ++j)
			{
				double r = dot(Q.col(i),Q.col(j));
				Q.col(i) -= r * Q.col(j);
			}

			double nor = norm(Q.col(i),2);

			if(nor < EPS)
			{
				for(std::size_t k = i; k < Q.n_cols; ++k)
					Q.col(k) -= Q.col(k);
				return;
			}
			Q.col(i) /= nor;
		}
	}

	void LUfraction(mat& A,mat &result){
		mat L,U,P;
		lu(L,U,P,A);
		result=P*L;
	}
	void LUfraction(sp_mat& A,mat &result){
		mat L,U,P;
		mat B=mat(A);
		lu(L,U,P,B);
		result=P*L;
	}

	void eigSVD(mat& A,mat &U,vec &S,mat &V){
		V=A.t()*A;
		eig_sym(S,V,V);
		mat V1=V;
		std::size_t i,j;
//		cout<<"in eigSVD"<<endl;
//#pragma omp parallel shared(V1,S) private(i,j)
//    {
//    #pragma omp for
		for(i=0;i<V1.n_cols;i++){
			S(i)=sqrt(S[i]);
			for(j=0;j<V1.n_rows;j++)
				V1(j,i)/=S[i];
		}
//}
//		cout<<"out eigSVD"<<endl;
		U=A*V1;
	}
	void eigSVD(sp_mat &A,mat &U,vec &S,mat &V){
		V=A.t()*A;
		eig_sym(S,V,V);
		mat V1=V;
		std::size_t i,j;
//#pragma omp parallel shared(V1,S) private(i,j)
//    {
//    #pragma omp for
		for(i=0;i<V1.n_cols;i++){
			S(i)=sqrt(S[i]);
			for(j=0;j<V1.n_rows;j++)
				V1(j,i)/=S[i];
		}
//}
		U=A*V1;
	}

	void matrix_get_selected_cols(mat &M,uvec &inds,mat &Mc){
	  std::size_t i;
//		cout<<"getin"<<endl;
//		cout<<"M: " << M.n_rows<<M.n_cols<<endl;
//		cout<<"Mc: "<<Mc.n_rows<<Mc.n_cols<<endl;
//		cout<<inds<<endl;
		for(i=0;i<Mc.n_cols;i++){
			Mc.col(i)=M.col(inds(i));
		}
	}
	class frPCA
	{
	public:
		frPCA() {}
		frPCA(sp_mat& A, const int k, const int q)
		{
//			cout<<"kaishi"<<endl;
			int s=5;
			//umat location=randi<arma::umat>(2,A.n_rows);
			umat location=randi<arma::umat>(2,A.n_rows,distr_param(0,k+s-1));
			for(std::size_t i=0;i<A.n_rows;i++)
				location(0,i)=i;
			vec values = randi<vec>(A.n_rows,distr_param(0,1))*2-1;
			umat location2=randi<arma::umat>(2,A.n_cols,distr_param(0,k+s-1));
			for(std::size_t i=0;i<A.n_cols;i++)
				location2(0,i)=i;
			vec values2 = randi<vec>(A.n_cols,distr_param(0,1))*2-1;
			sp_mat Q(location,values),Qt(location2,values2);
			mat tempQ,tempQt,P;
			vec SS=zeros<mat>(k+s);
			mat VV=zeros<mat>(k+s,k+s);
			mat Qnext,Qtnext;
//			cout<<"yuliu"<<endl;
			int flag=0;
//			struct timeval start,end,st1,en1;
//			clock_t start2,end2,st2,en2;
//			gettimeofday(&start,NULL);
//			start2=clock();
			if(q%2==0){
				Qnext=A*Qt;
				if(q==2){
//					cout<<"q%2 eigSVD"<<endl;
					eigSVD(Qnext,Qnext,SS,VV);
				}else{
//					cout<<"q%2 LU"<<endl;
				//	LUfraction(Q,Qnext);
					lu(Qnext,U,Qnext);
		//			Qnext=P*Qnext;
				}
				flag=1;
//				cout<<"position2"<<endl;
			}
//			end2=clock();
//			gettimeofday(&end,NULL);
//			cout<<"step1:"<<get_seconds_frac(start,end)<<" thread:"<<(end2-start2)/(double)CLOCKS_PER_SEC<<endl;

			int niter=(q-1)/2,i;
//			cout<<"position3"<<endl;
			for(i=1;i<=niter;i++){
//			gettimeofday(&start,NULL);
//			start2=clock();
				if(flag==0){
//					cout<<"position3.5"<<endl;
//					gettimeofday(&st1,NULL);
//					st2=clock();
					Qtnext=A.t()*Q;
//					en2=clock();
//					gettimeofday(&en1,NULL);
//					cout<<"position3.6time:"<<get_seconds_frac(st1,en1)<<" thread:"<<(en2-st2)/(double)CLOCKS_PER_SEC<<endl;
//					gettimeofday(&st1,NULL);
//					st2=clock();
					//Q=A*Qt;
					Q=(Qtnext.t()*A.t()).t();
//					en2=clock();
//					gettimeofday(&en1,NULL);
//                                        cout<<"position4time:"<<get_seconds_frac(st1,en1)<<" thread:"<<(en2-st2)/(double)CLOCKS_PER_SEC<<endl;
					if(i==niter){
						eigSVD(Q,Qnext,SS,VV);
					}else{
				//		LUfraction(Q,Qnext);
						lu(Qnext,U,mat(Q));
		//				Qnext=P*Qnext;
//					cout<<"position5"<<endl;
					}
					flag=1;
				}else if(flag==1){
//					cout<<"po4.45"<<endl;
//                                        gettimeofday(&st1,NULL);
//					st2=clock();
					//Qtnext=A.t()*Qnext;
					Qtnext=(Qnext.t()*A).t();
//					en2=clock();
//                                        gettimeofday(&en1,NULL);
//                                        cout<<"position4.6time:"<<get_seconds_frac(st1,en1)<<" thread:"<<(en2-st2)/(double)CLOCKS_PER_SEC<<endl;
//                                        gettimeofday(&st1,NULL);
//					st2=clock();
					//Qnext=A*Qtnext;
					Qnext=(Qtnext.t()*A.t()).t();
//					en2=clock();
//                                        gettimeofday(&en1,NULL);
//					cout<<"position4.5time:"<<get_seconds_frac(st1,en1)<<" thread:"<<(en2-st2)/(double)CLOCKS_PER_SEC<<endl;
//					cout<<"position4.5"<<endl;
					//cout<<Q.n_rows<<","<<Q.n_cols<<endl;
//					gettimeofday(&st1,NULL);
//					st2=clock();
					if(i==niter){
						eigSVD(Qnext,Qnext,SS,VV);
					}else{
//					cout<<"po5.25"<<endl;
//						LUfraction(Qnext,Qnext);
						lu(Qnext,U,Qnext);
//						Qnext=P*Qnext;
					}
//					en2=clock();
//					gettimeofday(&en1,NULL);
//					cout<<"position5.5time:"<<get_seconds_frac(st1,en1)<<" thread:"<<(en2-st2)/(double)CLOCKS_PER_SEC<<endl;

				}
//			end2=clock();
//			gettimeofday(&end,NULL);
//			cout<<"step"<<i<<":"<<get_seconds_frac(start,end)<<" thread:"<<(end2-start2)/(double)CLOCKS_PER_SEC<<endl;
			}
			uvec inds(k);
			for(i=s;i<k+s;i++)
				inds(i-s)=i;
			S=zeros<vec>(k);

			//cout<<"INDS"<<k<<"k"<<S<<endl;
			if(flag==0){
//				cout<<"position6"<<endl;
//				gettimeofday(&start,NULL);
//				start2=clock();
				Qt=A.t()*Q;
//				end2=clock();
//				gettimeofday(&end,NULL);
//				cout<<"step out:"<<get_seconds_frac(start,end)<<" thread:"<<(end2-start2)/(double)CLOCKS_PER_SEC<<endl;
				mat UU;
//				gettimeofday(&start,NULL);
//				start2=clock();
				eigSVD(Qt,UU,SS,VV);
//				end2=clock();
//                                gettimeofday(&end,NULL);
//                                cout<<"step final:"<<get_seconds_frac(start,end)<<" thread:"<<(end2-start2)/(double)CLOCKS_PER_SEC<<endl;
				//cout<<"0SS"<<SS<<endl;
				mat VV2=zeros<mat>(k+s,k);
				matrix_get_selected_cols(UU,inds,V);
				matrix_get_selected_cols(VV,inds,VV2);
				for(i=0;i<k;i++)
					S(i)=SS[s+k-i-1];
				U=Q*VV2;
			}else{
//				cout<<"position7"<<endl;
				mat UU;
//				gettimeofday(&start,NULL);
//				start2=clock();
				Qtnext=A.t()*Qnext;
//                                gettimeofday(&end,NULL);
//				end2=clock();
//                                cout<<"step out:"<<get_seconds_frac(start,end)<<" thread:"<<(end2-start2)/(double)CLOCKS_PER_SEC<<endl;

//				gettimeofday(&start,NULL);
//				start2=clock();
				eigSVD(Qtnext,UU,SS,VV);
//				end2=clock();
//                                gettimeofday(&end,NULL);
//                                cout<<"step final:"<<get_seconds_frac(start,end)<<" thread:"<<(end2-start2)/(double)CLOCKS_PER_SEC<<endl;
				//cout<<"1SS:"<<SS<<endl;
				mat VV2=zeros<mat>(k+s,k);
//				cout<<"position8"<<endl;
//				cout<<"UU:"<<UU.n_rows<<","<<UU.n_cols<<endl;
				V=zeros<mat>(A.n_cols,k);
				matrix_get_selected_cols(UU,inds,V);
				V=V.t();
//				cout<<"V:"<<V.n_rows<<","<<V.n_cols<<endl;

				for(i=0;i<k;i++)
					S(i)=SS[s+k-i-1];
//				cout<<"position9"<<endl;
				matrix_get_selected_cols(VV,inds,VV2);
				U=Q*VV2;
//				cout<<"U:"<<U.n_rows<<","<<U.n_cols<<endl;
			}
//			cout<<"end"<<endl;
		}
		mat matrixU() const
		{
			return U;
		}

		vec singularValues() const
		{
			return S;
		}

		mat matrixV() const
		{
			return V;
		}

	private:
		mat U;
		vec S;
		mat V;
	};
}

