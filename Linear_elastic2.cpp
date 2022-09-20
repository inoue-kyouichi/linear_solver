#include <iostream>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Geometry>

using namespace Eigen;
using namespace std;

typedef Matrix<double,6826*3,6826*3> MatrixNd;
typedef Matrix<double,6826,6826> Matrix6826d;
typedef Matrix<double,4,4> Matrix4d;
typedef Matrix<double,3,3> Matrix3d;
typedef Matrix<double,3,4> Matrix3_4d;
typedef Matrix<double,4,3> Matrix4_3d;
typedef Matrix<double,6826*3,1> Vector6826Nd;

const double mu = 0.3,  //[-]
             E = 205000.0, //[N/mm^2]
             F = 10000, //[N]
             W_tet = 1.0/4.0,
             W_tr = 1.0/3.0;

const double C[3][3] ={
                        {E/(1-2*mu), E*mu/((1+mu)*(1-2*mu)), E*mu/((1+mu)*(1-2*mu))},
                        {E*mu/((1+mu)*(1-2*mu)), E*(1-mu)/((1+mu)*(1-2*mu)), E*mu/((1+mu)*(1-2*mu)),},
                        {E*mu/((1+mu)*(1-2*mu)), E*mu/((1+mu)*(1-2*mu)), E*(1-mu)/((1+mu)*(1-2*mu)),}
                      },
             L1[3] = {0,0.5,0.5},
             L2[3] = {0.5,0,0.5},
             L3[3] = {0.5,0.5,0},
            // U[6826*3][6826*3] = {0};

int main()
{/*
    FILE *fp;
    char file_read_node[6826] = "node.dat",
         file_read_volume[35262] = "element.dat";
         //file_read_surface[234] = "roundbar_element_surf.dat";
/*
    int N = 6826, V = 35262, S = 234;
    double x[N], y[N], z[N];
    int   V_1[V], V_2[V], V_3[V], V_4[V],
          S_1[S], S_2[S], S_3[S];

    fp = fopen(file_read_node,"r");
    for (int i=0; i<N; i++){
        fscanf(fp, "%lf %lf %lf\n", &x[i], &y[i], &z[i]);
    }
    fclose(fp);

    fp = fopen(file_read_volume,"r");
    for (int i=0; i<V; i++){
        fscanf(fp, "%d %d %d %d\n", &V_1[i], &V_2[i], &V_3[i], &V_4[i]);
    }
    fclose(fp);
/*
    fp = fopen(file_read_surface,"r");
    for (int i=0; i<S; i++){
        fscanf(fp, "%d %d %d\n", &S_1[i], &S_2[i], &S_3[i]);
    }
    fclose(fp);*/
    
    //creating K
    Matrix3_4d NS,Nxi;
    Matrix3d J;
    Matrix4d ke;
    Matrix4_3d xe;
    //MatrixNd K;
    Vector4d Nx,Ny,Nz,N1,N2;
    int I[4] = {0};
    double detJ;
    NS << -1.0, 1.0, 0.0, 0.0,
          -1.0, 0.0, 1.0, 0.0,
          -1.0, 0.0, 0.0, 1.0;
    
    //K = MatrixNd::Zero();
/*
    for(int k = 0; k < V; k++){

        xe(0,0) = x[V_1[k]-1]; xe(0,1) = y[V_1[k]-1]; xe(0,2) = z[V_1[k]-1];
        xe(1,0) = x[V_2[k]-1]; xe(1,1) = y[V_2[k]-1]; xe(1,2) = z[V_2[k]-1];
        xe(2,0) = x[V_3[k]-1]; xe(2,1) = y[V_3[k]-1]; xe(2,2) = z[V_3[k]-1];
        xe(3,0) = x[V_4[k]-1]; xe(3,1) = y[V_4[k]-1]; xe(3,2) = z[V_4[k]-1];

        J = NS*xe;

        detJ = J.determinant();
        FullPivLU<Matrix3d> lu(J);
        Nxi = lu.solve(NS);
        cout<<detJ<<endl;
        for(int i; i<4; i++){
            Nx(i) = Nxi(0,i);
        }
        for(int i; i<4; i++){
            Ny(i) = Nxi(1,i);
        }
        for(int i; i<4; i++){
            Nz(i) = Nxi(2,i);
        }

        I[0] = V_1[k]-1; I[1] = V_2[k]-1; I[2] = V_3[k]-1; I[3] = V_4[k]-1;

        for(int i = 0; i<3; i++){
            switch(i){
                case 0: N1 = Nx; break;
                case 1: N1 = Ny; break;
                case 2: N1 = Nz; break;
            }

            for(int j = 0; j<3; j++){
                switch(j){
                    case 0: N2 = Nx; break;
                    case 1: N2 = Ny; break;
                    case 2: N2 = Nz; break; 
                }

                for(int p = 0; p<4; p++){
                    for(int q = 0; q<4; q++){
                        ke(p,q) = N1[p]*C[i][j]*N2[q]*detJ;

                        if(i == 0 && j == 0){
                            K(I[p],I[q]) += ke(p,q);
                        }
                        if(i == 0 && j == 1){
                            K(I[p],I[q]+N) += ke(p,q);
                        }
                        if(i == 0 && j == 2){
                            K(I[p],I[q]+2*N) += ke(p,q);
                        }
                        if(i == 1 && j == 0){
                            K(I[p]+N,I[q]) += ke(p,q);
                        }
                        if(i == 1 && j == 1){
                            K(I[p]+N,I[q]+N) += ke(p,q);
                        }
                        if(i == 1 && j == 2){
                            K(I[p]+N,I[q]+2*N) += ke(p,q);
                        }
                        if(i == 2 && j == 0){
                            K(I[p]+2*N,I[q]) += ke(p,q);
                        }
                        if(i == 2 && j == 1){
                            K(I[p]+2*N,I[q]+N) += ke(p,q);
                        }
                        if(i == 2 && j == 2){
                            K(I[p]+2*N,I[q]+2*N) += ke(p,q);
                        }
                    }
                }                
            }
        }
    }*/

    //creating G
   /* Vector3d A,B,cross;
    Matrix3d Ge,L;
    Matrix157d G157;
    MatrixNd G;
    double Area;
    int T[3] = {0};

    G157 = Matrix157d::Zero();
    G = MatrixNd::Zero();

    for(int i = 0; i < S; i++){
        A(0) = x[S_3[i]-1]-x[S_1[i]-1]; A(1) = y[S_3[i]-1]-y[S_1[i]-1]; A(2) = z[S_3[i]-1]-z[S_1[i]-1];
        B(0) = x[S_2[i]-1]-x[S_1[i]-1]; B(1) = y[S_2[i]-1]-y[S_1[i]-1]; B(2) = z[S_2[i]-1]-z[S_1[i]-1];
        cross = A.cross(B);
        Area = cross.norm();

        T[0] = S_1[i]-1; T[1] = S_2[i]-1; T[2] = S_3[i]-1;

        for(int p = 0; p<3; p++){
            for(int q = 0; q<3; q++){
                for(int m = 0; m<3; m++){
                    L(0,0) = pow(L1[m],2.0); L(0,1) = L1[m]*L2[m]; L(0,2) = L1[m]*L3[m];
                    L(1,0) = L1[m]*L2[m]; L(1,1) = pow(L2[m],2.0); L(1,2) = L2[m]*L3[m];
                    L(2,0) = L1[m]*L3[m]; L(2,1) = L2[m]*L3[m]; L(2,2) = pow(L3[m],2.0);

                    Ge(p,q) += W_tr*L(p,q)*Area;
                }

                G157(T[p],T[q]) += Ge(p,q);
            }
        }
    }
    for(int i = 0; i<N*3; i++){
        for(int j = 0; j<N*3; j++){
            if(i < N && j < N){
                G(i,j) = G157(i,j);
            }
            if(i >= N && i < 2*N){
                if(j >= N && j < 2*N){
                    G(i,j) = G157(i-N,j-N);
                }
            }
            if(i >= 2*N && j >= 2*N){
                G(i,j) = G157(i-2*N,j-2*N);
            }
        }
    }*/

    //boundary condition

    //displacement

    /*ofstream fout ("out.dat");
    if(fout.fail()){  
        cout << "出力ファイルをオープンできません" << endl;
    }
    fout<< G <<endl;*/
    
 /* for(int i = 0; i<V; i++){
            if(V_1[i] == 0){
            cout<< x[V_1[i]-1] <<endl;
            }
    }*/

    return 0;

}




    
