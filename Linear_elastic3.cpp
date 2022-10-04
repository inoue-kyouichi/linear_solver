#include <iostream>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Geometry>
#include <vector>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

using namespace Eigen;
using namespace std;

typedef Triplet<double> T;
typedef Matrix<double,4,4> Matrix4d;
typedef Matrix<double,3,3> Matrix3d;
typedef Matrix<double,3,4> Matrix3_4d;
typedef Matrix<double,4,3> Matrix4_3d;
typedef Matrix<double,6826*3,1> VectorNd;


const double mu = 0.3,  //[-]
             E = 205000000000.0, //[N/m^2] //iron
             //F = 10000000.0, //[N]
             F = 100.0,
             W_tet = 1.0/4.0,
             W_tr = 1.0/3.0;

const double C[6][6] ={
                        {E*(1.0-mu)/((1.0+mu)*(1.0-2.0*mu)), E*mu/((1.0+mu)*(1.0-2.0*mu)), E*mu/((1.0+mu)*(1.0-2.0*mu)), 0.0, 0.0, 0.0},
                        {E*mu/((1.0+mu)*(1.0-2.0*mu)), E*(1.0-mu)/((1.0+mu)*(1.0-2.0*mu)), E*mu/((1.0+mu)*(1.0-2.0*mu)), 0.0, 0.0, 0.0},
                        {E*mu/((1.0+mu)*(1.0-2.0*mu)), E*mu/((1.0+mu)*(1.0-2.0*mu)), E*(1.0-mu)/((1.0+mu)*(1.0-2.0*mu)), 0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0, E/(2.0*(1.0+mu)), 0.0, 0.0},
                        {0.0, 0.0, 0.0, 0.0, E/(2.0*(1.0+mu)), 0.0},
                        {0.0, 0.0, 0.0, 0.0, 0.0, E/(2.0*(1.0+mu))}
                      },
             L1[3] = {0.0,0.5,0.5},
             L2[3] = {0.5,0.0,0.5},
             L3[3] = {0.5,0.5,0.0};

int main()
{
    FILE *fp;
    char file_read_node[6826] = "node.dat",
         file_read_volume[35262] = "element.dat";
         //file_read_surface[234] = "roundbar_element_surf.dat";

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
    vector<vector<double>> Ka(6826*3,vector<double>(6826*3));
    Vector4d Nx,Ny,Nz,N1,N2;
    SparseMatrix<double> K(6826*3,6826*3);
    int I[4] = {0};
    double detJ;
    NS << -1.0, 1.0, 0.0, 0.0,
          -1.0, 0.0, 1.0, 0.0,
          -1.0, 0.0, 0.0, 1.0;

    J = Matrix3d::Zero();

    for(int k = 0; k < V; k++){

        xe(0,0) = x[V_1[k]]; xe(0,1) = y[V_1[k]]; xe(0,2) = z[V_1[k]];
        xe(1,0) = x[V_2[k]]; xe(1,1) = y[V_2[k]]; xe(1,2) = z[V_2[k]];
        xe(2,0) = x[V_3[k]]; xe(2,1) = y[V_3[k]]; xe(2,2) = z[V_3[k]];
        xe(3,0) = x[V_4[k]]; xe(3,1) = y[V_4[k]]; xe(3,2) = z[V_4[k]];

        J = NS*xe;

        detJ = J.determinant();
        FullPivLU<Matrix3d> lu(J);
        Nxi = lu.solve(NS);

        for(int i = 0; i<4; i++){
            Nx(i) = Nxi(0,i);
        }
        for(int i = 0; i<4; i++){
            Ny(i) = Nxi(1,i);
        }
        for(int i = 0; i<4; i++){
            Nz(i) = Nxi(2,i);
        }

        I[0] = V_1[k]; I[1] = V_2[k]; I[2] = V_3[k]; I[3] = V_4[k];

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
                            Ka.at(I[p]).at(I[q]) += ke(p,q);
                        }
                        if(i == 0 && j == 1){
                            Ka.at(I[p]).at(I[q]+N) += ke(p,q);
                        }
                        if(i == 0 && j == 2){
                            Ka.at(I[p]).at(I[q]+2*N) += ke(p,q);
                        }
                        if(i == 1 && j == 0){
                            Ka.at(I[p]+N).at(I[q]) += ke(p,q);
                        }
                        if(i == 1 && j == 1){
                            Ka.at(I[p]+N).at(I[q]+N) += ke(p,q);
                        }
                        if(i == 1 && j == 2){
                            Ka.at(I[p]+N).at(I[q]+2*N) += ke(p,q);
                        }
                        if(i == 2 && j == 0){
                            Ka.at(I[p]+2*N).at(I[q]) += ke(p,q);
                        }
                        if(i == 2 && j == 1){
                            Ka.at(I[p]+2*N).at(I[q]+N) += ke(p,q);
                        }
                        if(i == 2 && j == 2){
                            Ka.at(I[p]+2*N).at(I[q]+2*N) += ke(p,q);
                        }
                    }
                }                
            }
        }
    }
    //boundary condition
    for(int i = 0; i<=107; i++){
        Ka.at(i).at(i) = 1.0;
        for(int j = 0; j<3*N; j++){
             if( i == j){
             }
             else{
                Ka.at(i).at(j) = 0.0;
             }
        }
    }
    for(int i = N; i<=N+107; i++){
        Ka.at(i).at(i) = 1.0;
        for(int j = 0; j<3*N; j++){
             if( i == j){
             }
             else{
                Ka.at(i).at(j) = 0.0;
             }
        }
    }
    for(int i = 2*N; i<=2*N+107; i++){
        Ka.at(i).at(i) = 1.0;
        for(int j = 0; j<3*N; j++){
             if( i == j){
             }
             else{
                Ka.at(i).at(j) = 0.0;
             }
        }
    }
    /*for(int i = 108; i<=215; i++){
        Ka.at(i).at(i) = 1.0;
        for(int j = 0; j<3*N; j++){
            if(i == j){
            }
            else{
                Ka.at(i).at(j) = 0.0;
            }
        }
    }
    for(int i = N+108; i<=N+215; i++){
        Ka.at(i).at(i) = 1.0;
        for(int j = 0; j<3*N; j++){
            if(i == j){
            }
            else{
                Ka.at(i).at(j) = 0.0;
            }
        }
    }*/

    vector<T> tripletVec_K;
    for(int i = 0; i<3*N; i++){
        for(int j = 0; j<3*N; j++){
            if( Ka.at(i).at(j) != 0 ){
                tripletVec_K.push_back( T(i,j,Ka.at(i).at(j)) );
            }
        }
    }
    K.setFromTriplets(tripletVec_K.begin(), tripletVec_K.end());

    //creating G
    Vector3d A1, A2, A3, A4, B1, B2, B3, B4, cross1, cross2, cross3, cross4;
    Matrix3d Ge1,Ge2,Ge3,Ge4,L;
    vector<vector<double>> G157(6826, vector<double>(6826));
    SparseMatrix<double> G(6826*3,6826*3);
    double Area1, Area2, Area3, Area4;
    int T1[3], T2[3], T3[3], T4[3];

    for(int i = 0; i < V; i++){
        Ge1 = Matrix3d::Zero();
        Ge2 = Matrix3d::Zero();
        Ge3 = Matrix3d::Zero();
        Ge4 = Matrix3d::Zero();
        //triangle123
        A1(0) = x[V_3[i]]-x[V_1[i]]; A1(1) = y[V_3[i]]-y[V_1[i]]; A1(2) = z[V_3[i]]-z[V_1[i]];
        B1(0) = x[V_2[i]]-x[V_1[i]]; B1(1) = y[V_2[i]]-y[V_1[i]]; B1(2) = z[V_2[i]]-z[V_1[i]];
        cross1 = A1.cross(B1);
        Area1 = cross1.norm();
        T1[0] = V_1[i]; T1[1] = V_2[i]; T1[2] = V_3[i];
        //triangle124
        A2(0) = x[V_4[i]]-x[V_1[i]]; A2(1) = y[V_4[i]]-y[V_1[i]]; A2(2) = z[V_4[i]]-z[V_1[i]];
        B2(0) = x[V_2[i]]-x[V_1[i]]; B2(1) = y[V_2[i]]-y[V_1[i]]; B2(2) = z[V_2[i]]-z[V_1[i]];
        cross2 = A2.cross(B2);
        Area2 = cross2.norm();
        T2[0] = V_1[i]; T2[1] = V_2[i]; T2[2] = V_4[i];
        //triangle134
        A3(0) = x[V_4[i]]-x[V_1[i]]; A3(1) = y[V_4[i]]-y[V_1[i]]; A3(2) = z[V_4[i]]-z[V_1[i]];
        B3(0) = x[V_3[i]]-x[V_1[i]]; B3(1) = y[V_3[i]]-y[V_1[i]]; B3(2) = z[V_3[i]]-z[V_1[i]];
        cross3 = A3.cross(B3);
        Area3 = cross3.norm();
        T3[0] = V_1[i]; T3[1] = V_3[i]; T3[2] = V_4[i];
        //triangle234
        A4(0) = x[V_4[i]]-x[V_2[i]]; A4(1) = y[V_4[i]]-y[V_2[i]]; A4(2) = z[V_4[i]]-z[V_2[i]];
        B4(0) = x[V_3[i]]-x[V_2[i]]; B4(1) = y[V_3[i]]-y[V_2[i]]; B4(2) = z[V_3[i]]-z[V_2[i]];
        cross4 = A4.cross(B4);
        Area4 = cross4.norm();
        T4[0] = V_2[i]; T4[1] = V_3[i]; T4[2] = V_4[i];

        for(int p = 0; p<3; p++){
            for(int q = 0; q<3; q++){
                for(int m = 0; m<3; m++){
                    L(0,0) = pow(L1[m],2.0); L(0,1) = L1[m]*L2[m]; L(0,2) = L1[m]*L3[m];
                    L(1,0) = L1[m]*L2[m]; L(1,1) = pow(L2[m],2.0); L(1,2) = L2[m]*L3[m];
                    L(2,0) = L1[m]*L3[m]; L(2,1) = L2[m]*L3[m]; L(2,2) = pow(L3[m],2.0);

                    Ge1(p,q) += W_tr*L(p,q)*Area1;
                    Ge2(p,q) += W_tr*L(p,q)*Area2;
                    Ge3(p,q) += W_tr*L(p,q)*Area3;
                    Ge4(p,q) += W_tr*L(p,q)*Area4;
                }

                G157.at(T1[p]).at(T1[q]) += Ge1(p,q);
                G157.at(T2[p]).at(T2[q]) += Ge2(p,q);
                G157.at(T3[p]).at(T3[q]) += Ge3(p,q);
                G157.at(T4[p]).at(T4[q]) += Ge4(p,q);
            }
        }
    }
    vector<T> tripletVec_G;
    for(int i = 0; i<N*3; i++){
        for(int j = 0; j<N*3; j++){
            if( i < N && j < N){
                if(G157.at(i).at(j) != 0){
                    tripletVec_G.push_back( T(i,j,G157.at(i).at(j)) );
                }
            }
            if(i >= N && i < 2*N){
                if(j >= N && j < 2*N){
                    if(G157.at(i-N).at(j-N) != 0){
                        tripletVec_G.push_back( T(i,j,G157.at(i-N).at(j-N)) );
                    }
                }
            }
            if(i >= 2*N && j >= 2*N){
                if(G157.at(i-2*N).at(j-2*N) != 0){
                    tripletVec_G.push_back( T(i,j,G157.at(i-2*N).at(j-2*N)) );
                }
            }
        }
    }
    G.setFromTriplets(tripletVec_G.begin(), tripletVec_G.end()); 

    //boundary condition
    VectorNd u,t,Gt;
    t = VectorNd::Zero();
    Gt = VectorNd::Zero();
    u = VectorNd::Zero();
    //t(213+2*N,0) = F;   //t213_z is load point
    for(int i = 108; i<=215; i++){
        t(i+2*N,0) = F;
    }
    
    Gt = G*t;
    for(int i = 0; i<=107; i++){
        Gt(i,0) = 0.0;
        Gt(i+N,0) = 0.0;
        Gt(i+2*N,0) = 0.0;
    }
    /*for(int i = 108; i<=215; i++){
        Gt(i,0) = 0.0;
        Gt(i+N,0) = 0.0;
    }*/

    //displacement
    SparseLU<SparseMatrix<double>> solver;
    solver.compute(K);
    u = solver.solve(Gt);



    ofstream fout ("out.dat");
    if(fout.fail()){  
        cout << "出力ファイルをオープンできません" << endl;
    }
    fout<< u <<endl;

    /*for( int i = 0; i<4; i++){
        for(int j= 0; j<4; j++){
          if( k44_y(i,j) > 0){
            //cout<< Ka.at(i).at(j) << endl;
            cout << i << " " << j << endl;
          }
        }
    }*/
    /*for(int i = 0; i<N; i++){
        if( detJ != 0){
            cout << i << endl;
        }
    }*/

    return 0;

}




    
