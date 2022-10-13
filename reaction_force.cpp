//#include <iostream>
//#include <fstream>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Geometry>
#include <vector>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <bits/stdc++.h>
#include <cmath>

using namespace Eigen;
using namespace std;

typedef Triplet<double, int64_t> T;
typedef Matrix<double,3,4> Matrix3_4d;
typedef Matrix<double,4,3> Matrix4_3d;
typedef Matrix<double,6826*3,1> VectorNd;

const double nu = 0.3,  //[-]
             E = 205000000000.0, //[N/m^2] //iron
             Yp = 235000000.0, //[Pa] yield point
             W_tet = 1.0/4.0,
             W_tr = 1.0/3.0,
             r = 0.01, //[m]
             A = M_PI*pow(r,2.0),
             L = 0.1; //[m]
             //st = 0.00114;

const double L1[3] = {0.0,0.5,0.5},
             L2[3] = {0.5,0.0,0.5},
             L3[3] = {0.5,0.5,0.0};

void output_vtu(VectorNd F, double x[6826], double y[6826], double z[6826], int V_1[35262], int V_2[35262], int V_3[35262], int V_4[35262]) {
    int N_NODE = 6826, N_ELEM = 35262;
    ofstream ofs("reaction_force6.vtu");
    ofs << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
    ofs << "<UnstructuredGrid>\n";
    ofs << "<Piece NumberOfPoints= \"" << N_NODE << "\" NumberOfCells= \"" << N_ELEM << "\">\n";
    ofs << "<Points>\n";
    ofs << "<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(int i = 0; i < N_NODE; i++) {
        ofs << x[i] << " " << y[i] << " " << z[i] << "\n";
    }
    ofs << "</DataArray>\n";
    ofs << "</Points>\n";
    ofs << "<Cells>\n";
    ofs << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for(int i = 0; i < N_ELEM; i++) {
        ofs << V_1[i] << " " << V_2[i] << " " << V_3[i] << " " << V_4[i] <<"\n";
    }
    ofs << "</DataArray>\n";
    ofs << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for(int i = 0; i < N_ELEM; i++) {
        ofs << 4*(i+1) << "\n";
    }
    ofs << "</DataArray>\n";
    ofs << "<DataArray type=\"UInt32\" Name=\"types\" format=\"ascii\">";
    for(int i = 0; i < N_ELEM; i++) {
        ofs << "10\n";
    }
    ofs << "</DataArray>\n";
    ofs << "</Cells>\n";
    ofs << "<PointData>\n";
    ofs << "<DataArray type=\"Float64\" Name=\"force\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(int i = 0; i < N_NODE; i++) {
        ofs << F(i,0) << " " << F(i+N_NODE,0) << " " << F(i+2*N_NODE,0) << "\n";
    }
    ofs << "</DataArray>\n";
    ofs << "</PointData>\n";
    ofs << "</Piece></UnstructuredGrid></VTKFile>";
    ofs.close();
    return;
}

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
    vector<vector<double>> Ka(6826*3,vector<double>(6826*3,0.0));
    SparseMatrix<double> K(6826*3,6826*3), K1(6826*3,6826*3);
    int I[4] = {0};
    double detJ;
    double C[3][3][3][3] = {0.0};
    NS << -1.0, 1.0, 0.0, 0.0,
          -1.0, 0.0, 1.0, 0.0,
          -1.0, 0.0, 0.0, 1.0;

    C[0][0][0][0] = E*(1.0-nu)/((1.0+nu)*(1.0-2.0*nu)); C[0][0][1][1] = E*nu/((1.0+nu)*(1.0-2.0*nu)); C[0][0][2][2] = E*nu/((1.0+nu)*(1.0-2.0*nu));
    C[1][1][0][0] = E*nu/((1.0+nu)*(1.0-2.0*nu)); C[1][1][1][1] = E*(1.0-nu)/((1.0+nu)*(1.0-2.0*nu)); C[1][1][2][2] = E*nu/((1.0+nu)*(1.0-2.0*nu));
    C[2][2][0][0] = E*nu/((1.0+nu)*(1.0-2.0*nu)); C[2][2][1][1] = E*nu/((1.0+nu)*(1.0-2.0*nu)); C[2][2][2][2] = E*(1.0-nu)/((1.0+nu)*(1.0-2.0*nu));
    C[0][1][0][1] = E/(2.0*(1.0+nu)); C[0][2][0][2] = E/(2.0*(1.0+nu)); C[1][2][1][2] = E/(2.0*(1.0+nu));
    C[1][0][0][1] = E/(2.0*(1.0+nu)); C[0][1][1][0] = E/(2.0*(1.0+nu));
    C[2][0][0][2] = E/(2.0*(1.0+nu)); C[0][2][2][0] = E/(2.0*(1.0+nu));
    C[2][1][1][2] = E/(2.0*(1.0+nu)); C[1][2][2][1] = E/(2.0*(1.0+nu));
    C[1][0][1][0] = E/(2.0*(1.0+nu)); C[2][0][2][0] = E/(2.0*(1.0+nu)); C[2][1][2][1] = E/(2.0*(1.0+nu));

    J = Matrix3d::Zero(); 

    for(int t = 0; t < V; t++){

        xe(0,0) = x[V_1[t]]; xe(0,1) = y[V_1[t]]; xe(0,2) = z[V_1[t]];
        xe(1,0) = x[V_2[t]]; xe(1,1) = y[V_2[t]]; xe(1,2) = z[V_2[t]];
        xe(2,0) = x[V_3[t]]; xe(2,1) = y[V_3[t]]; xe(2,2) = z[V_3[t]];
        xe(3,0) = x[V_4[t]]; xe(3,1) = y[V_4[t]]; xe(3,2) = z[V_4[t]];

        J = NS*xe;

        detJ = J.determinant();
        FullPivLU<Matrix3d> lu(J);
        Nxi = lu.solve(NS);

        I[0] = V_1[t]; I[1] = V_2[t]; I[2] = V_3[t]; I[3] = V_4[t];

        for(int i = 0; i<3; i++){
            for(int k = 0; k<3; k++){
                ke = Matrix4d::Zero();

                for(int p = 0; p<4; p++){
                    for(int q = 0; q<4; q++){
                        for(int j = 0; j<3; j++){
                            for(int l = 0; l<3; l++){
                                ke(p,q) += Nxi(j,p)*C[i][j][k][l]*Nxi(l,q)*detJ/6.0;
                            }
                        }
                        if(i == 0 && k == 0){
                            Ka.at(I[p]).at(I[q]) += ke(p,q);
                        }
                        if(i == 0 && k == 1){
                            Ka.at(I[p]).at(I[q]+N) += ke(p,q);
                        }
                        if(i == 0 && k == 2){
                            Ka.at(I[p]).at(I[q]+2*N) += ke(p,q);
                        }
                        if(i == 1 && k == 0){
                            Ka.at(I[p]+N).at(I[q]) += ke(p,q);
                        }
                        if(i == 1 && k == 1){
                            Ka.at(I[p]+N).at(I[q]+N) += ke(p,q);
                        }
                        if(i == 1 && k == 2){
                            Ka.at(I[p]+N).at(I[q]+2*N) += ke(p,q);
                        }
                        if(i == 2 && k == 0){
                            Ka.at(I[p]+2*N).at(I[q]) += ke(p,q);
                        }
                        if(i == 2 && k == 1){
                            Ka.at(I[p]+2*N).at(I[q]+N) += ke(p,q);
                        }
                        if(i == 2 && k == 2){
                            Ka.at(I[p]+2*N).at(I[q]+2*N) += ke(p,q);
                        }
                    }
                }                
            }
        }
    }

    vector<T> tripletVec_K1;
    for(int i = 0; i<3*N; i++){
        for(int j = 0; j<3*N; j++){
            if( Ka.at(i).at(j) != 0 ){
                tripletVec_K1.push_back( T(i,j,Ka.at(i).at(j)) );
            }
        }
    }
    K1.setFromTriplets(tripletVec_K1.begin(), tripletVec_K1.end());

    //boundary condition
    for(int i = 0; i<108; i++){
        Ka.at(i).at(i) = 1.0;
        for(int j = 0; j<3*N; j++){
             if( i == j){
             }
             else{
                Ka.at(i).at(j) = 0.0;
             }
        }
    }
    for(int i = N; i<N+108; i++){
        Ka.at(i).at(i) = 1.0;
        for(int j = 0; j<3*N; j++){
             if( i == j){
             }
             else{
                Ka.at(i).at(j) = 0.0;
             }
        }
    }
    for(int i = 2*N; i<2*N+108; i++){
        Ka.at(i).at(i) = 1.0;
        for(int j = 0; j<3*N; j++){
             if( i == j){
             }
             else{
                Ka.at(i).at(j) = 0.0;
             }
        }
    }
    for(int i = 2*N+108; i<=2*N+215; i++){
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
             if( i == j){
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
   /* Vector3d A1, A2, A3, A4, B1, B2, B3, B4, cross1, cross2, cross3, cross4;
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
    G.setFromTriplets(tripletVec_G.begin(), tripletVec_G.end()); */

    //boundary condition
    //VectorNd u,F;
    //double p=0.0,sum;
    //u = VectorNd::Zero();
    //t(213+2*N,0) = F;   //t213_z is load point
    //ofstream fout ("out.dat");
    //for(int j=0; j<75; j++){
        //sum = 0;
       /* for(int i = 108; i<=215; i++){
            t(i+2*N,0) = p;
        }*/

        /*for(int i = 108; i<=215; i++){
            Gt(i,0) = 0.0;
            Gt(i+N,0) = 0.0;
        }*/

        //displacement
        /*SparseLU<SparseMatrix<double>> solver;
        solver.compute(K);
        u = solver.solve(Gt);*/

        //fout << u(144+2*N,0) << " " << u(178+2*N,0) << " " << u(182+2*N,0) << " " << u(186+2*N) << " " << u(213+2*N,0) <<endl;
    /*for(int i = 108; i<=215; i++){
        sum += Gt(i+2*N,0);
    }
    fout << sum << endl;
        p += Yp/75.0;
    }*/


    /*ofstream fout ("out.dat");
    if(fout.fail()){  
        cout << "出力ファイルをオープンできません" << endl;
    }
    fout<< u <<endl;*/

    //reaction force
    ofstream fout ("out.dat");
    double sum1, sum2, st = 0.0;
    VectorNd u,F;
    for(int j = 0; j<58 ;j++){
        sum1 = 0.0; sum2 = 0.0;
        u = VectorNd::Zero();
        F = VectorNd::Zero();
        
        for(int i = 0; i<108; i++){
            F(i,0) = 0.0;
            F(i+N,0) = 0.0;
            F(i+2*N,0) = 0.0;
        }
        for(int i = 108; i<=215; i++){
            F(i+2*N,0) = st*L;
        }

        SparseLU<SparseMatrix<double>> solver;
        solver.compute(K);
        u = solver.solve(F);

        F = K1*u;

        for(int i = 108; i<=215; i++){
            sum1 += F(i+2*N,0);
        }
        
        fout<<sum1<<endl;
        st += 0.00002;
    }

    /*ofstream fout ("out.dat");
    if(fout.fail()){  
        cout << "出力ファイルをオープンできません" << endl;
    }
    fout<< u <<endl;*/

   /*double sum1 = 0, sum2 = 0;
    for(int i = 108; i<=215; i++){
        sum1 += F(i+2*N,0);
    }
    for(int i = 0; i<108; i++){
        sum2 += F(i+2*N,0);
    }
    cout<<abs(sum1+sum2)<<endl;
    cout<<sum1<<endl;
    cout<<sum2<<endl;*/

    output_vtu(F,x,y,z,V_1,V_2,V_3,V_4);    

return 0;

}




    
