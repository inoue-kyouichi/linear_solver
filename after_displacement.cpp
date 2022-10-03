#include <iostream>
#include <fstream>

using namespace std;

int main(){
    FILE *fp;
    char file_read_node[6826] = "node.dat",
         file_read_volume[35262] = "element.dat",
         file_read_displacement[6826*3] = "displacement.dat";

    int N = 6826, V = 35262, D = 6826*3;
    double  x[N], y[N], z[N], d[D];
    int  V_1[V], V_2[V], V_3[V], V_4[V];
    
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

    fp = fopen(file_read_displacement,"r");
    for (int i=0; i<D; i++){
        fscanf(fp, "%lf\n", &d[i]);
    }
    fclose(fp);

    double x_result[N], y_result[N], z_result[N];

    for(int i = 0; i<N; i++){
        x_result[i] = d[i] + x[i];
        y_result[i] = d[i+N] + y[i];
        z_result[i] = d[i+2*N] + z[i];
    }

    ofstream fout ("out.dat");
    if(fout.fail()){  
        cout << "出力ファイルをオープンできません" << endl;
    }
    for(int i = 0; i<N; i++){
        fout<< x_result[i] << " " << y_result[i] << " " << z_result[i] << endl;
    }

}