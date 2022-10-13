#include <iostream>
#include <fstream>

using namespace std;

int main(){
    FILE *fp;
    char file_read_F[58] = "reaction_force6.dat",
         file_read_theoretical_F[58] = "theoretical_value.dat";

    int N = 58;
    double F[N], Ft[N];

    fp = fopen(file_read_F,"r");
    for (int i=0; i<N; i++){
        fscanf(fp, "%lf\n", &F[i]);
    }
    fclose(fp);

    fp = fopen(file_read_theoretical_F,"r");
    for (int i=0; i<N; i++){
        fscanf(fp, "%lf\n", &Ft[i]);
    }
    fclose(fp);
    
    double st = 0.0;
    ofstream fout ("plotdata.dat");
    for(int i = 0; i<N; i++){
        fout << st << " " << F[i] << " " << Ft[i] << endl;
        st += 0.00002;
    }
    return 0;
}