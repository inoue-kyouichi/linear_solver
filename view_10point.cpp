#include <stdio.h>
#include <stdlib.h>
#include <vtkSmartPointer.h>
#include <vtkTetra.h>
#include <vtkCellArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

int main(void)
{
     FILE *fp;
    char file_read_node[934] = "roundbar_node_10point.dat",
         file_read_volume[504] = "roundbar_element_vol_10point.dat",
         file_read_surface[234] = "roundbar_element_surf_10point.dat";

    int N = 934, E = 504, S = 234;
    double x[N], y[N], z[N];
    int    E_1[E], E_2[E], E_3[E], E_4[E], E_5[E], E_6[E], E_7[E], E_8[E], E_9[E], E_10[E],
           S_1[S], S_2[S], S_3[S], S_4[S], S_5[S], S_6[S], E_n[E], S_n[S];

    fp = fopen(file_read_node,"r");
    for (int i=0; i<N; i++){
        fscanf(fp, "%lf %lf %lf\n", &x[i], &y[i], &z[i]);
    }
    fclose(fp);

    fp = fopen(file_read_volume,"r");
    for (int i=0; i<E; i++){
        fscanf(fp, "%d %d %d %d %d %d %d %d %d %d %d\n", &E_n[i], &E_1[i], &E_2[i], &E_3[i], &E_4[i], &E_5[i], &E_6[i], &E_7[i], &E_8[i], &E_9[i], &E_10[i]);
    }
    fclose(fp);

    fp = fopen(file_read_surface,"r");
    for (int i=0; i<S; i++){
        fscanf(fp, "%d %d %d %d %d %d %d\n", &S_n[i], &S_1[i], &S_2[i], &S_3[i], &S_4[i], &S_5[i], &S_6[i]);
    }
    fclose(fp);

    cout<< E_1[0]<<endl;
    int sum = 0;
    //// Points ////
    auto Points =
    vtkSmartPointer<vtkPoints>::New();
    Points->SetDataTypeToFloat();
    Points->SetNumberOfPoints(934);
    for(int i=0; i<N; i++){
        Points->SetPoint(i, x[i], y[i], z[i]);
    }
    

    int sum0 = 0, sum1 = 1, sum2 = 2, sum3 = 3, sum4 = 4, sum5 = 5, sum6 = 6, sum7 = 7, sum8 = 8, sum9 = 9, sum10 = 10;
    //// Connectivity ////
    auto Connectivity =
    vtkSmartPointer<vtkIdTypeArray>::New();
    Connectivity->SetNumberOfTuples(5544);
    for (int i = 0; i<E; i++){
        
        Connectivity->SetComponent(sum0,0, 10);
        Connectivity->SetComponent(sum1,0, E_1[i]-1);
        Connectivity->SetComponent(sum2,0, E_2[i]-1);
        Connectivity->SetComponent(sum3,0, E_3[i]-1);
        Connectivity->SetComponent(sum4,0, E_4[i]-1);
        Connectivity->SetComponent(sum5,0, E_5[i]-1);
        Connectivity->SetComponent(sum6,0, E_6[i]-1);
        Connectivity->SetComponent(sum7,0, E_7[i]-1);
        Connectivity->SetComponent(sum8,0, E_8[i]-1);
        Connectivity->SetComponent(sum9,0, E_9[i]-1);
        Connectivity->SetComponent(sum10,0, E_10[i]-1); //Program numbers begin with one.
        sum0 = sum0+11;
        sum1 = sum1+11;
        sum2 = sum2+11;
        sum3 = sum3+11;
        sum4 = sum4+11;
        sum5 = sum5+11;
        sum6 = sum6+11;
        sum7 = sum7+11;
        sum8 = sum8+11;
        sum9 = sum9+11;
        sum10 = sum10+11;
    }

    

    //// CellArray ////
    auto CellArray =
    vtkSmartPointer<vtkCellArray>::New();
    CellArray->SetCells(504, Connectivity);

    //// Grid ////
    auto Grid =
    vtkSmartPointer<vtkUnstructuredGrid>::New();
    Grid->SetPoints(Points);
    Grid->SetCells(VTK_TETRA, CellArray);


    //// Writer ////
    auto Writer =
    vtkSmartPointer<vtkXMLUnstructuredGridWriter>
    ::New();
    Writer->SetInputData(Grid);
    Writer->SetFileName("roundbar.vtu");
    Writer->SetCompressorTypeToNone();
    Writer->SetDataModeToAscii();
    Writer->Write();

    return 0;
}