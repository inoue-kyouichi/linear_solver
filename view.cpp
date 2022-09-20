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
    char file_read_node[157] = "roundbar_node.dat",
         file_read_volume[504] = "roundbar_element_vol.dat",
         file_read_surface[234] = "roundbar_element_surf.dat";

    int N = 157, E = 504, S = 234;
    double x[N], y[N], z[N];
    int    E_1[E], E_2[E], E_3[E], E_4[E],
           S_1[S], S_2[S], S_3[S], E_n[E], S_n[S];

    fp = fopen(file_read_node,"r");
    for (int i=0; i<N; i++){
        fscanf(fp, "%lf %lf %lf\n", &x[i], &y[i], &z[i]);
    }
    fclose(fp);

    fp = fopen(file_read_volume,"r");
    for (int i=0; i<E; i++){
        fscanf(fp, "%d %d %d %d %d\n", &E_n[i], &E_1[i], &E_2[i], &E_3[i], &E_4[i]);
    }
    fclose(fp);

    fp = fopen(file_read_surface,"r");
    for (int i=0; i<S; i++){
        fscanf(fp, "%d %d %d %d\n", &S_n[i], &S_1[i], &S_2[i], &S_3[i]);
    }
    fclose(fp);

    int sum = 0;
    //// Points ////
    auto Points =
    vtkSmartPointer<vtkPoints>::New();
    Points->SetDataTypeToFloat();
    Points->SetNumberOfPoints(157);
    for(int i=0; i<N; i++){
        Points->SetPoint(i, x[i], y[i], z[i]);
    }
    

    int sum0 = 0, sum1 = 1, sum2 = 2, sum3 = 3, sum4 = 4;
    //// Connectivity ////
    auto Connectivity =
    vtkSmartPointer<vtkIdTypeArray>::New();
    Connectivity->SetNumberOfTuples(2520);
    for (int i = 0; i<E; i++){
        
        Connectivity->SetComponent(sum0,0, 4);
        Connectivity->SetComponent(sum1,0, E_1[i]-1);
        Connectivity->SetComponent(sum2,0, E_2[i]-1);
        Connectivity->SetComponent(sum3,0, E_3[i]-1);
        Connectivity->SetComponent(sum4,0, E_4[i]-1); //Program numbers begin with one.
        sum0 = sum0+5;
        sum1 = sum1+5;
        sum2 = sum2+5;
        sum3 = sum3+5;
        sum4 = sum4+5;
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