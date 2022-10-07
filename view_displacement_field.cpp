#include <stdio.h>
#include <stdlib.h>
#include <vtkSmartPointer.h>
#include <vtkTetra.h>
#include <vtkCellArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>

int main(void)
{
    FILE *fp;
    char file_read_node[6826] = "result.dat",
         file_read_volume[35262] = "element.dat",
         file_read_displacement[6826*3] = "displacement.dat";
    int N = 6826, E = 35262, S = 234, D = 6826*3;
    double x[N], y[N], z[N],D_1[D];
    int    E_1[E], E_2[E], E_3[E], E_4[E],
           S_1[S], S_2[S], S_3[S];


    fp = fopen(file_read_node,"r");
    for (int i=0; i<N; i++){
        fscanf(fp, "%lf %lf %lf\n", &x[i], &y[i], &z[i]);
    }
    fclose(fp);

    fp = fopen(file_read_volume,"r");
    for (int i=0; i<E; i++){
        fscanf(fp, "%d %d %d %d\n", &E_1[i], &E_2[i], &E_3[i], &E_4[i]);
    }
    fclose(fp);

    fp = fopen(file_read_displacement,"r");
    for (int i=0; i<D; i++){
        fscanf(fp, "%lf\n", &D_1[i]);
    }
    fclose(fp);

    int sum = 0;
    //// Points ////
    auto Points =
    vtkSmartPointer<vtkPoints>::New();
    Points->SetDataTypeToFloat();
    Points->SetNumberOfPoints(N);
    for(int i=0; i<N; i++){
        Points->SetPoint(i, x[i], y[i], z[i]);
    }
    

    int sum0 = 0, sum1 = 1, sum2 = 2, sum3 = 3, sum4 = 4;
    //// Connectivity ////
    auto Connectivity =
    vtkSmartPointer<vtkIdTypeArray>::New();
    Connectivity->SetNumberOfTuples(5*E);
    for (int i = 0; i<E; i++){
        
        Connectivity->SetComponent(sum0,0, 4);
        Connectivity->SetComponent(sum1,0, E_1[i]);
        Connectivity->SetComponent(sum2,0, E_2[i]);
        Connectivity->SetComponent(sum3,0, E_3[i]);
        Connectivity->SetComponent(sum4,0, E_4[i]); //Program numbers begin with one.
        sum0 = sum0+5;
        sum1 = sum1+5;
        sum2 = sum2+5;
        sum3 = sum3+5;
        sum4 = sum4+5;
    }

    //// CellArray ////
    auto CellArray =
    vtkSmartPointer<vtkCellArray>::New();
    CellArray->SetCells(E, Connectivity);

    //// Grid ////
    auto Grid =
    vtkSmartPointer<vtkUnstructuredGrid>::New();
    Grid->SetPoints(Points);
    Grid->SetCells(VTK_TETRA, CellArray);

     //// PointData ////
    auto PointData =
    vtkSmartPointer<vtkFloatArray>::New();
    PointData->SetNumberOfTuples(N);
    PointData->SetName("Point Scalar Data");
    for(int i=0; i<N; i++){
        PointData->SetComponent(i,0, D_1[i+2*N]);
    }
    Grid->GetPointData()->AddArray(PointData);

    //// PointData2 ////
    auto PointData2 =
    vtkSmartPointer<vtkFloatArray>::New();
    PointData2->SetNumberOfComponents(3);
    PointData2->SetNumberOfTuples(N);
    PointData2->SetName("Point Vector Data");
    for(int i=0; i<N; i++){
        PointData2->SetComponent(i,0, D_1[i]);
        PointData2->SetComponent(i,1, D_1[i+N]);
        PointData2->SetComponent(i,2, D_1[i+2*N]);
    }
    Grid->GetPointData()->AddArray(PointData2);

    //// CellData2 ////
    auto CellData2 =
    vtkSmartPointer<vtkFloatArray>::New();
    CellData2->SetNumberOfComponents(3);
    CellData2->SetNumberOfTuples(E);
    CellData2->SetName("Cell Vector Data");
    for(int i=0; i<E; i++){
        CellData2->SetComponent(i,0, 0.25*(D_1[E_1[i]+2*N]+D_1[E_2[i]+2*N]+D_1[E_3[i]+2*N]+D_1[E_4[i]+2*N]));
    }
    Grid->GetCellData()->AddArray(CellData2);


    //// Writer ////
    auto Writer =
    vtkSmartPointer<vtkXMLUnstructuredGridWriter>
    ::New();
    Writer->SetInputData(Grid);
    Writer->SetFileName("sample-c.vtu");
    Writer->SetCompressorTypeToNone();
    Writer->SetDataModeToAscii();
    Writer->Write();

    return 0;
}