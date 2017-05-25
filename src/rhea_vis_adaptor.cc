/*
 */

#include <rhea_vis_adaptor.h>
#include <rhea_base.h>

#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPProcessor.h>
#include <vtkCPPythonScriptPipeline.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>

namespace
{
vtkCPProcessor* Processor = NULL;
vtkUnstructuredGrid* VTKGrid;

void BuildVTKGrid(unsigned int numberOfPoints, double* pointsData,
                  unsigned int numberOfCells, unsigned int* cellToCoordMap,
                  int order)
{
  const int           nNodesPerElem = (order + 1) * (order + 1) * (order + 1);
  VTKCellType         cellType;

  switch (order) {
  case 1:
    cellType = VTK_HEXAHEDRON;
    break;
  case 2:
    cellType = VTK_QUADRATIC_HEXAHEDRON;
    break;
  default: /* unsupported order */
    RHEA_ABORT_NOT_REACHED ();
  }

  // create the points information
  vtkNew<vtkDoubleArray> pointArray;
  pointArray->SetNumberOfComponents(3);
  pointArray->SetArray(pointsData, static_cast<vtkIdType>(numberOfPoints * 3),
                       1);
  vtkNew<vtkPoints> points;
  points->SetData(pointArray.GetPointer());
  VTKGrid->SetPoints(points.GetPointer());

  // create the cells
  vtkIdType* coordId = RHEA_ALLOC (vtkIdType, nNodesPerElem);
  VTKGrid->Allocate(static_cast<vtkIdType>(numberOfCells * (nNodesPerElem+1)));
  for (unsigned int cell = 0; cell < numberOfCells; cell++)
  {
    unsigned int* cellPoints = cellToCoordMap + nNodesPerElem * cell;
    for (int node = 0; node < nNodesPerElem; node++)
    {
      coordId[node] = (vtkIdType) cellPoints[node];
    }
    VTKGrid->InsertNextCell(cellType, nNodesPerElem, coordId);
  }
  RHEA_FREE (coordId);
}

void UpdateVTKAttributes(unsigned int numberOfPoints, double* velocityData,
                         unsigned int numberOfCells, double* pressureData)
{
  if (VTKGrid->GetPointData()->GetNumberOfArrays() == 0)
  {
    // velocity array
    vtkNew<vtkDoubleArray> velocity;
    velocity->SetName("velocity");
    velocity->SetNumberOfComponents(3);
    velocity->SetNumberOfTuples(static_cast<vtkIdType>(numberOfPoints));
    VTKGrid->GetPointData()->AddArray(velocity.GetPointer());
  }
  if (VTKGrid->GetCellData()->GetNumberOfArrays() == 0)
  {
    // pressure array
    vtkNew<vtkDoubleArray> pressure;
    pressure->SetName("pressure");
    pressure->SetNumberOfComponents(1);
    VTKGrid->GetCellData()->AddArray(pressure.GetPointer());
  }

  vtkDoubleArray* velocity =
    vtkDoubleArray::SafeDownCast(VTKGrid->GetPointData()->GetArray("velocity"));
  // The velocity array is ordered as
  //   vx0,vx1,vx2,..,vy0,vy1,vy2,..,vz0,vz1,vz2,..
  // so we need to create a full copy of it with VTK's ordering of
  //   vx0,vy0,vz0,vx1,vy1,vz1,..
  for (unsigned int i = 0; i < numberOfPoints; i++)
  {
    double values[3] = { velocityData[i], velocityData[i + numberOfPoints],
      velocityData[i + 2 * numberOfPoints] };
    velocity->SetTypedTuple(i, values);
  }

  vtkDoubleArray* pressure =
    vtkDoubleArray::SafeDownCast(VTKGrid->GetCellData()->GetArray("pressure"));
  // The pressure array is a scalar array so we can reuse
  // memory as long as we ordered the points properly.
  pressure->SetArray(pressureData, static_cast<vtkIdType>(numberOfCells), 1);
}

void BuildVTKDataStructures(unsigned int numberOfPoints, double* points,
                            unsigned int numberOfCells, unsigned int* cells,
                            int order,
                            double* velocity, double* pressure)
{
  if (VTKGrid == NULL)
  {
    // The grid structure isn't changing so we only build it
    // the first time it's needed. If we needed the memory
    // we could delete it and rebuild as necessary.
    VTKGrid = vtkUnstructuredGrid::New();
    BuildVTKGrid(numberOfPoints, points, numberOfCells, cells, order);
  }
  UpdateVTKAttributes(numberOfPoints, velocity, numberOfCells, pressure);
}
}

void CatalystInitialize(int numScripts, char* scripts[])
{
  if (Processor == NULL)
  {
    Processor = vtkCPProcessor::New();
    Processor->Initialize();
  }
  else
  {
    Processor->RemoveAllPipelines();
  }

  for (int i = 1; i < numScripts; i++)
  {
    vtkNew<vtkCPPythonScriptPipeline> pipeline;
    pipeline->Initialize(scripts[i]);
    Processor->AddPipeline(pipeline.GetPointer());
  }
}

void CatalystFinalize()
{
  if (Processor)
  {
    Processor->Delete();
    Processor = NULL;
  }
  if (VTKGrid)
  {
    VTKGrid->Delete();
    VTKGrid = NULL;
  }
}

void CatalystCoProcess(unsigned int numberOfPoints, double* pointsData,
                       unsigned int numberOfCells, unsigned int* cellsData,
                       int order,
                       double* velocityData, double* pressureData)
                       //double time, unsigned int timeStep, int lastTimeStep)
{
  vtkNew<vtkCPDataDescription> dataDescription;
  dataDescription->AddInput("input");
#if 0
  dataDescription->SetTimeData(time, timeStep);
  if (lastTimeStep == true)
  {
    // assume that we want to all the pipelines to execute if it
    // is the last time step.
    dataDescription->ForceOutputOn();
  }
#else
  dataDescription->ForceOutputOn();
#endif
  if (Processor->RequestDataDescription(dataDescription.GetPointer()) != 0)
  {
    BuildVTKDataStructures(numberOfPoints, pointsData,
                           numberOfCells, cellsData,
                           order,
                           velocityData, pressureData);
    dataDescription->GetInputDescriptionByName("input")->SetGrid(VTKGrid);
    Processor->CoProcess(dataDescription.GetPointer());
  }
}
