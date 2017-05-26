/*
 */

#include <rhea_vis_adaptor.h>
#include <rhea_base.h> //TODO del

//TODO does this need ifdef?
//#ifdef USE_CATALYST
#if 1

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

/* global variables */
vtkCPProcessor* Processor = NULL;
vtkUnstructuredGrid* VTKGrid;

void
BuildVTKGrid(unsigned int order,
             unsigned int nPoints, const double *pointCoords,
             unsigned int nCells, unsigned int *cellToCoordIdx)
{
  /* set vtk cell type */
  VTKCellType cellType;
  vtkIdType nNodesPerCell;
  switch (order) {
  case 1:
    cellType = VTK_HEXAHEDRON;
    nNodesPerCell = 8;
    break;
  case 2:
    cellType = VTK_TRIQUADRATIC_HEXAHEDRON;
    nNodesPerCell = 27;
    break;
  default: /* unsupported order */
    RHEA_ABORT_NOT_REACHED (); //TODO do not use RHEA_ macro
  }

  /* create coordinates array */
  //TODO rather use SetTuple3 of (x,y,z) coordinates
  vtkNew<vtkDoubleArray> coordsArray;
  coordsArray->SetNumberOfComponents(3);
  coordsArray->SetArray(const_cast<double*>(pointCoords),
                        static_cast<vtkIdType>(nPoints*3),
                        1 /* do not delete */);

  /* create grid points */
  vtkNew<vtkPoints> points;
  points->SetData(coordsArray.GetPointer());
  VTKGrid->SetPoints(points.GetPointer());

  /* create cells */
  //TODO can maybe be optimized; don't use RHEA_ALLOC/RHEA_FREE
  vtkIdType* coordId = RHEA_ALLOC (vtkIdType, nNodesPerCell);
  VTKGrid->Allocate( static_cast<vtkIdType>(nCells * (nNodesPerCell + 1)) );
  for (unsigned int cell = 0; cell < nCells; cell++) {
    unsigned int *cellPoints = cellToCoordIdx + nNodesPerCell * cell;
    for (int node = 0; node < nNodesPerCell; node++) {
      coordId[node] = (vtkIdType) cellPoints[node];
    }
    VTKGrid->InsertNextCell(cellType, nNodesPerCell, coordId);
  }
  RHEA_FREE (coordId);
}

void
CreateArrays(unsigned int nPoints, unsigned int nCells)
{
  /* create velocity array */
  if (VTKGrid->GetPointData()->GetNumberOfArrays() == 0) {
    vtkNew<vtkDoubleArray> velocity;
    velocity->SetName("velocity");
    velocity->SetNumberOfComponents(3);
    velocity->SetNumberOfTuples(static_cast<vtkIdType>(nPoints));
    VTKGrid->GetPointData()->AddArray(velocity.GetPointer());
  }

  /* create pressure array */
  if (VTKGrid->GetCellData()->GetNumberOfArrays() == 0) {
    vtkNew<vtkDoubleArray> pressure;
    pressure->SetName("pressure");
    pressure->SetNumberOfComponents(1);
    VTKGrid->GetCellData()->AddArray(pressure.GetPointer());
  }
}

void
SetArrays(const double *velocityData, unsigned int nPoints,
          const double *pressureData, unsigned int nCells)
{
  /* fill velocity array */
  vtkDoubleArray* velocity =
    vtkDoubleArray::SafeDownCast(VTKGrid->GetPointData()->GetArray("velocity"));
//TODO del
//for (unsigned int i = 0; i < nPoints; i++) {
//  velocity->SetTypedTuple(i, &(velocityData[3*i]));
//}
  velocity->SetArray(const_cast<double*>(velocityData),
                     static_cast<vtkIdType>(nPoints*3), 1 /* do not delete */);

  /* fill pressure array
   * Note: The pressure array is a scalar array so we can reuse memory as long
   *       as we ordered the points properly. */
  vtkDoubleArray* pressure =
    vtkDoubleArray::SafeDownCast(VTKGrid->GetCellData()->GetArray("pressure"));
  pressure->SetArray(const_cast<double*>(pressureData),
                     static_cast<vtkIdType>(nCells), 1 /* do not delete */);
}

void
CreateVTKData(unsigned int order,
              unsigned int nPoints, const double *pointCoords,
              unsigned int nCells, unsigned int *cellToCoordIdx,
              const double *velocityData, const double *pressureData)
{
  /* build grid
   * Note: The grid structure isn't changing so we only build it the first time
   *       it's needed. If we needed the memory we could delete it and rebuild
   *       as necessary. */
  if (VTKGrid == NULL) {
    VTKGrid = vtkUnstructuredGrid::New();
    BuildVTKGrid(order, nPoints, pointCoords, nCells, cellToCoordIdx);
  }

  /* create arrays */
  CreateArrays(nPoints, nCells);

  /* set data of arrays */
  SetArrays(velocityData, nPoints, pressureData, nCells);
}

} /* end namespace */

void
rhea_vis_initialize (const char *catalyst_scripts[],
                     const int n_catalyst_scripts)
{
  /* return if nothing to do */
  if (n_catalyst_scripts <= 0) {
    return;
  }

  /* create vis processor */
  if (Processor == NULL) {
    Processor = vtkCPProcessor::New();
    Processor->Initialize();
  }
  else {
    Processor->RemoveAllPipelines();
  }

  /* pass scripts to vis processor */
  for (int i = 0; i < n_catalyst_scripts; i++) {
    vtkNew<vtkCPPythonScriptPipeline> pipeline;
    pipeline->Initialize(catalyst_scripts[i]);
    Processor->AddPipeline(pipeline.GetPointer());
  }
}

void
rhea_vis_finalize ()
{
  if (Processor != NULL) {
    Processor->Delete();
    Processor = NULL;
  }
  if (VTKGrid != NULL) {
    VTKGrid->Delete();
    VTKGrid = NULL;
  }
}

void
rhea_vis_process (unsigned int order,
                  unsigned int nPoints, const double *pointCoords,
                  unsigned int nCells, unsigned int *cellToCoordIdx,
                  const double *velocityData, const double *pressureData)
{
  vtkNew<vtkCPDataDescription> dataDescription;
  dataDescription->AddInput("input");

  /* force output */
//TODO del
//dataDescription->SetTimeData(time, timeStep);
//if (lastTimeStep == true) {
//  // assume that we want to all the pipelines to execute if it
//  // is the last time step.
//  dataDescription->ForceOutputOn();
//}
  dataDescription->ForceOutputOn();

  /* process visualization */
  if (Processor->RequestDataDescription(dataDescription.GetPointer()) != 0) {
    CreateVTKData(order, nPoints, pointCoords, nCells, cellToCoordIdx,
                  velocityData, pressureData);
    dataDescription->GetInputDescriptionByName("input")->SetGrid(VTKGrid);
    Processor->CoProcess(dataDescription.GetPointer());
  }
}

#endif /* USE_CATALYST */
