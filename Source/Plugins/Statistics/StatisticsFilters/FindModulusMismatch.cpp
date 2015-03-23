/* ============================================================================
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "DREAM3DLib/FilterParameters/AbstractFilterParametersReader.h"
#include "DREAM3DLib/FilterParameters/AbstractFilterParametersWriter.h"

#include "FindModulusMismatch.h"

#ifdef DREAM3D_USE_PARALLEL_ALGORITHMS
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/partitioner.h>
#include <tbb/task_scheduler_init.h>
#endif

class FindModulusMismatchImpl
{
    int32_t* m_Labels;
    float* m_Moduli;
    float* m_DeltaModuli;

  public:
    FindModulusMismatchImpl(int32_t* labels, float* mod, float* deltaMod) :
      m_Labels(labels),
      m_Moduli(mod),
      m_DeltaModuli(deltaMod)
    {}

    virtual ~FindModulusMismatchImpl() {}

    void generate(size_t start, size_t end) const
    {
      int feature1, feature2;
      float deltaE;

      for (size_t i = start; i < end; i++)
      {
        feature1 = m_Labels[2 * i];
        feature2 = m_Labels[2 * i + 1];
        deltaE = 0.0f;
        if(feature1 > 0 && feature2 > 0)
        {
          deltaE = m_Moduli[feature1] - m_Moduli[feature2];
          if(deltaE < 0)
            deltaE = -deltaE;
        }
        m_DeltaModuli[i] = deltaE;
      }
    }

#ifdef DREAM3D_USE_PARALLEL_ALGORITHMS
    void operator()(const tbb::blocked_range<size_t>& r) const
    {
      generate(r.begin(), r.end());
    }
#endif
};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FindModulusMismatch::FindModulusMismatch()  :
  AbstractFilter(),
  m_ModuliArrayPath(DREAM3D::Defaults::DataContainerName, DREAM3D::Defaults::CellFeatureAttributeMatrixName, "ElasticModuli"),
  m_SurfaceMeshFaceLabelsArrayPath(DREAM3D::Defaults::DataContainerName, DREAM3D::Defaults::FaceAttributeMatrixName, DREAM3D::FaceData::SurfaceMeshFaceLabels),
  m_SurfaceMeshDeltaModulusArrayName("SurfaceMeshDeltaModulus"),
  m_ModuliArrayName("ElasticModuli"),
  m_Moduli(NULL),
  m_SurfaceMeshFaceLabelsArrayName(DREAM3D::FaceData::SurfaceMeshFaceLabels),
  m_SurfaceMeshFaceLabels(NULL)
{
  setupFilterParameters();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FindModulusMismatch::~FindModulusMismatch()
{
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindModulusMismatch::setupFilterParameters()
{
  FilterParameterVector parameters;
  parameters.push_back(FilterParameter::New("Required Information", "", FilterParameterWidgetType::SeparatorWidget, "", true));
  parameters.push_back(FilterParameter::New("Moduli", "ModuliArrayPath", FilterParameterWidgetType::DataArraySelectionWidget, getModuliArrayPath(), true, ""));
  parameters.push_back(FilterParameter::New("SurfaceMeshFaceLabels", "SurfaceMeshFaceLabelsArrayPath", FilterParameterWidgetType::DataArraySelectionWidget, getSurfaceMeshFaceLabelsArrayPath(), true, ""));
  parameters.push_back(FilterParameter::New("Created Information", "", FilterParameterWidgetType::SeparatorWidget, "", true));
  parameters.push_back(FilterParameter::New("SurfaceMeshDeltaModulus", "SurfaceMeshDeltaModulusArrayName", FilterParameterWidgetType::StringWidget, getSurfaceMeshDeltaModulusArrayName(), true, ""));
  setFilterParameters(parameters);
}
// -----------------------------------------------------------------------------
void FindModulusMismatch::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setSurfaceMeshDeltaModulusArrayName(reader->readString("SurfaceMeshDeltaModulusArrayName", getSurfaceMeshDeltaModulusArrayName() ) );
  setSurfaceMeshFaceLabelsArrayPath(reader->readDataArrayPath("SurfaceMeshFaceLabelsArrayPath", getSurfaceMeshFaceLabelsArrayPath() ) );
  setModuliArrayPath(reader->readDataArrayPath("ModuliArrayPath", getModuliArrayPath() ) );
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int FindModulusMismatch::writeFilterParameters(AbstractFilterParametersWriter* writer, int index)
{
  writer->openFilterGroup(this, index);
  DREAM3D_FILTER_WRITE_PARAMETER(FilterVersion)
  DREAM3D_FILTER_WRITE_PARAMETER(SurfaceMeshDeltaModulusArrayName)
  DREAM3D_FILTER_WRITE_PARAMETER(SurfaceMeshFaceLabelsArrayPath)
  DREAM3D_FILTER_WRITE_PARAMETER(ModuliArrayPath)
  writer->closeFilterGroup();
  return ++index; // we want to return the next index that was just written to
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindModulusMismatch::dataCheckVoxel()
{
  setErrorCondition(0);

  QVector<size_t> dims(1, 1);
  m_ModuliPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<float>, AbstractFilter>(this, getModuliArrayPath(), dims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if( NULL != m_ModuliPtr.lock().get() ) /* Validate the Weak Pointer wraps a non-NULL pointer to a DataArray<T> object */
  { m_Moduli = m_ModuliPtr.lock()->getPointer(0); } /* Now assign the raw pointer to data from the DataArray<T> object */
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindModulusMismatch::dataCheckSurfaceMesh()
{
  DataArrayPath tempPath;
  setErrorCondition(0);

  DataContainer::Pointer sm = getDataContainerArray()->getPrereqDataContainer<AbstractFilter>(this, getSurfaceMeshFaceLabelsArrayPath().getDataContainerName());
  if(getErrorCondition() < 0 || NULL == sm.get()) { return; }

  QVector<size_t> dims(1, 2);
  m_SurfaceMeshFaceLabelsPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<int32_t>, AbstractFilter>(this, getSurfaceMeshFaceLabelsArrayPath(), dims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if( NULL != m_SurfaceMeshFaceLabelsPtr.lock().get() ) /* Validate the Weak Pointer wraps a non-NULL pointer to a DataArray<T> object */
  { m_SurfaceMeshFaceLabels = m_SurfaceMeshFaceLabelsPtr.lock()->getPointer(0); } /* Now assign the raw pointer to data from the DataArray<T> object */
  dims[0] = 1;
  tempPath.update(m_SurfaceMeshFaceLabelsArrayPath.getDataContainerName(), m_SurfaceMeshFaceLabelsArrayPath.getAttributeMatrixName(), getSurfaceMeshDeltaModulusArrayName() );
  m_SurfaceMeshDeltaModulusPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>, AbstractFilter, float>(this,  tempPath, 180.0, dims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if( NULL != m_SurfaceMeshDeltaModulusPtr.lock().get() ) /* Validate the Weak Pointer wraps a non-NULL pointer to a DataArray<T> object */
  { m_SurfaceMeshDeltaModulus = m_SurfaceMeshDeltaModulusPtr.lock()->getPointer(0); } /* Now assign the raw pointer to data from the DataArray<T> object */
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindModulusMismatch::preflight()
{
  setInPreflight(true);
  emit preflightAboutToExecute();
  emit updateFilterParameters(this);
  dataCheckVoxel();
  dataCheckSurfaceMesh();
  emit preflightExecuted();
  setInPreflight(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindModulusMismatch::execute()
{
  setErrorCondition(0);
  dataCheckVoxel();
  if(getErrorCondition() < 0) { return; }
  dataCheckSurfaceMesh();
  if(getErrorCondition() < 0) { return; }

#ifdef DREAM3D_USE_PARALLEL_ALGORITHMS
  tbb::task_scheduler_init init;
  bool doParallel = true;
#endif

  int64_t numTriangles = m_SurfaceMeshFaceLabelsPtr.lock()->getNumberOfTuples();


#ifdef DREAM3D_USE_PARALLEL_ALGORITHMS
  if (doParallel == true)
  {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, numTriangles),
                      FindModulusMismatchImpl(m_SurfaceMeshFaceLabels, m_Moduli, m_SurfaceMeshDeltaModulus), tbb::auto_partitioner());
  }
  else
#endif
  {
    FindModulusMismatchImpl serial(m_SurfaceMeshFaceLabels, m_Moduli, m_SurfaceMeshDeltaModulus);
    serial.generate(0, numTriangles);
  }

  notifyStatusMessage(getHumanLabel(), "Completed");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer FindModulusMismatch::newFilterInstance(bool copyFilterParameters)
{
  FindModulusMismatch::Pointer filter = FindModulusMismatch::New();
  if(true == copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindModulusMismatch::getCompiledLibraryName()
{ return Statistics::StatisticsBaseName; }


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindModulusMismatch::getGroupName()
{ return DREAM3D::FilterGroups::StatisticsFilters; }


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindModulusMismatch::getSubGroupName()
{ return DREAM3D::FilterSubGroups::CrystallographicFilters; }


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindModulusMismatch::getHumanLabel()
{ return "Find Elastic Modulus Mismatch"; }

