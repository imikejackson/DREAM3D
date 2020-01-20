/* ============================================================================
 * Copyright (c) 2020-2020 BlueQuartz Software, LLC
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice, this
 * list of conditions and the following disclaimer in the documentation and/or
 * other materials provided with the distribution.
 *
 * Neither the name of BlueQuartz Software, the US Air Force, nor the names of its
 * contributors may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 * USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "HierarchicalSmoothing.h"

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataContainerSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatFilterParameter.h"
#include "SIMPLib/FilterParameters/UInt64FilterParameter.h"
#include "SIMPLib/Geometry/TriangleGeom.h"

#include "SurfaceMeshing/SurfaceMeshingConstants.h"
#include "SurfaceMeshing/SurfaceMeshingVersion.h"
#include "SurfaceMeshing/SurfaceMeshingFilters/HierarchicalSmooth/VolumeSolver.h"

namespace
{
constexpr size_t k_TrianglesDimY = 3;
constexpr size_t k_VerticesDimY = 3;
constexpr size_t k_FaceLabelsDimY = 2;
constexpr size_t k_NodeTypesDimY = 1;
} // namespace

struct HierarchicalSmoothing::Impl
{
  UInt64ArrayType::WeakPointer m_TriList;
  FloatArrayType::WeakPointer m_VertexList;
  Int32ArrayType::ConstWeakPointer m_FaceLabelsList;
  Int32ArrayType::ConstWeakPointer m_NodeTypesList;

  Impl() = default;

  ~Impl() = default;

  Impl(const Impl&) = delete;
  Impl(Impl&&) = delete;
  Impl& operator=(const Impl&) = delete;
  Impl& operator=(Impl&&) = delete;

  void reset()
  {
    m_TriList.reset();
    m_VertexList.reset();
    m_FaceLabelsList.reset();
    m_NodeTypesList.reset();
  }
};

// -----------------------------------------------------------------------------
HierarchicalSmoothing::HierarchicalSmoothing()
: p_Impl(std::make_unique<Impl>())
, m_Iterations(53)
, m_Threshold(0.001f)
{
  initialize();

  setupFilterParameters();
}

// -----------------------------------------------------------------------------
HierarchicalSmoothing::~HierarchicalSmoothing() = default;

// -----------------------------------------------------------------------------
void HierarchicalSmoothing::initialize()
{
  clearErrorCode();
  clearWarningCode();
  setCancel(false);
}

// -----------------------------------------------------------------------------
void HierarchicalSmoothing::setupFilterParameters()
{
  FilterParameterVectorType parameters;

  parameters.push_back(SIMPL_NEW_UINT64_FP("Max Number of Iterations", Iterations, FilterParameter::Category::Parameter, HierarchicalSmoothing));

  parameters.push_back(SIMPL_NEW_FLOAT_FP("Threshold", Threshold, FilterParameter::Category::Parameter, HierarchicalSmoothing));

  {
    DataContainerSelectionFilterParameter::RequirementType req;
    req.dcGeometryTypes.push_back(IGeometry::Type::Triangle);
    parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("Input Geometry", DataContainerPath, FilterParameter::Category::RequiredArray, HierarchicalSmoothing, req));
  }

  {
    DataArraySelectionFilterParameter::RequirementType req =
        DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int32, k_FaceLabelsDimY, AttributeMatrix::Type::Face, IGeometry::Type::Triangle);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Face Labels", FaceLabelsPath, FilterParameter::Category::RequiredArray, HierarchicalSmoothing, req));
  }

  {
    DataArraySelectionFilterParameter::RequirementType req =
        DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int32, k_NodeTypesDimY, AttributeMatrix::Type::Vertex, IGeometry::Type::Triangle);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Node Types", NodeTypesPath, FilterParameter::Category::RequiredArray, HierarchicalSmoothing, req));
  }

  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
void HierarchicalSmoothing::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
void HierarchicalSmoothing::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  p_Impl->reset();

  if(!(m_Threshold > 0.0f))
  {
    setErrorCondition(-10, QObject::tr("Threshold value must be positive"));
  }

  auto dca = getDataContainerArray();

  if(dca == nullptr)
  {
    setErrorCondition(-11, QObject::tr("Failed to obtain DataContainerArray"));
    return;
  }

  TriangleGeom::Pointer triGeom = dca->getPrereqGeometryFromDataContainer<TriangleGeom, AbstractFilter>(this, m_DataContainerPath);

  if(triGeom == nullptr)
  {
    return;
  }

  auto triangles = triGeom->getTriangles();
  if(triangles == nullptr)
  {
    setErrorCondition(-12, QObject::tr("Failed to obtain Triangle DataArray"));
    return;
  }

  auto vertices = triGeom->getVertices();
  if(vertices == nullptr)
  {
    setErrorCondition(-13, QObject::tr("Failed to obtain Vertices DataArray"));
    return;
  }

  const std::vector<size_t> faceLabelsComponentDims{k_FaceLabelsDimY};
  const std::vector<size_t> nodeTypesComponentDims{k_NodeTypesDimY};

  auto faceLabels = dca->getPrereqArrayFromPath<Int32ArrayType, AbstractFilter>(this, m_FaceLabelsPath, faceLabelsComponentDims);
  if(faceLabels == nullptr)
  {
    return;
  }

  auto nodeTypes = dca->getPrereqArrayFromPath<Int32ArrayType, AbstractFilter>(this, m_NodeTypesPath, nodeTypesComponentDims);
  if(nodeTypes == nullptr)
  {
    return;
  }

  const std::vector<size_t> trianglesComponentDims{k_TrianglesDimY};
  const std::vector<size_t> verticesComponentDims{k_VerticesDimY};

  if(triangles->getComponentDimensions() == trianglesComponentDims)
  {
    setErrorCondition(-14, QObject::tr("Triangles data array has the wrong component dimensions."));
    return;
  }

  if(vertices->getComponentDimensions() == verticesComponentDims)
  {
    setErrorCondition(-15, QObject::tr("Vertices data array has the wrong component dimensions."));
    return;
  }

  if(triangles->getNumberOfTuples() != faceLabels->getNumberOfTuples())
  {
    setErrorCondition(-16, QObject::tr("Face labels data array and triangles data array must have the same number of tuples."));
    return;
  }

  if(vertices->getNumberOfTuples() != nodeTypes->getNumberOfTuples())
  {
    setErrorCondition(-17, QObject::tr("Node types data array and vertices data array must have the same number of tuples."));
    return;
  }

  p_Impl->m_TriList = triangles;
  p_Impl->m_VertexList = vertices;
  p_Impl->m_FaceLabelsList = faceLabels;
  p_Impl->m_NodeTypesList = nodeTypes;
}

// -----------------------------------------------------------------------------
void HierarchicalSmoothing::preflight()
{
  // These are the REQUIRED lines of CODE to make sure the filter behaves correctly
  setInPreflight(true);              // Set the fact that we are preflighting.
  emit preflightAboutToExecute();    // Emit this signal so that other widgets can do one file update
  emit updateFilterParameters(this); // Emit this signal to have the widgets push their values down to the filter
  dataCheck();                       // Run our DataCheck to make sure everthing is setup correctly
  emit preflightExecuted();          // We are done preflighting this filter
  setInPreflight(false);             // Inform the system this filter is NOT in preflight mode anymore.
}

// -----------------------------------------------------------------------------
void HierarchicalSmoothing::execute()
{
  initialize();
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  if(getCancel())
  {
    return;
  }

  auto triList = p_Impl->m_TriList.lock();
  if(triList == nullptr)
  {
    setErrorCondition(-18, QObject::tr("Failed to obtain %1").arg(triList->getDataArrayPath().serialize()));
    return;
  }
  auto vertexList = p_Impl->m_VertexList.lock();
  if(vertexList == nullptr)
  {
    setErrorCondition(-19, QObject::tr("Failed to obtain %1").arg(vertexList->getDataArrayPath().serialize()));
    return;
  }
  auto faceLabelList = p_Impl->m_FaceLabelsList.lock();
  if(faceLabelList == nullptr)
  {
    setErrorCondition(-20, QObject::tr("Failed to obtain %1").arg(faceLabelList->getDataArrayPath().serialize()));
    return;
  }
  auto nodeTypesList = p_Impl->m_NodeTypesList.lock();
  if(nodeTypesList == nullptr)
  {
    setErrorCondition(-21, QObject::tr("Failed to obtain %1").arg(nodeTypesList->getDataArrayPath().serialize()));
    return;
  }

  auto smoothedVertexList = FloatArrayType::CreateArray(vertexList->getNumberOfTuples(), vertexList->getComponentDimensions(), vertexList->getName(), true);
  if(smoothedVertexList == nullptr)
  {
    setErrorCondition(-22, QObject::tr("Failed to create output DataArray"));
    return;
  }

  {
    using namespace HierarchicalSmooth;

    auto triangles = Eigen::Map<TriMesh>(triList->data(), triList->getNumberOfTuples(), k_TrianglesDimY);
    auto vertices = Eigen::Map<const MeshNode>(vertexList->data(), vertexList->getNumberOfTuples(), k_VerticesDimY);
    auto faceLabels = Eigen::Map<const FaceLabel>(faceLabelList->data(), faceLabelList->getNumberOfTuples(), k_FaceLabelsDimY);
    auto nodeTypes = Eigen::Map<const NodeType>(nodeTypesList->data(), nodeTypesList->getNumberOfTuples(), k_NodeTypesDimY);

    auto logFunc = [this](const std::string& message) { notifyStatusMessage(QString::fromStdString(message)); };

    auto smoothedVertices = Eigen::Map<MeshNode>(smoothedVertexList->data(), smoothedVertexList->getNumberOfTuples(), k_VerticesDimY);

    hierarchicalSmooth(triangles, vertices, faceLabels, nodeTypes, smoothedVertices, m_Threshold, m_Iterations, logFunc);

    smoothedVertexList->copyIntoArray(vertexList);
  }
}

// -----------------------------------------------------------------------------
AbstractFilter::Pointer HierarchicalSmoothing::newFilterInstance(bool copyFilterParameters) const
{
  HierarchicalSmoothing::Pointer filter = HierarchicalSmoothing::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
QString HierarchicalSmoothing::getCompiledLibraryName() const
{
  return SurfaceMeshingConstants::SurfaceMeshingBaseName;
}

// -----------------------------------------------------------------------------
QString HierarchicalSmoothing::getBrandingString() const
{
  return SurfaceMeshingConstants::SurfaceMeshingBaseName;
}

// -----------------------------------------------------------------------------
QString HierarchicalSmoothing::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << SurfaceMeshing::Version::Major() << "." << SurfaceMeshing::Version::Minor() << "." << SurfaceMeshing::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
QString HierarchicalSmoothing::getGroupName() const
{
  return SIMPL::FilterGroups::Unsupported;
}

// -----------------------------------------------------------------------------
QString HierarchicalSmoothing::getSubGroupName() const
{
  return SurfaceMeshingConstants::SurfaceMeshingPluginDisplayName;
}

// -----------------------------------------------------------------------------
QString HierarchicalSmoothing::getHumanLabel() const
{
  return "Hierarchical Smoothing";
}

// -----------------------------------------------------------------------------
QUuid HierarchicalSmoothing::getUuid() const
{
  return QUuid("{7873b4b5-6b2a-551c-a683-2bdc10308ebf}");
}

// -----------------------------------------------------------------------------
HierarchicalSmoothing::Pointer HierarchicalSmoothing::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
HierarchicalSmoothing::Pointer HierarchicalSmoothing::New()
{
  struct make_shared_enabler : public HierarchicalSmoothing
  {
  };
  return std::make_shared<make_shared_enabler>();
}

// -----------------------------------------------------------------------------
QString HierarchicalSmoothing::getNameOfClass() const
{
  return ClassName();
}

// -----------------------------------------------------------------------------
QString HierarchicalSmoothing::ClassName()
{
  return QString("HierarchicalSmoothing");
}

// -----------------------------------------------------------------------------
void HierarchicalSmoothing::setIterations(uint64_t value)
{
  m_Iterations = value;
}

// -----------------------------------------------------------------------------
uint64_t HierarchicalSmoothing::getIterations() const
{
  return m_Iterations;
}

// -----------------------------------------------------------------------------
void HierarchicalSmoothing::setThreshold(float value)
{
  m_Threshold = value;
}

// -----------------------------------------------------------------------------
float HierarchicalSmoothing::getThreshold() const
{
  return m_Threshold;
}

// -----------------------------------------------------------------------------
void HierarchicalSmoothing::setFaceLabelsPath(const DataArrayPath& value)
{
  m_FaceLabelsPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath HierarchicalSmoothing::getFaceLabelsPath() const
{
  return m_FaceLabelsPath;
}

// -----------------------------------------------------------------------------
void HierarchicalSmoothing::setNodeTypesPath(const DataArrayPath& value)
{
  m_NodeTypesPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath HierarchicalSmoothing::getNodeTypesPath() const
{
  return m_NodeTypesPath;
}

// -----------------------------------------------------------------------------
void HierarchicalSmoothing::setDataContainerPath(const DataArrayPath& value)
{
  m_DataContainerPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath HierarchicalSmoothing::getDataContainerPath() const
{
  return m_DataContainerPath;
}
