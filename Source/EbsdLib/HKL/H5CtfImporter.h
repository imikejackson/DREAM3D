/* ============================================================================
* Copyright (c) 2009-2015 BlueQuartz Software, LLC
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
* The code contained herein was partially funded by the followig contracts:
*    United States Air Force Prime Contract FA8650-07-D-5800
*    United States Air Force Prime Contract FA8650-10-D-5210
*    United States Prime Contract Navy N00173-07-C-2068
*
* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef _h5ctfimporter_h_
#define _h5ctfimporter_h_

#include "hdf5.h"

#include <QtCore/QVector>
#include <QtCore/QString>

#include "EbsdLib/EbsdLib.h"
#include "EbsdLib/EbsdSetGetMacros.h"
#include "EbsdLib/EbsdImporter.h"
#include "EbsdLib/HKL/CtfReader.h"
#include "EbsdLib/HKL/CtfPhase.h"

/**
 * @class H5CtfImporter H5CtfImporter.h EbsdLib/HKL/H5CtfImporter.h
 * @brief This class will read a series of .ctf files and store the values into
 * an HDF5 file according to the .h5ebsd specification
 *
 * @date March 23, 2011
 * @version 1.2
 *
 */
class EbsdLib_EXPORT H5CtfImporter : public EbsdImporter
{

  public:
    EBSD_SHARED_POINTERS(H5CtfImporter)
    EBSD_TYPE_MACRO_SUPER(H5CtfImporter, EbsdImporter)
    EBSD_STATIC_NEW_SUPERCLASS(EbsdImporter, H5CtfImporter)

    virtual ~H5CtfImporter();

    /**
     * @brief Imports a specific file into the HDF5 file
     * @param fileId The valid HDF5 file Id for an already open HDF5 file
     * @param index The slice index for the file
     * @param angFile The absolute path to the input .ang file
     */
    int importFile(hid_t fileId, int64_t index, const QString& angFile);

    /**
     * @brief Writes the phase data into the HDF5 file
     * @param reader Valid AngReader instance
     * @param gid Valid HDF5 Group ID for the phases.
     * @return error condition
     */
    int writePhaseData(CtfReader& reader, hid_t gid);

    /**
     * @brief Returns the dimensions for the EBSD Data set
     * @param x Number of X Voxels (out)
     * @param y Number of Y Voxels (out)
     */
    virtual void getDims(int64_t& x, int64_t& y);

    /**
     * @brief Returns the x and y resolution of the voxels
     * @param x The x resolution (out)
     * @param y The y resolution (out)
     */
    virtual void getResolution(float& x, float& y);

    /**
     * @brief Return the number of slices imported
     * @return
     */
    virtual int numberOfSlicesImported();

    /**
     * @brief This function sets the version of the H5Ebsd file that will be written.
     * @param version
     * @return
     */
    virtual void setFileVersion(uint32_t version);

  protected:
    H5CtfImporter();

    int writeSliceData(hid_t fileId, CtfReader& reader, int z, int actualSlice);

  private:
    int64_t xDim;
    int64_t yDim;
    int64_t zDim;
    float xRes;
    float yRes;
    float zRes;
    int m_NumSlicesImported;
    int   m_FileVersion;

    H5CtfImporter(const H5CtfImporter&); // Copy Constructor Not Implemented
    void operator=(const H5CtfImporter&); // Operator '=' Not Implemented
};

#endif /* _H5CTFIMPORTER_H_ */

