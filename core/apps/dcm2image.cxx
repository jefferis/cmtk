/*
//
//  Copyright 2014, 2019 Google LLC
//
//  Copyright 1997-2012 Torsten Rohlfing
//
//  Copyright 2004-2014, 2022 SRI International
//
//  This file is part of the Computational Morphometry Toolkit.
//
//  http://www.nitrc.org/projects/cmtk/
//
//  The Computational Morphometry Toolkit is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  The Computational Morphometry Toolkit is distributed in the hope that it
//  will be useful, but WITHOUT ANY WARRANTY; without even the implied
//  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with the Computational Morphometry Toolkit.  If not, see
//  <http://www.gnu.org/licenses/>.
//
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#include <cmtkconfig.h>

#include <System/cmtkConsole.h>
#include <System/cmtkDebugOutput.h>
#include <System/cmtkCommandLine.h>
#include <System/cmtkStrUtility.h>
#include <System/cmtkExitException.h>

#include <Base/cmtkMetaInformationObject.h>

#include <IO/cmtkImageFileDICOM.h>
#include <IO/cmtkImageStackDICOM.h>
#include <IO/cmtkImageOperationXMLOptions.h>

#ifndef CMTK_USE_DCMTK
#error Build system is broken: this application should not be build if CMTK_USE_DCMTK is not set.
#endif

#ifdef _MSC_VER
#include <windows.h>
#else
#include <dirent.h>
#include <fnmatch.h>
#include <unistd.h>
#endif

#include <ctype.h>

#include <sys/types.h>
#include <sys/stat.h>

#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include <sstream>

#include <iostream>
#include <fstream>

#include <dcmtk/dcmdata/dctk.h>

void AddPatternToMap(std::map<DcmTagKey, std::string> &map, const char *pattern)
{
    const char *equals = strchr(pattern, '=');
    if (equals)
    {
        const std::string tagName(pattern, static_cast<size_t>(equals - pattern));

        const DcmDictEntry *dictEntry = dcmDataDict.rdlock().findEntry(tagName.c_str());
#if PACKAGE_VERSION_NUMBER < 364
        dcmDataDict.unlock();
#else
        dcmDataDict.rdunlock();
#endif

        if (dictEntry)
        {
            map.insert(std::pair<DcmTagKey, std::string>(dictEntry->getKey(), equals + 1));
        }
        else
        {
            cmtk::StdErr << "WARNING: DCMTK data dictionary does not have a tag entry named '" << tagName << "' - ignoring this\n";
        }
    }
}

// Map of file includion patterns
std::map<DcmTagKey, std::string> IncludePatterns;

void CallbackAddInclude(const char *pattern)
{
    AddPatternToMap(IncludePatterns, pattern);
}

// Map of file exclusion patterns
std::map<DcmTagKey, std::string> ExcludePatterns;

void CallbackAddExclude(const char *pattern)
{
    AddPatternToMap(ExcludePatterns, pattern);
}

const char *OutPathPattern = "image%n.nii";
std::vector<std::string> SearchRootDirVector;

std::ofstream cnull("/dev/null");

bool Progress = true;
bool Recursive = false;

/// Enum type to select sort order key.
typedef enum
{
    /// No sorting - take files in order stored in file system.
    SORT_NONE = 0,
    /// Sort by file name.
    SORT_FILENAME = 1,
    /// Sort by instance number.
    SORT_INSTANCE_NUMBER = 2
} SortKeyEnum;

SortKeyEnum SortFiles = SORT_INSTANCE_NUMBER;

bool WriteImage = true;
bool WriteXML = false;
bool StrictXML = false;
bool IncludeIdentifiers = false;
bool IncludeNDARIdentifiers = false;

bool WriteSingleSlices = false;

bool DisableOrientationCheck = false;
double Tolerance = 1e-5;

bool IgnoreAcquisitionNumber = false;

using cmtk::ImageFileDICOM;
using cmtk::ImageOperationXMLOptions;
using cmtk::ImageStackDICOM;

/// Selector for embedded image information.
ImageStackDICOM::EmbedInfoEnum EmbedInfo = ImageStackDICOM::EMBED_NONE;

/// Class handling a list of image stacks, i.e., volumes.
class VolumeList : public std::vector<ImageStackDICOM::SmartPtr>
{
public:
    /// This class.
    typedef VolumeList Self;

    /// Constructor.
    VolumeList(const cmtk::Types::Coordinate tolerance = 0 /*!< Tolerance for floating point comparisons, e.g., when testing for uniform pixel/slice spacings.*/) : m_Tolerance(tolerance) {}

    /// Build a volume list by traversing a directory.
    void AppendFromDirectory(const std::string &path, const char *wildcard);

    /** Write all volumes in this list.
     *\return Return code - 0, no problems; 3, at least one volume was not written correctly (e.g., due to non-uniform slice spacing).
     */
    int WriteVolumes();

private:
    /// Stored floating point tolerance.
    cmtk::Types::Coordinate m_Tolerance;

    /// Add a new image file to the correct stack or create a new stack.
    void AddImageFile(ImageFileDICOM::SmartConstPtr &newImage);
};

std::string
PutNumberAndSanitize(const std::string &path, const std::string &numberString)
{
    if (numberString.empty())
        return cmtk::StrReplace(cmtk::StrReplace(path, "%n", ""), "%N", "");
    else
        return cmtk::StrReplace(cmtk::StrReplace(path, "%n", numberString), "%N", "-" + numberString);
}

int VolumeList::WriteVolumes()
{
    int returnCode = 0;

    size_t cntSingleImages = 0;

    std::map<std::string, Self> pathToVolumeMap;
    for (const_iterator it = this->begin(); it != this->end(); ++it)
    {
        if (((*it)->size() > 1) || (*(*it)->begin())->m_IsMultislice || WriteSingleSlices)
        {
            // replace place holders
            std::string path(OutPathPattern);
            path = cmtk::StrReplace(path, "%D", cmtk::StrMakeLegalInPath((*it)[0][0]->GetTagValue(DCM_SeriesDescription)));
            path = cmtk::StrReplace(path, "%R", cmtk::StrMakeLegalInPath((*it)[0][0]->GetTagValue(DCM_RepetitionTime)));
            path = cmtk::StrReplace(path, "%E", cmtk::StrMakeLegalInPath((*it)[0][0]->GetTagValue(DCM_EchoTime)));
            path = cmtk::StrReplace(path, "%T", (*it)[0][0]->m_RawDataType);

            // replace %0 with file name
            path = cmtk::StrReplace(path, "%0", (*it)[0][0]->m_FileName);

            // and replace %1 through %9 with directory components of decreasing depth
            char replaceThis[3] = {'%', '0', 0};            // this will be modified to be the current replacement
            std::string directory = (*it)[0][0]->m_FileDir; // start with full directory
            for (size_t dirIdx = 1; (dirIdx < 10) && !directory.empty(); ++dirIdx)
            {
                replaceThis[1] = '0' + dirIdx;
                const size_t lastSlash = directory.rfind(CMTK_PATH_SEPARATOR);
                if (lastSlash == std::string::npos)
                {
                    path = cmtk::StrReplace(path, replaceThis, directory);
                    directory.clear();
                }
                else
                {
                    path = cmtk::StrReplace(path, replaceThis, directory.substr(lastSlash + 1));
                    directory = directory.substr(0, lastSlash);

                    // trim trailing slashes, both forward and back
                    const size_t lastNotSlash = directory.find_last_not_of("/\\");
                    if (lastNotSlash != std::string::npos)
                    {
                        directory.erase(lastNotSlash + 1);
                    }
                    else
                    {
                        directory.clear();
                    }
                }
            }

            if (path.length() > PATH_MAX)
                cmtk::StdErr << "ERROR: output path exceeds maximum path length";
            else
                pathToVolumeMap[path].push_back(*it);
        }
        else
        {
            ++cntSingleImages;
        }
    }

    if (cntSingleImages)
    {
        cmtk::DebugOutput(1) << "\n====================================================\n";
        cmtk::DebugOutput(1) << "WARNING: " << cntSingleImages << " single image(s) could not be assigned to multi-image stacks:\n\n";
        for (const_iterator it = begin(); it != end(); ++it)
        {
            if (((*it)->size() == 1) && !(*(*it)->begin())->m_IsMultislice)
            {
                (*(*it))[0]->Print();
                cmtk::DebugOutput(1) << "\n";
            }
        }
        cmtk::DebugOutput(1) << "\n====================================================\n";
    }

    for (std::map<std::string, Self>::const_iterator it = pathToVolumeMap.begin(); it != pathToVolumeMap.end(); ++it)
    {
        const size_t nVolumes = it->second.size();
        const ImageOperationXMLOptions xmlOptions(StrictXML, IncludeNDARIdentifiers, IncludeIdentifiers);

        // if there is only one volume with the given output path, just write it
        if (nVolumes == 1)
        {
            // if there's a "number" tag, get rid of it.
            std::string uniquePath = PutNumberAndSanitize(it->first, "");

            cmtk::UniformVolume::SmartConstPtr volume = it->second[0]->WriteImage(uniquePath.c_str(), EmbedInfo, !WriteImage);

            if (volume)
            {
                if (WriteXML)
                {
                    it->second[0]->WriteXML((uniquePath + ".xml").c_str(), *volume, xmlOptions);
                }
            }
            else
            {
                returnCode = 3;
            }
        }
        else
        {
            // otherwise, make unique paths for each of them
            for (size_t i = 0; i < nVolumes; ++i)
            {
                std::ostringstream numberString;
                numberString.width(1 + static_cast<int>(log((double)nVolumes) / M_LN10));
                numberString.fill('0');
                numberString << std::right << 1 + i;

                std::string uniquePath = PutNumberAndSanitize(it->first, numberString.str());

                cmtk::UniformVolume::SmartConstPtr volume = it->second[i]->WriteImage(uniquePath.c_str(), EmbedInfo, !WriteImage);

                if (volume)
                {
                    if (WriteXML)
                    {
                        it->second[i]->WriteXML((uniquePath + ".xml").c_str(), *volume, xmlOptions);
                    }
                }
                else
                {
                    returnCode = 3;
                }
            }
        }
    }

    return returnCode;
}

void VolumeList::AddImageFile(ImageFileDICOM::SmartConstPtr &newImage)
{
    if (!empty())
    {
        const_iterator it = this->begin();
        while (it != this->end())
        {
            ImageStackDICOM::SmartPtr study = *it;
            if (study->Match(*newImage, this->m_Tolerance, DisableOrientationCheck, IgnoreAcquisitionNumber))
            {
                study->AddImageFile(newImage);
                return;
            }
            else
            {
                ++it;
            }
        }
    }

    ImageStackDICOM::SmartPtr newImageStack(new ImageStackDICOM(this->m_Tolerance));
    newImageStack->AddImageFile(newImage);
    push_back(newImageStack);
}

void VolumeList::AppendFromDirectory(const std::string &path, const char *wildcard)
{
    std::vector<ImageFileDICOM::SmartConstPtr> fileList;

#ifdef _MSC_VER
    static int progress = 0;
    const char progress_chars[] = "-\\|/";

    WIN32_FIND_DATA fData;
    const std::string pattern = path + "\\" + wildcard;
    HANDLE hFind = FindFirstFile(pattern.c_str(), &fData);
    do
    {
        const std::string fullname = path + "\\" + fData.cFileName;
        if (fData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
        {
            if (Recursive && (fData.cFileName[0] != '.'))
            {
                this->AppendFromDirectory(fullname, wildcard);
            }
        }
        else
        {
            try
            {
                ImageFileDICOM::SmartPtr newFile(new ImageFileDICOM(fullname));

                if (newFile->MatchAllPatterns(IncludePatterns) && !newFile->MatchAnyPattern(ExcludePatterns))
                {
                    newFile->ReleaseDocument();
                    fileList.push_back(newFile);
                }

                if (Progress)
                {
                    (cmtk::StdErr << "\r" << progress_chars[++progress % 4]).flush();
                }
            }
            catch (...)
            {
                // not a valid DICOM file
            }
        }
    } while (FindNextFile(hFind, &fData) != 0);
#else
    DIR *dir_pointer = opendir(path.c_str());
    if (dir_pointer != NULL)
    {
        static int progress = 0;
        const char progress_chars[] = "-\\|/";
        struct dirent *entry_pointer;

        while ((entry_pointer = readdir(dir_pointer)))
        {
            const std::string fullname = path + "/" + entry_pointer->d_name;

            struct stat entry_status;
            if (!stat(fullname.c_str(), &entry_status))
            {
                if (S_ISDIR(entry_status.st_mode) && Recursive && (entry_pointer->d_name[0] != '.'))
                {
                    this->AppendFromDirectory(fullname + "/", wildcard);
                }
                else
                {
                    if (!fnmatch(wildcard, entry_pointer->d_name, FNM_PERIOD))
                    {
                        try
                        {
                            ImageFileDICOM::SmartPtr newFile(new ImageFileDICOM(fullname));

                            if (newFile->MatchAllPatterns(IncludePatterns) && !newFile->MatchAnyPattern(ExcludePatterns))
                            {
                                newFile->ReleaseDocument();
                                fileList.push_back(newFile);
                            }

                            if (Progress)
                            {
                                (cmtk::StdErr << "\r" << progress_chars[++progress % 4]).flush();
                            }
                        }
                        catch (...)
                        {
                            // not a valid DICOM file
                        }
                    }
                }
            }
        }
        (void)closedir(dir_pointer);
    }
#endif

    switch (SortFiles)
    {
    case SORT_NONE:
    default:
        break;
    case SORT_FILENAME:
        std::sort(fileList.begin(), fileList.end(), ImageFileDICOM::lessFileName);
        break;
    case SORT_INSTANCE_NUMBER:
        std::sort(fileList.begin(), fileList.end(), ImageFileDICOM::lessInstanceNumber);
        break;
    }

    for (std::vector<ImageFileDICOM::SmartConstPtr>::iterator it = fileList.begin(); it != fileList.end(); ++it)
    {
        try
        {
            this->AddImageFile(*it);
        }
        catch (int)
        {
        }
    }

    if (Progress)
    {
        cmtk::StdErr << "\r";
    }
}

int doMain(const int argc, const char *argv[])
{
    if (!dcmDataDict.isDictionaryLoaded())
    {
#ifdef CMTK_DCMDICTPATH
        if (dcmDataDict.wrlock().loadDictionary(CMTK_DCMDICTPATH))
        {
#if PACKAGE_VERSION_NUMBER < 364
            dcmDataDict.unlock();
#else
            dcmDataDict.wrunlock();
#endif
        }
        else
#endif
#ifdef CMTK_DCMDICTPATH_INSTALL
            if (dcmDataDict.wrlock().loadDictionary(CMTK_DCMDICTPATH_INSTALL))
        {
#if PACKAGE_VERSION_NUMBER < 364
            dcmDataDict.unlock();
#else
            dcmDataDict.wrunlock();
#endif
        }
        else
#endif
        {
            cmtk::StdErr << "Data dictionary not avaliable. Please set DCMDICTPATH variable as path to dicom.dic file.\n";
            throw cmtk::ExitException(1);
        }
    }

    try
    {
        cmtk::CommandLine cl;
        cl.SetProgramInfo(cmtk::CommandLine::PRG_TITLE, "DICOM to Image");
        cl.SetProgramInfo(cmtk::CommandLine::PRG_DESCR, "Combine sets of DICOM slices to 3D image stacks");

        typedef cmtk::CommandLine::Key Key;
        cl.AddSwitch(Key("no-progress"), &Progress, false, "Disable progress reporting.");

        cl.BeginGroup("Input", "Input Options");
        cl.AddSwitch(Key('r', "recurse"), &Recursive, true, "Recurse into directories");
        cl.EndGroup();

        cl.BeginGroup("Output", "Output Options");
        cl.AddOption(Key('O', "out-pattern"), &OutPathPattern, "Output image path pattern. Use the following substitutions: "
                                                               "printf-style %d variante (image number); "
                                                               "%n (image number with automatic number of digits); "
                                                               "%N (like %n, but with a hyphen '-' before number if there is more than one image); "
                                                               "%D (DICOM SeriesDescription); "
                                                               "%R (DICOM RepetitionTime - MRI only); "
                                                               "%E (DICOM EchoTime - MRI only); "
                                                               "%T (RawDataType - vendor-specific, currently GE MRI only)");

        cl.AddSwitch(Key('x', "xml"), &WriteXML, true, "Write XML sidecar file for each created image.");
        cl.AddSwitch(Key('d', "xml-only"), &WriteImage, false, "Do not write image file, only write the XML sidecar file.");
        cl.AddSwitch(Key("strict-xml"), &StrictXML, true, "Use strict XML when writing sidecar file.");
        cl.AddSwitch(Key("include-identifiers"), &IncludeIdentifiers, true, "Include potentially protected identifying information (e.g., UIDs, device serial numbers, dates) in the created XML sidecar files.");
        cl.AddSwitch(Key("include-ndar"), &IncludeNDARIdentifiers, true, "Include additional NDAR identifiers");

        cmtk::CommandLine::EnumGroup<cmtk::ImageStackDICOM::EmbedInfoEnum>::SmartPtr embedGroup = cl.AddEnum("embed", &EmbedInfo, "Embed DICOM information into output images as 'description' (if supported by output file format).");
        embedGroup->AddSwitch(Key("StudyID_StudyDate"), cmtk::ImageStackDICOM::EMBED_STUDYID_STUDYDATE, "StudyID, tag (0020,0010), then underscore, followed by StudyDate, tag (0008,0020). "
                                                                                                        "Date is appended because StudyID is four digits only and will repeat sooner or later.");
        embedGroup->AddSwitch(Key("PatientName"), cmtk::ImageStackDICOM::EMBED_PATIENTNAME, "Patient name, tag (0010,0010)");
        embedGroup->AddSwitch(Key("SeriesDescription"), cmtk::ImageStackDICOM::EMBED_SERIESDESCR, "Series description, tag (0008,103e)");
        embedGroup->AddSwitch(Key("None"), cmtk::ImageStackDICOM::EMBED_NONE, "Embed no information - leave 'description' field empty.");
        cl.EndGroup();

        cl.BeginGroup("Filtering", "Filtering Options")->SetProperties(cmtk::CommandLine::PROPS_ADVANCED);
        cl.AddCallback(Key("filter"), &CallbackAddInclude, "Filter DICOM files and include only those matching the given pattern of the form 'TagName=text', such that the value of the DICOM tag with the given name contains the given "
                                                           "text. If multiple filter patterns are provided via repeated use of this option, only files that match ALL patterns are included.");
        cl.AddCallback(Key("exclude"), &CallbackAddExclude, "Exclude all DICOM files matching the given pattern of the form 'TagName=text', such that the value of the DICOM tag with the given name contains the given text. "
                                                            "If multiple exclusion patterns are provided, all files are excluded that match ANY of the patterns.");
        cl.EndGroup();

        cl.BeginGroup("Sorting", "Sorting Options")->SetProperties(cmtk::CommandLine::PROPS_ADVANCED);
        cl.AddSwitch(Key("no-sort"), &SortFiles, SORT_NONE, "Do NOT sort files by file name (sorting determines image stack order when resolving spatial collisions)");
        cl.AddSwitch(Key("sort-by-name"), &SortFiles, SORT_FILENAME, "Sort files lexicographically by file name. Use this when instance numbers are non-unique.");
        cl.AddSwitch(Key("sort-by-instance"), &SortFiles, SORT_INSTANCE_NUMBER, "Sort files by image instance number. Use this when file names are different lengths, etc.");
        cl.EndGroup();

        cl.BeginGroup("Stacking", "Stacking Options")->SetProperties(cmtk::CommandLine::PROPS_ADVANCED);
        cl.AddSwitch(Key("write-single-slices"), &WriteSingleSlices, true, "Also write output images for single-slice DICOM files that could not be assigned to any 3D stacks. By default, these are skipped.");
        cl.AddSwitch(Key("ignore-acq-number"), &IgnoreAcquisitionNumber, true, "Ignore 'AcquisitionNumber' tag for image grouping, i.e., do not separate stacks based on this tag.");
        cl.AddSwitch(Key("no-orientation-check"), &DisableOrientationCheck, true, "Disable checking of image orientations (to avoid rounding issues)");
        cl.AddOption(Key("tolerance"), &Tolerance, "Tolerance for floating-point comparisons (must be >= 0; 0 = exact matches only; default: 1e-5). "
                                                   "If one or more volumes cannot be stacked because of non-uniform slice spacing that exceeds this threshold, dcm2image will return an exit code of 3.");
        cl.EndGroup();

        cl.AddParameterVector(&SearchRootDirVector, "SearchDirList", "List of directories to search for DICOM files. Subdirectories are also search if '--recurse' option is used.");

        cl.Parse(argc, const_cast<const char **>(argv));
    }
    catch (const cmtk::CommandLine::Exception &ex)
    {
        cmtk::StdErr << ex << "\n";
        throw cmtk::ExitException(1);
    }

    if (!dcmDataDict.rdlock().findEntry("RawDataType_ImageType"))
    {
#if PACKAGE_VERSION_NUMBER < 364
        dcmDataDict.unlock();
#else
        dcmDataDict.rdunlock();
#endif
        dcmDataDict.wrlock().addEntry(new DcmDictEntry(0x0043, 0x102f, EVR_SS, "RawDataType_ImageType", 1, 1, NULL, OFFalse, "GE"));
#if PACKAGE_VERSION_NUMBER < 364
        dcmDataDict.unlock();
#else
        dcmDataDict.wrunlock();
#endif
    }

    VolumeList volumeList(Tolerance);
    for (std::vector<std::string>::const_iterator it = SearchRootDirVector.begin(); it != SearchRootDirVector.end(); ++it)
    {
        volumeList.AppendFromDirectory(*it, "*");
    }

    return volumeList.WriteVolumes();
}

#include "cmtkSafeMain"
