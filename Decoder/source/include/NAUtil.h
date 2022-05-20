//////////////////////////////////////////////////////////////////
//
// Utility
//
// Satoshi Higashino
// satoshi.higashino@cern.ch / higashino@people.kobe-u.ac.jp
//
//////////////////////////////////////////////////////////////////
#ifndef NA_UTIL_H
#define NA_UTIL_H

// sys
#include <sys/types.h>
#include <sys/stat.h>

// C++ basics
#include <iostream>
#include <cstring>
#include <sstream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdio>

// STL
#include <set>
#include <vector>
#include <map>
#include <list>

// File operation
#include <dirent.h>

// io manipulator
#include <iomanip>

// ---------------------------------------------------------------
// Constant parameter definition
// ---------------------------------------------------------------

// // for debug
// #define DEBUG(val) std::cout<<"Debugging : "<<val<<std::endl;

class NAUtil
{
public:

    //////////////////////////////////////////////////////////////////
    //
    // Converter
    //
    //////////////////////////////////////////////////////////////////

    // ---------------------------------------------------------------
    // -> string data converter
    template <typename T> static std::string ToString( const T& val );

    // ---------------------------------------------------------------
    // -> string data converter
    static int StrToInt( const std::string& str );

    // ---------------------------------------------------------------
    // -> string data converter
    static double StrToDouble( const std::string& str );

    // ---------------------------------------------------------------
    // Convert angle unit
    // ---------------------------------------------------------------
    // ---------------------------------------------------------------
    // Radian -> Degree
    static double RadToDeg( double rad );

    // ---------------------------------------------------------------
    // Degree -> Radian
    static double DegToRad( double deg );

    // ---------------------------------------------------------------
    // Convert timee unit
    // ---------------------------------------------------------------
    // ---------------------------------------------------------------
    // hour -> minute
    static double HourToMin( double hour );

    // ---------------------------------------------------------------
    // hour -> second
    static double HourToSec( double hour );

    // ---------------------------------------------------------------
    // minute -> hour
    static double MinToHour( double min );

    // ---------------------------------------------------------------
    // minute -> second
    static double MinToSec( double min );

    // ---------------------------------------------------------------
    // second -> hour
    static double SecToHour( double sec );

    // ---------------------------------------------------------------
    // second -> minute
    static double SecToMin( double sec );

    // ---------------------------------------------------------------
    // Console Output
    // ---------------------------------------------------------------
    // ---------------------------------------------------------------
    // only applied for string value.
    // eg. MuUitl::Cout( "Example" );
    // template <typename T> static void Cout( const T& val );
    static void Cout( const std::string& val );

    // ---------------------------------------------------------------
    // output error message
    // template <typename T> static void Cerr( const T& val );
    static void Cerr( const std::string& val );

    // ---------------------------------------------------------------
    // output error message
    // template <typename T> static void Cwarn( const T& val );
    static void Cwarn( const std::string& val );

    // ---------------------------------------------------------------
    // output error message
    // template <typename T> static void Cinfo( const T& val );
    static void Cinfo( const std::string& val );

    //////////////////////////////////////////////////////////////////
    //
    // Comparetor
    //
    //////////////////////////////////////////////////////////////////
    // ---------------------------------------------------------------
    // Comparing "double" value
    // at 10e-6 (the significance of float value)
    // return :
    //     1 for first > second
    //    -1 for first < second
    //     0 for first = second
    static int CompareDouble( const double& first, const double& second );

    // ---------------------------------------------------------------
    // Comparing "double" value to Zero
    // at 10e-6 (the significance of float value)
    static bool IsZero( const double& val );

    //////////////////////////////////////////////////////////////////
    //
    // File Checker
    //
    //////////////////////////////////////////////////////////////////
    // ---------------------------------------------------------------
    // return:
    //    true ...exist 
    //    false...NOT exist 
    static bool ExistFile( const std::string& filePath );

    // ---------------------------------------------------------------
    // Check directry
    // return:
    //    true ...exist 
    //    false...NOT exist 
    // comment:
    //    faster to use S_ISDIR( ) instead of opendir( )?
    static bool ExistDir( const std::string& dirPath );

    // ---------------------------------------------------------------
    // Check directry: if the directory is not found, create it
    static void ExistCreateDir( const std::string& dirPath );

    // ---------------------------------------------------------------
    // Check directry belonging to file path
    // return:
    //    true ...exist 
    //    false...NOT exist 
    static bool ExistFilePathDir( const std::string& filePath );

    // ---------------------------------------------------------------
    // Get file name from full path
    // return: file name
    static std::string GetFileName( const std::string& path );

    // ---------------------------------------------------------------
    // Extract extension from file name
    // return: base file name (without extension)
    static std::string GetBaseFileName( const std::string& path );

    // ---------------------------------------------------------------
    // get file path list/array from input directory
    // argument:
    //     [input]  inputDir           ... input directory (Don't include "/" at the end of directory name!!!)
    //     [output] pFilePathList(Arr) ... relative file path from input directory
    //     [input]  recursive          ... if true, output file path recursively
    //     [input]  fileKey            ... picking up only containing "fileKey" phrase in a file name
    //     [input]  dirKey             ... picking up only containing "dirKey" phrase in a directory name
    // return: 
    //     true ... success
    //     false... failure
    static bool GetFilePathList( const std::string&        inputDir,
                                 std::list< std::string >* pFilePathList,
                                 const bool&               recursive    = true,
                                 const std::string&        fileKey      = "",
                                 const std::string&        dirKey       = "" );

    static bool GetFilePathArr( const std::string&          inputDir,
                                std::vector< std::string >* pFilePathArr,
                                const bool&                 recursive    = true,
                                const std::string&          fileKey      = "",
                                const std::string&          dirKey       = "" );
    
    // ---------------------------------------------------------------
    // get lines from a text file
    // argument:
    //     [input]  inputFile  ... input filename
    //     [output] pList      ... line list
    //     [input]  commSyntax ... comment syntax (default: #)
    // return: 
    //     true ... success
    //     false... failure
    static bool GetLines( const std::string&        inputDir,
                          std::list< std::string >* pList,
                          const std::string&        commSyntax = "#" );

    //////////////////////////////////////////////////////////////////
    //
    // Calculation
    //
    //////////////////////////////////////////////////////////////////
    // ---------------------------------------------------------------
    // Get Sum of Squares
    // 
    static double GetSumOfSqrt( const double& val1, const double& val2 );
    static double GetSumOfSqrt( const double& val1, const double& val2, const double& val3 );

    //////////////////////////////////////////////////////////////////
    //
    // std::string Operation
    //
    //////////////////////////////////////////////////////////////////
    // ---------------------------------------------------------------
    // Replacement
    // 
    static std::string Replace( const std::string& inputStr, const std::string& targetStr, const std::string& replacedStr );

    // ---------------------------------------------------------------
    // Split (by white spaces)
    // 
    static bool Split( const std::string& inputStr, const char delim, std::vector< std::string >* pArr );

    //////////////////////////////////////////////////////////////////
    //
    // Progress Bar
    //
    //////////////////////////////////////////////////////////////////
    static void PrintProgressBar( const int& index, const int& total );
    
};

#endif // NA_UTIL_H
