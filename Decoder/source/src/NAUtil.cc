//////////////////////////////////////////////////////////////////
//
// Utility
//
// Satoshi Higashino
// satoshi.higashino@cern.ch / higashino@people.kobe-u.ac.jp
//
//////////////////////////////////////////////////////////////////
#include "NAUtil.h"

//////////////////////////////////////////////////////////////////
//
// Converter
//
//////////////////////////////////////////////////////////////////

// ---------------------------------------------------------------
// -> string data converter
template <typename T> std::string NAUtil::ToString( const T& val )
{
    std::string str = "";
    std::stringstream stream;

    stream << val;
    stream >> str;

    return str;
}

// ---------------------------------------------------------------
// -> string data converter
int NAUtil::StrToInt( const std::string& str )
{
    std::stringstream stream( str );
    int retVal = 0;
        
    // stream << str;
    stream >> retVal;

    return retVal;
}

// ---------------------------------------------------------------
// -> string data converter
double NAUtil::StrToDouble( const std::string& str )
{
    std::stringstream stream;
    double retVal = 0.0;
        
    stream << str;
    stream >> retVal;

    return retVal;
}

// ---------------------------------------------------------------
// Convert angle unit
// ---------------------------------------------------------------
// ---------------------------------------------------------------
// Radian -> Degree
double NAUtil::RadToDeg( double rad )
{
    return rad * 180.0 / 3.141593;
}

// ---------------------------------------------------------------
// Degree -> Radian
double NAUtil::DegToRad( double deg )
{
    return deg * 3.141593 / 180.0;
}

// ---------------------------------------------------------------
// Convert timee unit
// ---------------------------------------------------------------
// ---------------------------------------------------------------
// hour -> minute
double NAUtil::HourToMin( double hour )
{
    return hour * 60.0;
}

// ---------------------------------------------------------------
// hour -> second
double NAUtil::HourToSec( double hour )
{
    return hour * 3600.0;
}

// ---------------------------------------------------------------
// minute -> hour
double NAUtil::MinToHour( double min )
{
    return min / 60.0;
}

// ---------------------------------------------------------------
// minute -> second
double NAUtil::MinToSec( double min )
{
    return min * 60.0;
}

// ---------------------------------------------------------------
// second -> hour
double NAUtil::SecToHour( double sec )
{
    return sec * 3600.0;
}

// ---------------------------------------------------------------
// second -> minute
double NAUtil::SecToMin( double sec )
{
    return sec * 60.0;
}

// ---------------------------------------------------------------
// Console Output
// ---------------------------------------------------------------
// ---------------------------------------------------------------
// only applied for string value.
// eg. NAUtil::Cout( "Example" );
void NAUtil::Cout( const std::string& val )
{
    std::cout << val << std::endl;
}

// // ---------------------------------------------------------------
// // output error message
void NAUtil::Cerr( const std::string& val )
{
    std::cout << "NAUtil Error: " << val << std::endl;
}

// ---------------------------------------------------------------
// output error message
void NAUtil::Cwarn( const std::string& val )
{
    std::cout << "NAUtil Warning: " << val << std::endl;
}

// ---------------------------------------------------------------
// output error message
void NAUtil::Cinfo( const std::string& val )
{
    std::cout << "NAUtil Information: " << val << std::endl;
}

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
int NAUtil::CompareDouble( const double& first, const double& second )
{
    if( fabs( first - second ) < 0.000001 ) {
        return 0;
    }
    else {
        if( first > second ) return 1;
        else                 return -1;
    }
}

// ---------------------------------------------------------------
// Comparing "double" value to Zero
// at 10e-6 (the significance of float value)
bool NAUtil::IsZero( const double& val )
{
    if( CompareDouble( val, 0.0 ) == 0 ) return true;
    return false;
}

//////////////////////////////////////////////////////////////////
//
// File Checker
//
//////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------
// return:
//    true ...exist 
//    false...NOT exist 
bool NAUtil::ExistFile( const std::string& filePath )
{
    bool exist = false;
    FILE* fp = fopen( filePath.c_str( ), "r" );
    if( fp != nullptr ){
        exist = true;
        fclose( fp );
    }
    return exist;
}

// ---------------------------------------------------------------
// Check directry
// return:
//    true ...exist 
//    false...NOT exist 
// comment:
//    faster to use S_ISDIR( ) instead of opendir( )?
bool NAUtil::ExistDir( const std::string& dirPath )
{
    bool exist = false;
    DIR* dp = opendir( dirPath.c_str( ) );
    if( dp != nullptr ) {
        exist = true;
        closedir( dp );
    }
    return exist;
}

// ---------------------------------------------------------------
// Check directry: if the directory is not found, create it
void NAUtil::ExistCreateDir( const std::string& dirPath )
{
    if( ExistDir( dirPath ) == true ) return;

    // create directory
    if( mkdir( dirPath.c_str( ), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH ) == 0 ) {
        std::string output = "create directory : " + dirPath;
        Cinfo( output );
    }
    return;
}

// ---------------------------------------------------------------
// Check directry belonging to file path
// return:
//    true ...exist 
//    false...NOT exist 
bool NAUtil::ExistFilePathDir( const std::string& filePath )
{
    int index = filePath.rfind( "/", filePath.size( ) - 1 );
    std::string dirPath = filePath.substr( 0, index );

    return NAUtil::ExistDir( dirPath );
}

// ---------------------------------------------------------------
// Get file name from full path
// return: file name
std::string NAUtil::GetFileName( const std::string& path )
{
    size_t pos1 = path.rfind( '/' );
    if( pos1 != std::string::npos ) {
        return path.substr( pos1 + 1, path.size( ) - pos1 - 1 );
    }
        
    return path;
}

// ---------------------------------------------------------------
// Extract extension from file name
// return: file name (without extension)
std::string NAUtil::ExtractPathWithoutExt( const std::string& path )
{
    std::string::size_type pos;
    if( ( pos = path.find_last_of( "." ) ) == std::string::npos ) {
        return path;
    }

    return path.substr( 0, pos );
}

// ---------------------------------------------------------------
// get file path list from input directory
// argument:
//     [input]  inputDir      ... input directory (Don't include "/" at the end of directory name!!!)
//     [output] pFilePathList ... relative file path from input directory
//     [input]  recursive     ... if true, output file path recursively
//     [input]  fileKey       ... picking up only containing "fileKey" phrase in a file name
//     [input]  dirKey        ... picking up only containing "dirKey" phrase in a directory name
// return: 
//     true ... success
//     false... failure
bool NAUtil::GetFilePathList( const std::string&        inputDir,
                              std::list< std::string >* pFilePathList,
                              const bool&               recursive,
                              const std::string&        fileKey,
                              const std::string&        dirKey )
{
    if( pFilePathList == nullptr ) return false;
    // if( ExistDir( inputDir ) == false ) return false;

    bool retVal = true;
    DIR* dp = opendir( inputDir.c_str( ) );
    if( dp != nullptr ) {
        struct dirent* entry;
        struct stat st;
        std::string name = "";
            
        while( (entry = readdir( dp ) ) != nullptr ) {
            name = entry->d_name;
            std::string filePath = inputDir + "/";
            filePath += name;
                
            // skip "." and ".."
            if( name.length( ) <= 0 || name == "." || name == ".." ) continue;

            // skip stat file
            if( stat( filePath.c_str( ), &st ) == -1 ) continue;

            if( S_ISDIR( st.st_mode ) == true ) {
                if( recursive == true ) {
                    if( GetFilePathList( filePath, pFilePathList, recursive, fileKey, dirKey ) == false )
                        retVal = false;
                }
            }
            else {
                if( fileKey.length( ) > 0 && name.find( fileKey ) == std::string::npos ) continue;
                if( dirKey.length( ) > 0 && filePath.find( dirKey ) == std::string::npos ) continue;
                pFilePathList->push_back( filePath );
                // std::cout << "Add file: " << filePath << std::endl;
            }
        }

        closedir( dp );
    }
    else {
        Cerr( "error!!!" );
    }

    return retVal;
}


bool NAUtil::GetFilePathArr( const std::string&          inputDir,
                             std::vector< std::string >* pFilePathArr,
                             const bool&                 recursive,
                             const std::string&          fileKey,
                             const std::string&          dirKey )
{
    if( pFilePathArr == nullptr ) return false;

    bool retVal = true;
    DIR* dp = opendir( inputDir.c_str( ) );
    if( dp != nullptr ) {
        struct dirent* entry;
        struct stat st;
        std::string name = "";
            
        while( (entry = readdir( dp ) ) != nullptr ) {
            name = entry->d_name;
            std::string filePath = inputDir + "/";
            filePath += name;
                
            // skip "." and ".."
            if( name.length( ) <= 0 || name == "." || name == ".." ) continue;

            // skip stat file
            if( stat( filePath.c_str( ), &st ) == -1 ) continue;

            if( S_ISDIR( st.st_mode ) == true ) {
                if( recursive == true ) {
                    if( GetFilePathArr( filePath, pFilePathArr, recursive, fileKey, dirKey ) == false )
                        retVal = false;
                }
            }
            else {
                if( fileKey.length( ) > 0 && name.find( fileKey ) == std::string::npos ) continue;
                if( dirKey.length( ) > 0 && filePath.find( dirKey ) == std::string::npos ) continue;
                pFilePathArr->push_back( filePath );
                // std::cout << "Add file: " << filePath << std::endl;
            }
        }

        closedir( dp );
    }
    else {
        Cerr( "error!!!" );
    }

    return retVal;
}


// ---------------------------------------------------------------
// get lines from a text file
// argument:
//     [input]  inputFile  ... input filename
//     [output] pList      ... line list
//     [input]  commSyntax ... comment syntax (default: #)
// return: 
//     true ... success
//     false... failure
bool NAUtil::GetLines( const std::string&        inputFile,
                       std::list< std::string >* pList,
                       const std::string&        commSyntax )
{
    if( pList == nullptr ) return false;

    std::ifstream ifs( inputFile );
    if( ifs.is_open( ) == false ) return false;

    while( !ifs.eof( ) ) {
        std::string line = "";
        std::getline( ifs, line );
        if( line.length( ) <= 0 || strncmp( line.c_str( ), commSyntax.c_str( ), 1 ) == 0 ) continue;
            
        pList->push_back( line );
    }

    if( pList->size( ) <= 0 ) Cwarn( "NAUtil::getLines( ) ... no lines are retrieved." );
        
    return true;
}



//////////////////////////////////////////////////////////////////
//
// Calculation
//
//////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------
// Get Sum of Squares
// 
double NAUtil::GetSumOfSqrt( const double& val1, const double& val2 )
{
    return sqrt( val1 * val1 + val2 * val2 );
}
double NAUtil::GetSumOfSqrt( const double& val1, const double& val2, const double& val3 )
{
    return sqrt( val1*val1 + val2*val2 + val3*val3 );
}

//////////////////////////////////////////////////////////////////
//
// String Operation
//
//////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------
// Replacement
// 
std::string NAUtil::Replace( const std::string& inputStr, const std::string& targetStr, const std::string& replacedStr )
{
    std::string retVal = "";
    size_t pos = inputStr.find( targetStr );
    if( pos != std::string::npos ) {
        retVal = inputStr;
        retVal.replace( pos, targetStr.length( ), replacedStr.c_str( ) );
    }
    else {
        retVal = inputStr;
    }

    return retVal;
}

// ---------------------------------------------------------------
// Split (by white spaces)
// 
bool NAUtil::Split( const std::string& inputStr, const char delim, std::vector< std::string >* pArr )
{
    if( pArr == nullptr ) return false;
    std::stringstream ss( inputStr );
    std::string element = "";
    
    while( std::getline( ss, element, delim ) ) pArr->push_back( element );
        
    return true;
}


//////////////////////////////////////////////////////////////////
//
// Progress Bar
//
//////////////////////////////////////////////////////////////////
void NAUtil::PrintProgressBar( const int& index, const int& total )
{
    if( index % 10 == 0 ) {
        std::string printBar = " [";
        double progress = static_cast< double >( index ) / static_cast< double >( total );
        for( int bar = 0; bar < 20; ++bar ) {
            double currentFraction = static_cast< double >( bar ) * 0.05;
            if( progress > currentFraction ) printBar += "/";
            else printBar += ".";
        }
        printBar += "] ";
        double percent = 100.0 * progress;
        std::stringstream percentSS;
        percentSS << std::setprecision( 2 ) << percent;
        std::string text = printBar + " ";
        text += percentSS.str( );
        std::cout << std::flush; 
        std::cout << text << "%\r" << std::flush; 
    }
    return;
}
    
