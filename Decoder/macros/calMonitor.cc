#include "NAPStyle.h"
#include "NIConfig.hh"
#include "NIConfig.cc"
#include "NAUtil.h"
#include "NAUtil.cc"

void calMonitor( const std::string& inputFileList, const std::string& configFile, const std::string& outputDir )
{
    std::list< std::string > filePathList;
    if( NAUtil::GetLines( inputFileList, &filePathList ) == false ) return;

    NIConfig* pConfig = new NIConfig( );
    if( pConfig->ReadConfigJSON( configFile ) == false ) {
        NAUtil::Cerr( "failed to read CONFIG file..." );
        return;
    }

    pConfig->PrintConfigJSON( );

    NAUtil::ExistCreateDir( outputDir );
    // TFile outFile( Form( "%s/out.root", outputDir.c_str( ) ), "RECREATE" );
    SetNAPStyle( );
    gStyle->SetMarkerStyle(7);
    // gStyle->SetMarkerSize( 0.04 );

    // debug
    // for( auto path : filePathList ) std::cout << path << std::endl;

    // analysis
    double calFactor = pConfig->cal_factor;

    int ev = 0;
    int fileID = 0;
    vector< double >* pVec_xz_x = nullptr;
    vector< double >* pVec_xz_z = nullptr;
    vector< double >* pVec_yz_y = nullptr;
    vector< double >* pVec_yz_z = nullptr;
    double ave_x = 0.0, ave_y = 0.0, ave_z = 0.0;
    double a_hg_sum_pulse_height = 0.0, a_hg_sum_charge = 0.0;

    double refEnergy = -100.0;
    TF1 gausFunc( "gausFunc", "gaus", 80.0, 180.0 );
    
    double perList[100]    = { };
    double perListErr[100] = { };
    double eneList[100]    = { };
    double eneListErr[100] = { };
    int nFiles = filePathList.size( );
    int idx = 0;
    for( auto filePath : filePathList ) {
        std::string fileName = NAUtil::GetFileName( filePath );
        std::string perIDStr = fileName.substr( 3, fileName.find( "_" ) - 3);
        int perID = std::stoi( perIDStr );
        std::cout << perID << std::endl;

        TFile file( filePath.c_str( ) );
        TTree* pTree = dynamic_cast< TTree* >( file.Get( "ana_tree" ) );
        if( pTree == nullptr ) continue;
        
        pTree->SetBranchAddress( "eventID",         &ev              );
        pTree->SetBranchAddress( "fileID",          &fileID          );
        pTree->SetBranchAddress( "xz_x",            &pVec_xz_x       );
        pTree->SetBranchAddress( "xz_z",            &pVec_xz_z       );
        pTree->SetBranchAddress( "yz_y",            &pVec_yz_y       );
        pTree->SetBranchAddress( "yz_z",            &pVec_yz_z       );
        pTree->SetBranchAddress( "a_hg_sum_charge", &a_hg_sum_charge );

        std::string baseFileName = NAUtil::GetBaseFileName( fileName );
        std::string histName = Form( "hist_energy_%s", baseFileName.c_str( ) );
        TH1D hist( histName.c_str( ), histName.c_str( ), 100, 0, 400 );
        
        int totEvt = pTree->GetEntries( );
        for( int evtIdx = 0; evtIdx < totEvt; ++evtIdx ) {
            pTree->GetEntry( evtIdx );
            
            double xz_x_min = 1000.0, xz_x_max = -1000.0;
            double xz_z_min = 1000.0, xz_z_max = -1000.0;
            for( int i = 0; i < pVec_xz_x->size( ); ++i ) {
                if( xz_x_min > pVec_xz_x->at( i ) ) {
                    xz_x_min = pVec_xz_x->at( i );
                    xz_z_min = pVec_xz_z->at( i );
                }
                if( xz_x_max < pVec_xz_x->at( i ) ) {
                    xz_x_max = pVec_xz_x->at( i );
                    xz_z_max = pVec_xz_z->at( i );
                }
            }
            double yz_y_min = 1000.0, yz_y_max = -1000.0;
            double yz_z_min = 1000.0, yz_z_max = -1000.0;
            for( int i = 0; i < pVec_yz_y->size( ); ++i ) {
                if( yz_y_min > pVec_yz_y->at( i ) ) {
                    yz_y_min = pVec_yz_y->at( i );
                    yz_z_min = pVec_yz_z->at( i );
                }
                if( yz_y_max < pVec_yz_y->at( i ) ) {
                    yz_y_max = pVec_yz_y->at( i );
                    yz_z_max = pVec_yz_z->at( i );
                }
            }

            double lengthXZ = sqrt( pow( xz_x_max - xz_x_min, 2) + pow( yz_y_max - yz_y_min, 2) + pow( xz_z_max - xz_z_min, 2) );
            if( lengthXZ > 1.0 && lengthXZ < 1.3 ) hist.Fill( a_hg_sum_charge * calFactor );
        }
        
        hist.Fit( &gausFunc, "LMI", "", 50.0, 150.0 );

        // hist.SetDirectory( &outFile );
        // hist.Write( );
        // hist.SetDirectory( nullptr );

        double peakVal    = gausFunc.GetParameter( 1 );
        double peakValErr = gausFunc.GetParError( 1 );
        if( refEnergy < 0.0 ) refEnergy = peakVal;
        
        perList[idx]    = perID;        
        perListErr[idx] = 0.5;
        eneList[idx]    = peakVal / refEnergy;        
        eneListErr[idx] = peakValErr / refEnergy; // tmp
        ++idx;
        file.Close( );
    }
    
    TGraphErrors graph( nFiles, perList, eneList, perListErr, eneListErr );
    TCanvas cvs( "cvs", "cvs", 800, 600 );
    graph.SetMarkerSize( 4 );
    graph.Draw( "AP" );
    graph.GetXaxis( )->SetTitle( "period" );
    graph.GetYaxis( )->SetTitle( "Energy scale (normalized by the first period) " );

    cvs.SaveAs( Form( "%s/ene.png", outputDir.c_str( ) ) );
    
    return;
}
