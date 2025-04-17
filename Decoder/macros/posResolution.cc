#include "NAPStyle.h"
#include "NIConfig.hh"
#include "NIConfig.cc"
#include "NAUtil.h"
#include "NAUtil.cc"

#define SAMPLING_HELZ 2.5e6 //Hz

bool getFileTable( const std::string& inputFiles, std::map< std::string, double >* pTable )
{
    if( pTable == nullptr ) return false;

    std::ifstream ifs( inputFiles );
    if( ifs.is_open( ) == false ) return false;
    while( ifs.eof( ) == false ) {
        std::string line = "";
        std::getline( ifs, line );
        if( line.length( ) <= 0 || strncmp( line.c_str( ), "#", 1 ) == 0 ) continue;

        std::string filename = "";
        double value = 0.0;
        std::stringstream ss( line );
        ss >> filename >> value;
        pTable->insert( std::make_pair( filename, value ) );
    }

    return true;
}

TLatex* CreateDrawText( const double&       x,
                        const double&       y,
                        const std::string&  text,
                        const double&       size = 0.05,
                        const Color_t&      color = 1 )
{
    if( text.size( ) <= 0 ) return nullptr;
    TLatex l;
    l.SetNDC( );
    l.SetTextColor( color );
    l.SetTextSize( size );
    return l.DrawLatex( x, y, text.c_str( ) );
}

void posResolution( const std::string& inputFiles, const std::string& configFile, const std::string& outputDir )
{
    SetNAPStyle( );

    std::map< std::string, double > inputFileTable;
    if( getFileTable( inputFiles, &inputFileTable ) == false ) return;
    // if( NAUtil::GetFilePathArr( inputFiles, &inputFileList ) == false ) return;

    NIConfig* pConfig = new NIConfig( );
    if( pConfig->ReadConfigJSON( configFile ) == false ) {
        NAUtil::Cerr( "failed to read CONFIG file..." );
        return;
    }

    pConfig->PrintConfigJSON( );

    NAUtil::ExistCreateDir( outputDir );
    
    bool   isAlphaCalib = pConfig->is_alpha_calib == 0 ? false : true;
    double calFactor    = pConfig->cal_factor;
    double driftV_sf6   = pConfig->driftV_main;
    double driftV_sf5   = pConfig->driftV_mino;
    int    wf_int_range_min = pConfig->wf_integral_range_min;
    int    wf_int_range_max = pConfig->wf_integral_range_max;
    double minority_ROI_range = pConfig->minority_ROI_range;
    double minority_ROI_offset = pConfig->minority_ROI_offset;

    int ev = 0, fileID = 0;
    std::vector< double >* pVec_xz_x = nullptr;
    std::vector< double >* pVec_xz_z = nullptr;
    std::vector< double >* pVec_yz_y = nullptr;
    std::vector< double >* pVec_yz_z = nullptr;
    double ave_x = 0.0, ave_y = 0.0, ave_z = 0.0;
    double a_hg_sum_pulse_height = 0.0, a_hg_sum_charge = 0.0;
    double c_hg_sum_pulse_height = 0.0, c_hg_sum_charge = 0.0;
    std::vector< double >* pVec_a_hg_mainrise_time = nullptr;
    std::vector< double >* pVec_a_hg_mainfall_time = nullptr;
    std::vector< double >* pVec_c_hg_mainrise_time = nullptr;
    std::vector< double >* pVec_c_hg_mainfall_time = nullptr;
    std::vector< double >* pVec_dt = nullptr;

    std::vector< double >* pVec_a_hg_mino_peak_wf = nullptr;
    int wfBin = wf_int_range_max - wf_int_range_min;
    double wf_time_min = (double)wf_int_range_min/SAMPLING_HELZ*1e6;
    double wf_time_max = (double)wf_int_range_max/SAMPLING_HELZ*1e6;

    std::map< TH1F*, double > histTable;

    for( auto pair : inputFileTable ) {
        std::string filename = pair.first;
        double sourcePos = pair.second;
        if( sourcePos < 10.0 ) continue;

        TFile file( filename.c_str( ) );
        TTree* pTree = dynamic_cast< TTree* >( file.Get( "ana_tree" ) );
        if( pTree == nullptr ) {
            NAUtil::Cerr( "innput file doesn't include TTree." );
            NAUtil::Cerr( filename );
            continue;
        }

        std::cout << "Processing " << filename << std::endl;

        std::string histZPosName = Form( "hist_sourcePos_%d", static_cast< int >( sourcePos ) );
        TH1F* pHistZPos = new TH1F( histZPosName.c_str( ), histZPosName.c_str( ), 32, 0.0, 160.0 );
        pHistZPos->SetDirectory( nullptr );
        histTable.insert( std::make_pair( pHistZPos, sourcePos ) );

        std::string histZPosAllStripName = Form( "hist_sourcePosAllStrip_%d", static_cast< int >( sourcePos ) );
        TH1F* pHistZPosAllStrip = new TH1F( histZPosAllStripName.c_str( ), histZPosAllStripName.c_str( ), 32, 0.0, 160.0 );
        pHistZPosAllStrip->SetDirectory( nullptr );
        histTable.insert( std::make_pair( pHistZPosAllStrip, sourcePos ) );

        std::string histZPosGausName = Form( "hist_sourcePosGaus_%d", static_cast< int >( sourcePos ) );
        TH1F* pHistZPosGaus = new TH1F( histZPosGausName.c_str( ), histZPosGausName.c_str( ), 32, 0.0, 160.0 );
        pHistZPosGaus->SetDirectory( nullptr );
        histTable.insert( std::make_pair( pHistZPosGaus, sourcePos ) );

        std::string histZPosSumName = Form( "hist_sourcePosSum_%d", static_cast< int >( sourcePos ) );
        // TH1F* pHistZPosSum = new TH1F( histZPosSumName.c_str( ), histZPosSumName.c_str( ), 32, 0.0, 160.0 );
        TH1F* pHistZPosSum = new TH1F( histZPosSumName.c_str( ), histZPosSumName.c_str( ), 40, 0.0, 160.0 );
        pHistZPosSum->SetDirectory( nullptr );
        histTable.insert( std::make_pair( pHistZPosSum, sourcePos ) );

        pTree->SetBranchAddress( "eventID", &ev );
        pTree->SetBranchAddress( "fileID", &fileID );
        pTree->SetBranchAddress( "xz_x", &pVec_xz_x );
        pTree->SetBranchAddress( "xz_z", &pVec_xz_z );
        pTree->SetBranchAddress( "yz_y", &pVec_yz_y );
        pTree->SetBranchAddress( "yz_z", &pVec_yz_z );
        pTree->SetBranchAddress( "ave_x", &ave_x );
        pTree->SetBranchAddress( "ave_y", &ave_y );
        pTree->SetBranchAddress( "ave_z", &ave_z );
        pTree->SetBranchAddress( "a_hg_sum_pulse_height", &a_hg_sum_pulse_height );
        pTree->SetBranchAddress( "a_hg_sum_charge",       &a_hg_sum_charge       );
        pTree->SetBranchAddress( "c_hg_sum_pulse_height", &c_hg_sum_pulse_height );
        pTree->SetBranchAddress( "c_hg_sum_charge",       &c_hg_sum_charge       );
        pTree->SetBranchAddress( "a_hg_mainrise_time",    &pVec_a_hg_mainrise_time );
        pTree->SetBranchAddress( "a_hg_mainfall_time",    &pVec_a_hg_mainfall_time );
        pTree->SetBranchAddress( "c_hg_mainrise_time",    &pVec_c_hg_mainrise_time );
        pTree->SetBranchAddress( "c_hg_mainfall_time",    &pVec_c_hg_mainfall_time );
        pTree->SetBranchAddress( "dt",                    &pVec_dt );
        pTree->SetBranchAddress( "a_hg_mino_peak_wf",     &pVec_a_hg_mino_peak_wf );

        int totEvt = pTree->GetEntries( );
        for( int evt = 0; evt < totEvt; ++evt ) {

            pVec_xz_x = nullptr;
            pVec_xz_z = nullptr;
            pVec_yz_y = nullptr;
            pVec_yz_z = nullptr;
            pVec_a_hg_mainrise_time = nullptr;
            pVec_a_hg_mainfall_time = nullptr;
            pVec_c_hg_mainrise_time = nullptr;
            pVec_c_hg_mainfall_time = nullptr;
            pVec_dt = nullptr;
            pVec_a_hg_mino_peak_wf = nullptr;

            pTree->GetEntry( evt );
            if( pVec_yz_y->size( ) <= 0 ) continue;


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
            double lengthYZ = sqrt( pow( xz_x_max - xz_x_min, 2) + pow( yz_y_max - yz_y_min, 2) + pow( yz_z_max - yz_z_min, 2) );
            double lengthXY = sqrt( pow( xz_x_max - xz_x_min, 2) + pow( yz_y_max - yz_y_min, 2) );
            double lengthX  = sqrt( pow( xz_x_max - xz_x_min, 2) );
            double lengthY  = sqrt( pow( yz_y_max - yz_y_min, 2) );

            // fiducial cut
            if( xz_x_max < 1.1 || xz_x_min > 0.0 ) continue;

            // energy cut
            if( a_hg_sum_charge * calFactor < 10.0 ) continue;


            double sumToT = 0.0;
            int nHit = pVec_a_hg_mainrise_time->size( );
            double mainrise_sum = 0.0, mainrise_sqsum = 0.0;
            if( pVec_a_hg_mainfall_time->size( ) == nHit ) {
                for( int i = 0; i < nHit; ++i ) {
                    sumToT += ( pVec_a_hg_mainfall_time->at( i ) - pVec_a_hg_mainrise_time->at( i ) );
                    mainrise_sum += pVec_a_hg_mainrise_time->at( i );
                    mainrise_sqsum += pow(pVec_a_hg_mainrise_time->at( i ), 2.0);
                }
            }
            double mainrise_mean = mainrise_sum / (double)nHit;
            double mainrise_sigma = sqrt(mainrise_sqsum / (double)nHit - pow(mainrise_mean, 2));
            if( lengthXZ < 1.3 ) {
                // at least 2 hits are required in cathode strips
                if( pVec_yz_y->size( ) <= 1 ) continue;
                if( mainrise_sigma > 4.0 ) continue;
                // if( mainrise_sigma > 8.0 ) continue;

                double averageDt = 0.0;
                if( pVec_dt->size( ) > 0 ) averageDt = std::accumulate( pVec_dt->begin( ), pVec_dt->end( ), 0.0 ) / pVec_dt->size( );

                double absZ = averageDt / ( (1.0 / driftV_sf6) - (1.0 / driftV_sf5) );
                pHistZPos->Fill( absZ * 10.0 ); // cm -> mm

                TH1D histZPosAllStripOneEvt( "histZPosAllStripOneEvt", "histZPosAllStripOneEvt", 160, 0.0, 160.0 );
                for( auto eachDt : *pVec_dt ) {
                    double absZEachDt = eachDt / ( (1.0 / driftV_sf6) - (1.0 / driftV_sf5) );
                    pHistZPosAllStrip->Fill( absZEachDt * 10.0 );
                    histZPosAllStripOneEvt.Fill( absZEachDt * 10.0 );
                }
                
                TF1 fitGaus( "fitGaus", "gaus", 0.0, 160.0 );
                fitGaus.SetParameter( 1, histZPosAllStripOneEvt.GetMean( ) );
                fitGaus.SetParameter( 2, histZPosAllStripOneEvt.GetRMS( ) );
                histZPosAllStripOneEvt.Fit( &fitGaus, "LMQ", "", 0.0, 160.0 );

                fitGaus.SetParLimits( 1, fitGaus.GetParameter( 1 ) - fitGaus.GetParameter( 2 )*1.2, fitGaus.GetParameter( 1 ) + fitGaus.GetParameter( 2 )*1.2 );
                histZPosAllStripOneEvt.Fit( &fitGaus, "LMQ", "", fitGaus.GetParameter( 1 ) - fitGaus.GetParameter( 2 )*1.2, fitGaus.GetParameter( 1 ) + fitGaus.GetParameter( 2 )*1.2 );
                
                fitGaus.SetParLimits( 1, fitGaus.GetParameter( 1 ) - fitGaus.GetParameter( 2 )*1.0, fitGaus.GetParameter( 1 ) + fitGaus.GetParameter( 2 )*1.0 );
                histZPosAllStripOneEvt.Fit( &fitGaus, "LMQ", "", fitGaus.GetParameter( 1 ) - fitGaus.GetParameter( 2 )*1.0, fitGaus.GetParameter( 1 ) + fitGaus.GetParameter( 2 )*1.0 );

                // fitGaus.SetParLimits( 1, fitGaus.GetParameter( 1 ) - fitGaus.GetParameter( 2 )*1.5, fitGaus.GetParameter( 1 ) + fitGaus.GetParameter( 2 )*1.5 );
                // histZPosAllStripOneEvt.Fit( &fitGaus, "LMQ", "", fitGaus.GetParameter( 1 ) - fitGaus.GetParameter( 2 )*2.0, fitGaus.GetParameter( 1 ) + fitGaus.GetParameter( 2 )*2.0 );
               
                pHistZPosGaus->Fill( fitGaus.GetParameter( 1 ) );

                // z reconstruction for SUM method
                TH1F histWF( "histWF", "histWF", wfBin, wf_time_min, wf_time_max );
                if( pVec_a_hg_mino_peak_wf->size( ) == wfBin ) {
                    for( int i = 0; i < wfBin; ++i ) { 
                        if( histWF.GetBinCenter( i + 1 ) < -1.0 * minority_ROI_offset + 5.0 && histWF.GetBinCenter( i + 1 ) > -1.0 * (minority_ROI_offset + minority_ROI_range) - 5.0 )
                            histWF.SetBinContent( i + 1, pVec_a_hg_mino_peak_wf->at( i ) );
                    }
                }

                TSpectrum tSpec( 10 );
                int numPeaks = tSpec.Search( &histWF, 0.01, "nodrawnobackground", 0.3 );
                // cout << numPeaks << endl;
                double* peakTimes   = tSpec.GetPositionX( );
                double* peakHeights = tSpec.GetPositionY( );
                double  peakMax = 0.0;
                double  timePeakMax = 0.0;
                for( int peaks = 0; peaks < numPeaks; ++peaks ) {
                    // ROI definition
                    if( peakTimes[peaks] > -1.0 * minority_ROI_offset || peakTimes[peaks] < -1.0 * (minority_ROI_offset + minority_ROI_range) ) continue;
                    // cout <<  peakTimes[peaks] << "\t" << peakMax << endl;

                    if( peakHeights[peaks] > peakMax ) {
                    // if( peakTimes[peaks] < timePeakMax ) {
                        timePeakMax = peakTimes[peaks];
                        peakMax = peakHeights[peaks];
                    }
                }
                
                double dtSum = fabs(timePeakMax);
                double absZSum = dtSum / ( (1.0 / driftV_sf6) - (1.0 / driftV_sf5) );

                // cout << dtSum << "\t" << peakMax << "\t" << absZSum << endl;
                // cout << "file ID: " << fileID << ",\teventID: " << ev << ",\tdt: " << dtSum << ",\tabsZ: " << absZSum << endl;

                pHistZPosSum->Fill( absZSum * 10.0 );
            }
        
        }
    }


    int graphNumPlot = inputFileTable.size( );
    // int graphNumPlot = histTable.size( );
    double graphSrcPos[10]               = { };
    double graphSrcPosErr[10]            = { };
    double graphSrcPosMeasure[10]        = { };
    double graphSrcPosMeasureErr[10]     = { };
    double graphSrcPosResoMeasure[10]    = { };
    double graphSrcPosResoMeasureErr[10] = { };
    
    double graphSrcPosAllStrip[10]           = { };
    double graphSrcPosAllStripErr[10]        = { };
    double graphSrcPosAllStripMeasure[10]    = { };
    double graphSrcPosAllStripMeasureErr[10] = { };

    double graphSrcPosGaus[10]               = { };
    double graphSrcPosGausErr[10]            = { };
    double graphSrcPosGausMeasure[10]        = { };
    double graphSrcPosGausMeasureErr[10]     = { };
    double graphSrcPosResoGausMeasure[10]    = { };
    double graphSrcPosResoGausMeasureErr[10] = { };

    double graphSrcPosSum[10]               = { };
    double graphSrcPosSumErr[10]            = { };
    double graphSrcPosSumMeasure[10]        = { };
    double graphSrcPosSumMeasureErr[10]     = { };
    double graphSrcPosResoSumMeasure[10]    = { };
    double graphSrcPosResoSumMeasureErr[10] = { };

    TCanvas cvs( "cvs", "cvs", 800, 600 );
    TF1 fitGaus( "gaus", "gaus", 0.0, 16.0 );
    int idx = 0, idxAllStrip = 0, idxGaus = 0, idxSum = 0;
    for( auto pair : histTable ) {
        TH1F* pHist = pair.first;
        double sourcePos = pair.second;
        if( pHist == nullptr ) continue;

        pHist->SetMarkerStyle( 20 );
        pHist->SetXTitle( "#it{z} [mm]" );
        pHist->SetYTitle( "Events" );
        pHist->GetYaxis()->SetRangeUser(0.0, pHist->GetMaximum( ) * 1.7 );
        pHist->Draw( "ep" );

        std::string histName = pHist->GetName( );
        
        double fitRangeMin = 0.0, fitRangeMax = 0.0;
        std::string saveFileName = "";
        if( histName.find( "hist_sourcePos_" ) != std::string::npos ) {
            fitGaus.SetParameter( 1, sourcePos*0.1 );
            fitGaus.SetParameter( 2, 1.0 );
            if( sourcePos < 100.0 ) {
                fitRangeMin = 55.0;
                fitRangeMax = 120.0;
            }
            else {
                fitRangeMin = 80.0;
                fitRangeMax = 140.0;
            }

            pHist->Fit( &fitGaus, "LMI", "", fitRangeMin, fitRangeMax );

            graphSrcPos[idx] = sourcePos;
            graphSrcPosErr[idx] = 0.000001;
            graphSrcPosMeasure[idx] = fitGaus.GetParameter( 1 );
            graphSrcPosMeasureErr[idx] = fitGaus.GetParError( 1 );
            graphSrcPosResoMeasure[idx] = fitGaus.GetParameter( 2 );
            graphSrcPosResoMeasureErr[idx] = fitGaus.GetParError( 2 );

            ++idx;
            saveFileName = Form( "%s/hist_%d.png", outputDir.c_str( ), static_cast< int >( sourcePos ) );
        }
        else if( histName.find( "hist_sourcePosAllStrip_" ) != std::string::npos ) {
            fitGaus.SetParameter( 1, sourcePos );
            fitGaus.SetParameter( 2, 1.0 );
            fitRangeMin = sourcePos * ( 1.0 - 0.2 );
            fitRangeMax = sourcePos * ( 1.0 + 0.2 );
            pHist->Fit( &fitGaus, "LMI", "", fitRangeMin, fitRangeMax );
            
            graphSrcPosAllStrip[idxAllStrip] = sourcePos;
            graphSrcPosAllStripErr[idxAllStrip] = 0.000001;
            graphSrcPosAllStripMeasure[idxAllStrip] = fitGaus.GetParameter( 1 );
            graphSrcPosAllStripMeasureErr[idxAllStrip] = fitGaus.GetParError( 1 );

            ++idxAllStrip;
            saveFileName = Form( "%s/histAllStrip_%d.png", outputDir.c_str( ), static_cast< int >( sourcePos ) );
        }
        else if( histName.find( "hist_sourcePosGaus_" ) != std::string::npos ) {
            fitGaus.SetParameter( 1, sourcePos );
            fitGaus.SetParameter( 2, 1.0 );
            fitRangeMin = sourcePos * ( 1.0 - 0.2 );
            fitRangeMax = sourcePos * ( 1.0 + 0.2 );
            pHist->Fit( &fitGaus, "LMI", "", fitRangeMin, fitRangeMax );
            
            graphSrcPosGaus[idxGaus] = sourcePos;
            graphSrcPosGausErr[idxGaus] = 0.000001;
            graphSrcPosGausMeasure[idxGaus] = fitGaus.GetParameter( 1 );
            graphSrcPosGausMeasureErr[idxGaus] = fitGaus.GetParError( 1 );
            graphSrcPosResoGausMeasure[idxGaus] = fitGaus.GetParameter( 2 );
            graphSrcPosResoGausMeasureErr[idxGaus] = fitGaus.GetParError( 2 );

            ++idxGaus;
            saveFileName = Form( "%s/histGaus_%d.png", outputDir.c_str( ), static_cast< int >( sourcePos ) );
        }
        else if( histName.find( "hist_sourcePosSum_" ) != std::string::npos ) {
            fitGaus.SetParameter( 1, sourcePos );
            fitGaus.SetParameter( 2, 1.0 );
            fitRangeMin = sourcePos * ( 1.0*0.9 - 0.3 );
            fitRangeMax = sourcePos * ( 1.0*0.9 + 0.3 );
            pHist->Fit( &fitGaus, "LMI", "", fitRangeMin, fitRangeMax );
            
            graphSrcPosSum[idxSum] = sourcePos;
            graphSrcPosSumErr[idxSum] = 0.000001;
            graphSrcPosSumMeasure[idxSum] = fitGaus.GetParameter( 1 );
            graphSrcPosSumMeasureErr[idxSum] = fitGaus.GetParError( 1 );
            graphSrcPosResoSumMeasure[idxSum] = fitGaus.GetParameter( 2 );
            graphSrcPosResoSumMeasureErr[idxSum] = fitGaus.GetParError( 2 );

            ++idxSum;
            saveFileName = Form( "%s/histSum_%d.png", outputDir.c_str( ), static_cast< int >( sourcePos ) );
        }
        else {
            continue;
        }

        CreateDrawText( 0.15, 0.85, "^{252}Cf run", 0.05 );
        CreateDrawText( 0.15, 0.79, "Z projection hit map", 0.05 );
        CreateDrawText( 0.15, 0.71, Form("Source position: %d mm", static_cast< int >( sourcePos ) ), 0.05 );

        CreateDrawText( 0.65, 0.79, Form("Mean: %0.2lf", fitGaus.GetParameter( 1 ) ), 0.05 );
        CreateDrawText( 0.65, 0.73, Form("Sigma: %0.2lf", fitGaus.GetParameter( 2 ) ), 0.05 );

        cvs.SaveAs( saveFileName.c_str( ) );
    }

    
    TGraphErrors zPosGraph( graphNumPlot, graphSrcPos, graphSrcPosMeasure, graphSrcPosErr, graphSrcPosMeasureErr );
    zPosGraph.SetMarkerStyle( 20 );
    zPosGraph.SetLineWidth( 2 );
    zPosGraph.Draw( "ap" );
    zPosGraph.GetXaxis( )->SetTitle( "^{241}Am source position [mm]" );
    zPosGraph.GetXaxis( )->SetRangeUser( 0.0, 140.0 );
    zPosGraph.GetXaxis( )->SetLimits( 0.0, 140.0 );
    zPosGraph.GetYaxis( )->SetTitle( "Reconstructed track position [mm]" );
    zPosGraph.GetYaxis( )->SetRangeUser( 0.0, 140.0 );
    cvs.SaveAs( Form( "%s/zpos.png", outputDir.c_str( ) ) );

    TGraphErrors zPosResoGraph( graphNumPlot, graphSrcPos, graphSrcPosResoMeasure, graphSrcPosErr, graphSrcPosResoMeasureErr );
    zPosResoGraph.SetMarkerStyle( 20 );
    zPosResoGraph.SetLineWidth( 2 );
    zPosResoGraph.GetXaxis( )->SetTitle( "^{241}Am source position [mm]" );
    // zPosResoGraph.GetXaxis( )->SetRangeUser( 0.0, 140.0 );
    // zPosResoGraph.GetXaxis( )->SetLimits( 0.0, 140.0 );
    zPosResoGraph.GetYaxis( )->SetTitle( "Reconstructed track position resolution [mm]" );
    // zPosResoGraph.GetYaxis( )->SetRangeUser( 0.0, 40.0 );
    zPosResoGraph.Draw( "ap" );
    cvs.SaveAs( Form( "%s/zposreso.png", outputDir.c_str( ) ) );

    TGraphErrors zPosAllStripGraph( graphNumPlot, graphSrcPosAllStrip, graphSrcPosAllStripMeasure, graphSrcPosAllStripErr, graphSrcPosAllStripMeasureErr );
    zPosAllStripGraph.SetMarkerStyle( 20 );
    zPosAllStripGraph.SetLineWidth( 2 );
    zPosAllStripGraph.Draw( "ap" );
    zPosAllStripGraph.GetXaxis( )->SetTitle( "^{241}Am source position [mm]" );
    zPosAllStripGraph.GetYaxis( )->SetTitle( "Reconstructed track position [mm]" );
    zPosAllStripGraph.GetXaxis( )->SetLimits( 0.0, 140.0 );
    zPosAllStripGraph.GetYaxis( )->SetRangeUser( 0.0, 140.0 );
    zPosAllStripGraph.Draw( "ap" );
    cvs.SaveAs( Form( "%s/zposAllStrip.png", outputDir.c_str( ) ) );


    TGraphErrors zPosGausGraph( graphNumPlot, graphSrcPosGaus, graphSrcPosGausMeasure, graphSrcPosGausErr, graphSrcPosGausMeasureErr );
    zPosGausGraph.SetMarkerStyle( 20 );
    zPosGausGraph.SetLineWidth( 2 );
    zPosGausGraph.Draw( "ap" );
    zPosGausGraph.GetXaxis( )->SetTitle( "^{241}Am source position [mm]" );
    zPosGausGraph.GetXaxis( )->SetRangeUser( 0.0, 140.0 );
    zPosGausGraph.GetXaxis( )->SetLimits( 0.0, 140.0 );
    zPosGausGraph.GetYaxis( )->SetTitle( "Reconstructed track position [mm]" );
    zPosGausGraph.GetYaxis( )->SetRangeUser( 0.0, 140.0 );
    cvs.SaveAs( Form( "%s/zposGaus.png", outputDir.c_str( ) ) );

    TGraphErrors zPosResoGausGraph( graphNumPlot, graphSrcPosGaus, graphSrcPosResoGausMeasure, graphSrcPosGausErr, graphSrcPosResoGausMeasureErr );
    zPosResoGausGraph.SetMarkerStyle( 20 );
    zPosResoGausGraph.SetLineWidth( 2 );
    zPosResoGausGraph.GetXaxis( )->SetTitle( "^{241}Am source position [mm]" );
    zPosResoGausGraph.GetXaxis( )->SetRangeUser( 0.0, 140.0 );
    zPosResoGausGraph.GetXaxis( )->SetLimits( 0.0, 140.0 );
    zPosResoGausGraph.GetYaxis( )->SetTitle( "Reconstructed track position resolution [mm]" );
    zPosResoGausGraph.GetYaxis( )->SetRangeUser( 0.0, 40.0 );
    zPosResoGausGraph.Draw( "ap" );
    cvs.SaveAs( Form( "%s/zposresoGaus.png", outputDir.c_str( ) ) );


    TGraphErrors zPosSumGraph( graphNumPlot, graphSrcPosSum, graphSrcPosSumMeasure, graphSrcPosSumErr, graphSrcPosSumMeasureErr );
    zPosSumGraph.SetMarkerStyle( 20 );
    zPosSumGraph.SetLineWidth( 2 );
    zPosSumGraph.Draw( "ap" );
    zPosSumGraph.GetXaxis( )->SetTitle( "^{241}Am source position [mm]" );
    zPosSumGraph.GetXaxis( )->SetRangeUser( 0.0, 140.0 );
    zPosSumGraph.GetXaxis( )->SetLimits( 0.0, 140.0 );
    zPosSumGraph.GetYaxis( )->SetTitle( "Reconstructed track position [mm]" );
    zPosSumGraph.GetYaxis( )->SetRangeUser( 0.0, 140.0 );
    zPosSumGraph.Draw( "ap" );

    TLine linearLine( 0.0, 0.0, 140.0, 140.0 );
    linearLine.SetLineWidth( 1 );
    linearLine.SetLineColor( kBlue+2 );
    linearLine.SetLineStyle( 2 );
    linearLine.Draw( );
    
    cvs.SaveAs( Form( "%s/zposSum.png", outputDir.c_str( ) ) );

    TGraphErrors zPosResoSumGraph( graphNumPlot, graphSrcPosSum, graphSrcPosResoSumMeasure, graphSrcPosSumErr, graphSrcPosResoSumMeasureErr );
    zPosResoSumGraph.SetMarkerStyle( 20 );
    zPosResoSumGraph.SetLineWidth( 2 );
    zPosResoSumGraph.GetXaxis( )->SetTitle( "^{241}Am source position [mm]" );
    zPosResoSumGraph.GetXaxis( )->SetRangeUser( 0.0, 140.0 );
    zPosResoSumGraph.GetXaxis( )->SetLimits( 0.0, 140.0 );
    zPosResoSumGraph.GetYaxis( )->SetTitle( "Reconstructed track position resolution [mm]" );
    zPosResoSumGraph.GetYaxis( )->SetRangeUser( 0.0, 40.0 );
    zPosResoSumGraph.Draw( "ap" );
    cvs.SaveAs( Form( "%s/zposresoSum.png", outputDir.c_str( ) ) );


    return;
}
