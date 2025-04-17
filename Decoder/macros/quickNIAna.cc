#include "NAPStyle.h"
#include "NIConfig.hh"
#include "NIConfig.cc"
#include "NAUtil.h"
#include "NAUtil.cc"

#include "TEfficiency.h" 

const double SIM_ENERGY_DEPOSIT = 128.0; // keV (3 board DAQ: 1.28 cm)

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

void quickNIAna( const std::string& inputFile, const std::string& configFile, const std::string& outputDir )
{
    // file check
    TFile file( inputFile.c_str( ) );
    TTree* pTree = dynamic_cast< TTree* >( file.Get( "ana_tree" ) );
    if( pTree == nullptr ) {
        NAUtil::Cerr( "failed to read ROOT file..." );
        return;
    }

    NIConfig* pConfig = new NIConfig( );
    if( pConfig->ReadConfigJSON( configFile ) == false ) {
        NAUtil::Cerr( "failed to read CONFIG file..." );
        return;
    }

    pConfig->PrintConfigJSON( );

    NAUtil::ExistCreateDir( outputDir );
    SetNAPStyle( );

    // analysis

    bool   isAlphaCalib = pConfig->is_alpha_calib == 0 ? false : true;
    double calFactor    = pConfig->cal_factor;
    double driftV_sf6   = pConfig->driftV_main;
    double driftV_sf5   = pConfig->driftV_mino;

    int ev = 0;
    int fileID = 0;
    vector< double >* pVec_xz_x = nullptr;
    vector< double >* pVec_xz_z = nullptr;
    vector< double >* pVec_yz_y = nullptr;
    vector< double >* pVec_yz_z = nullptr;
    double ave_x = 0.0, ave_y = 0.0, ave_z = 0.0;
    double a_hg_sum_pulse_height = 0.0, a_hg_sum_charge = 0.0;
    double c_hg_sum_pulse_height = 0.0, c_hg_sum_charge = 0.0;
    vector< double >* pVec_a_hg_mainrise_time = nullptr;
    vector< double >* pVec_a_hg_mainfall_time = nullptr;
    vector< double >* pVec_c_hg_mainrise_time = nullptr;
    vector< double >* pVec_c_hg_mainfall_time = nullptr;
    vector< double >* pVec_dt = nullptr;

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
    pTree->SetBranchAddress( "a_hg_mainrise_time", &pVec_a_hg_mainrise_time );
    pTree->SetBranchAddress( "a_hg_mainfall_time", &pVec_a_hg_mainfall_time );
    pTree->SetBranchAddress( "c_hg_mainrise_time", &pVec_c_hg_mainrise_time );
    pTree->SetBranchAddress( "c_hg_mainfall_time", &pVec_c_hg_mainfall_time );
    pTree->SetBranchAddress( "dt",                 &pVec_dt );

    TH1D* pHistEnergy    = new TH1D( "histEnergy", "histEnergy", 80, 0.0, 400.0 );
    TH1D* pHistEnergyPulse = new TH1D( "histEnergyPulse", "histEnergyPulse", 80, 0.0, 400.0 );
    TH1D* pHistEnergyCut = new TH1D( "histEnergyCut", "histEnergyCut", 80, 0.0, 400.0 );

    TH2D* pHistHitmapXY = new TH2D( "histHitmapXY", "histHitmapXY", 100, -1.5, 1.5, 100, -1.5,  1.5 );
    TH2D* pHistHitmapXZ = new TH2D( "histHitmapXZ", "histHitmapXZ", 100, -1.5, 1.5, 100,  0.0, 12.0 );
    TH2D* pHistHitmapYZ = new TH2D( "histHitmapYZ", "histHitmapYZ", 100, -1.5, 1.5, 100,  0.0, 12.0 );

    TH2D* pHistHitmapXYLenCut = new TH2D( "histHitmapXYLenCut", "histHitmapXYLenCut", 100, -1.5, 1.5, 100, -1.5,  1.5 );

    TH2D* pHistHitmapXMaxMin = new TH2D( "histHitmapXMaxMin", "histHitmapXMaxMin", 100, -1.5, 1.5, 100, -1.5,  1.5 );
    TH2D* pHistHitmapYMaxMin = new TH2D( "histHitmapYMaxMin", "histHitmapYMaxMin", 100, -1.5, 1.5, 100, -1.5,  1.5 );

    TH2D* pHistHitmapXMaxMinCut = new TH2D( "histHitmapXMaxMinCut", "histHitmapXMaxMinCut", 100, -1.5, 1.5, 100, -1.5,  1.5 );
    TH2D* pHistHitmapYMaxMinCut = new TH2D( "histHitmapYMaxMinCut", "histHitmapYMaxMinCut", 100, -1.5, 1.5, 100, -1.5,  1.5 );

    TH2D* pHistEnergyToT    = new TH2D( "histEnergyToT", "histEnergyToT", 100, 0.0, 400.0, 100, 0.0, 2000.0 );
    TH2D* pHistEnergyToTCut = new TH2D( "histEnergyToTCut",   "histEnergyToTCut",   100, 0.0, 400.0, 100, 0.0, 2000.0 );

    TH2D* pHistEnergyToTEn    = new TH2D( "histEnergyToTEn", "histEnergyToTEn", 100, 0.0, 400.0, 100, 0.0, 10.0 );
    TH2D* pHistEnergyToTEnCut = new TH2D( "histEnergyToTEnCut", "histEnergyToTEnCut", 100, 0.0, 400.0, 100, 0.0, 10.0 );

    TH2D* pHistEnergyLengthXZ = new TH2D( "histEnergyLengthXZ", "histEnergyLengthXZ", 40, 0.0, 400.0, 40, 0.0, 2.0 );
    TH2D* pHistEnergyLengthYZ = new TH2D( "histEnergyLengthYZ", "histEnergyLengthYZ", 40, 0.0, 400.0, 40, 0.0, 2.0 );
    TH2D* pHistEnergyLengthXY = new TH2D( "histEnergyLengthXY", "histEnergyLengthXY", 40, 0.0, 400.0, 40, 0.0, 2.0 );

    TH2D* pHistEnergyLengthXZCut = new TH2D( "histEnergyLengthXZCut", "histEnergyLengthXZCut", 40, 0.0, 400.0, 40, 0.0, 2.0 );
    TH2D* pHistEnergyLengthYZCut = new TH2D( "histEnergyLengthYZCut", "histEnergyLengthYZCut", 40, 0.0, 400.0, 40, 0.0, 2.0 );

    TH2D* pHistHitmapXYCut = new TH2D( "histHitmapXYCut", "histHitmapXYCut", 100, -1.5, 1.5, 100, -1.5,  1.5 );

    TH2D* pHistEnergyDt = new TH2D( "histEnergyDt", "histEnergyDt", 100, 0.0, 400.0, 100, 0.0, 500.0 );

    TH3D* pHistHitmapXYZ = new TH3D( "histHitmapXYZ", "histHitmapXYZ", 150, -0.75, 2.25, 150, -1.5,  1.5, 150, 0.0, 14.4 ); 
    TH2D* pHistHitmapAllSelXY = new TH2D( "histHitmapAllSelXY", "histHitmapAllSelXY", 1000, -0.75, 2.25, 1000, -1.5,  2.0 );
    TH2D* pHistHitmapAllSelXZ = new TH2D( "histHitmapAllSelXZ", "histHitmapAllSelXZ", 1000, -0.75, 2.25, 1000, 0.0,  16.0 );
    TH2D* pHistHitmapAllSelYZ = new TH2D( "histHitmapAllSelYZ", "histHitmapAllSelYZ", 1000, -1.5, 1.5, 1000, 0.0,  16.0 );
    TH1D* pHistHitmapAllSelZ  = new TH1D( "histHitmapZ", "histHitmapZ", 16, 0.0, 16.0 );

    TH2D* pHistRiseMeanSigma = new TH2D( "histRiseMeanSigma","histRiseMeanSigma",100,1000,1200,100,0,20 );
    TH2D* pHistRiseSigmaY    = new TH2D( "histRiseSigmaY","histRiseSigmaY",100,-1.5,1.5,100,0,20 );

    TH1D* pHistAllDt = new TH1D( "histAllDt","histAllDt",100, 0.0, 200.0 );

    ifstream dataF("SRIM_F_in_20torr_SF6.dat");
    double SRIMenergy_F[1000], SRIMlength_F[1000];
    int SRIMnum_F=0;
    while(dataF>>SRIMenergy_F[SRIMnum_F]>>SRIMlength_F[SRIMnum_F]) SRIMnum_F++;
    dataF.close();
    for(int i=0;i<79;i++) SRIMlength_F[i]=SRIMlength_F[i]*0.1;
    TGraph *g_SRIM_F = new TGraph(SRIMnum_F, SRIMenergy_F, SRIMlength_F);
    g_SRIM_F->SetName("g_SRIM_F");
    g_SRIM_F->SetTitle("g_SRIM_F");
    g_SRIM_F->SetLineColor(1);
    g_SRIM_F->SetLineWidth(2);

    ifstream dataHe("SRIM_He_in_20torr_SF6.dat");
    double SRIMenergy_He[1000], SRIMlength_He[1000];
    int SRIMnum_He=0;
    while(dataHe>>SRIMenergy_He[SRIMnum_He]>>SRIMlength_He[SRIMnum_He]) SRIMnum_He++;
    dataHe.close();
    for(int i=0;i<79;i++) SRIMlength_He[i]=SRIMlength_He[i]*0.1;
    TGraph *g_SRIM_He = new TGraph(SRIMnum_He, SRIMenergy_He, SRIMlength_He);
    g_SRIM_He->SetName("g_SRIM_He");
    g_SRIM_He->SetTitle("g_SRIM_He");
    g_SRIM_He->SetLineColor(2);
    g_SRIM_He->SetLineWidth(2);

    ifstream dataH("SRIM_H_in_20torr_SF6.dat");
    double SRIMenergy_H[1000], SRIMlength_H[1000];
    int SRIMnum_H=0;
    while(dataH>>SRIMenergy_H[SRIMnum_H]>>SRIMlength_H[SRIMnum_H]) SRIMnum_H++;
    dataH.close();
    for(int i=0;i<79;i++) SRIMlength_H[i]=SRIMlength_H[i]*0.1;
    TGraph *g_SRIM_H = new TGraph(SRIMnum_H, SRIMenergy_H, SRIMlength_H);
    g_SRIM_H->SetName("g_SRIM_H");
    g_SRIM_H->SetTitle("g_SRIM_H");
    g_SRIM_H->SetLineColor(4);
    g_SRIM_H->SetLineWidth(2);

    int totEvt = pTree->GetEntries( );
    int selEvtAll = 0, selEvt50 = 0, selEvt100 = 0, selEvt200 = 0, selEvt400 = 0;
    int selEvtAllMino = 0, selEvt50Mino = 0, selEvt100Mino = 0, selEvt200Mino = 0, selEvt400Mino = 0;

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

        // NAUtil::PrintProgressBar( evt, totEvt );
        pTree->GetEntry( evt );
        
        if( pVec_yz_y->size( ) <= 0 ) continue;

        pHistHitmapXY->Fill( ave_x, ave_y );
        pHistHitmapXZ->Fill( ave_x, ave_z );
        pHistHitmapYZ->Fill( ave_y, ave_z );


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

        pHistHitmapXMaxMin->Fill( xz_x_max, xz_x_min );
        pHistHitmapYMaxMin->Fill( yz_y_max, yz_y_min );

        // fiducial cut
        if( isAlphaCalib == true ) {
            if( xz_x_max < 1.1 || xz_x_min > 0.0 ) continue;
        }
        else {
            if( xz_x_max > 1.1 || xz_x_min < 0.1 ) continue;
            if( yz_y_max > 1.1 || yz_y_min < -1.1 ) continue;
        }

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
        if( isAlphaCalib == true ) {
            pHistEnergyToTCut->Fill( a_hg_sum_charge * calFactor, sumToT );
            pHistEnergyToTEnCut->Fill( a_hg_sum_charge * calFactor, sumToT/(a_hg_sum_charge * calFactor) );
            pHistEnergyLengthXZCut->Fill( a_hg_sum_charge * calFactor, lengthXZ );
            pHistEnergyLengthYZCut->Fill( a_hg_sum_charge * calFactor, lengthYZ );
            if( lengthXZ < 1.3 ) {
                // at least 2 hits are required in cathode strips
                if( pVec_yz_y->size( ) <= 1 ) continue;
                // if( mainrise_sigma > 8.0 ) continue;
                if( mainrise_sigma > 4.0 ) continue;

                pHistEnergyCut->Fill( a_hg_sum_charge * calFactor );
                pHistEnergyPulse->Fill( a_hg_sum_pulse_height * 127.8 / 16445.0 );

                double averageDt = 0.0;
                if( pVec_dt->size( ) > 0 ) averageDt = std::accumulate( pVec_dt->begin(),pVec_dt->end(),0.0)/pVec_dt->size();
                pHistEnergyDt->Fill( a_hg_sum_charge * calFactor, averageDt );

                for( auto eachDt : *pVec_dt ) pHistAllDt->Fill( eachDt );

                double absZ = averageDt / ( (1.0 / driftV_sf6) - (1.0 / driftV_sf5) );
                pHistHitmapXYZ->Fill( ave_x, ave_y, absZ );
                pHistHitmapAllSelXY->Fill( ave_x, ave_y );
                pHistHitmapAllSelXZ->Fill( ave_x, absZ );
                pHistHitmapAllSelYZ->Fill( ave_y, absZ );
                pHistHitmapAllSelZ->Fill( absZ );

                std::cout << "file ID: " << fileID << "\tevtID: " << ev << "\tenergy: " << a_hg_sum_charge * calFactor << "\tdt: " << absZ << "\tX: " << ave_x << "\tY: " <<  ave_y << std::endl;

                pHistRiseMeanSigma->Fill( mainrise_mean, mainrise_sigma );
                pHistRiseSigmaY->Fill( ave_y, mainrise_sigma );
            }

        }
        else {
            if( lengthXZ < 0.004 * a_hg_sum_charge * calFactor + 0.4 ) {
                pHistEnergyToTCut->Fill( a_hg_sum_charge * calFactor, sumToT );
                pHistEnergyToTEnCut->Fill( a_hg_sum_charge * calFactor, sumToT/(a_hg_sum_charge * calFactor) );
                pHistEnergyLengthXZCut->Fill( a_hg_sum_charge * calFactor, lengthXZ );
                pHistEnergyLengthYZCut->Fill( a_hg_sum_charge * calFactor, lengthYZ );
                if( a_hg_sum_charge * calFactor > 30.0 && a_hg_sum_charge * calFactor < 400.0 ) {
                    ++selEvtAll;
                    if     ( a_hg_sum_charge * calFactor < 50.0  ) ++selEvt50;
                    else if( a_hg_sum_charge * calFactor < 100.0 ) ++selEvt100;
                    else if( a_hg_sum_charge * calFactor < 200.0 ) ++selEvt200;
                    else if( a_hg_sum_charge * calFactor < 400.0 ) ++selEvt400;

                    // std::cout << std::endl << "file ID: " << fileID << "\tevtID: " << ev << std::endl;
                    // std::cout << "file ID: " << fileID << "\tevtID: " << ev << "\tenergy: " << a_hg_sum_charge * calFactor << std::endl;

                    // if( ( lengthXZ > 1.0 || lengthYZ > 1.0 ) && ( lengthXZ < 1.3 || lengthYZ < 1.3 ) ) {

                    pHistEnergyCut->Fill( a_hg_sum_charge * calFactor );
                    pHistHitmapXYLenCut->Fill( ave_x, ave_y );
                    double averageDt = 0.0;
                    if( pVec_dt->size( ) > 0 ) {
                        averageDt = std::accumulate( pVec_dt->begin(),pVec_dt->end(),0.0)/pVec_dt->size();
                        ++selEvtAllMino;
                        if     ( a_hg_sum_charge * calFactor < 50.0  ) ++selEvt50Mino;
                        else if( a_hg_sum_charge * calFactor < 100.0 ) ++selEvt100Mino;
                        else if( a_hg_sum_charge * calFactor < 200.0 ) ++selEvt200Mino;
                        else if( a_hg_sum_charge * calFactor < 400.0 ) ++selEvt400Mino;

                        // if( averageDt < 120 )
                        //     std::cout << "file ID: " << fileID << "\tevtID: " << ev << "\tenergy: " << a_hg_sum_charge * calFactor << "\tdt: " << averageDt << std::endl;

                    }
                    else {
                        // std::cout << "file ID: " << fileID << "\tevtID: " << ev << "\tenergy: " << a_hg_sum_charge * calFactor << std::endl;
                    }
                
                
                    pHistEnergyDt->Fill( a_hg_sum_charge * calFactor, averageDt );

                    double absZ = averageDt / ( (1.0 / driftV_sf6) - (1.0 / driftV_sf5) );
                    // std::cout << absZ << std::endl;
                    if( absZ > 0.1 && absZ < 16.0 ) {
                        pHistHitmapXYZ->Fill( ave_x, ave_y, absZ );
                        pHistHitmapAllSelXY->Fill( ave_x, ave_y );
                        pHistHitmapAllSelXZ->Fill( ave_x, absZ );
                        pHistHitmapAllSelYZ->Fill( ave_y, absZ );
                        pHistHitmapAllSelZ->Fill( absZ );
                        if( absZ > 14.0 )
                        // if( a_hg_sum_charge * calFactor < 100.0 && a_hg_sum_charge * calFactor > 50.0 )
                            std::cout << "file ID: " << fileID << "\tevtID: " << ev << "\tenergy: " << a_hg_sum_charge * calFactor << "\tdt: " << absZ << "\tX: " << ave_x << "\tY: " <<  ave_y << std::endl;

                    }
                }
            }
        }

        pHistHitmapXYCut->Fill( ave_x, ave_y );
        pHistHitmapXMaxMinCut->Fill( xz_x_max, xz_x_min );
        pHistHitmapYMaxMinCut->Fill( yz_y_max, yz_y_min );

        pHistEnergy->Fill( a_hg_sum_charge * calFactor );

        pHistEnergyLengthXZ->Fill( a_hg_sum_charge * calFactor, lengthXZ );
        pHistEnergyLengthXY->Fill( a_hg_sum_charge * calFactor, lengthXY );
        pHistEnergyLengthYZ->Fill( a_hg_sum_charge * calFactor, lengthYZ );

        pHistEnergyToT->Fill( a_hg_sum_charge * calFactor, sumToT );
        pHistEnergyToTEn->Fill( a_hg_sum_charge * calFactor, sumToT/(a_hg_sum_charge * calFactor) );
    }
    
    TCanvas cvs( "cvs", "cvs", 800, 600 );
    // TCanvas cvs( "cvs", "cvs", 800, 800 );

    gPad->SetRightMargin( 0.1 );

    double zMean = 0.0, zRes = 0.0;
    if( isAlphaCalib == true ) {
        TF1 zResGauss( "zResGauss", "gaus", 5, 12 );
        pHistHitmapAllSelZ->Fit( &zResGauss, "LMI", "", 5, 12 );
        zMean = zResGauss.GetParameter( 1 );
        zRes  = zResGauss.GetParameter( 2 );
    }

    pHistHitmapAllSelZ->SetMarkerStyle(20);
    pHistHitmapAllSelZ->GetXaxis()->SetTitle( "#it{z} [cm]" );
    pHistHitmapAllSelZ->GetYaxis()->SetTitle( "Events" );
    pHistHitmapAllSelZ->GetYaxis()->SetRangeUser(0.0, pHistHitmapAllSelZ->GetMaximum( ) * 1.7 );
    pHistHitmapAllSelZ->Draw("ep");
    CreateDrawText( 0.55, 0.85, "^{252}Cf run", 0.05 );
    CreateDrawText( 0.55, 0.79, "Z projection hit map", 0.05 );

    TLine fidLineZ( 14.4, 0.0, 14.4, 10.5 );
    fidLineZ.SetLineWidth( 2 );
    fidLineZ.SetLineColor( kBlue );
    // fidLineZ.SetLineStyle( 3 );
    fidLineZ.Draw( );
    TArrow fidArrowZ(14.4, 9.0, 13.0, 9.0, 0.02, "|>");
    fidArrowZ.SetLineWidth( 2 );
    fidArrowZ.SetLineColor( kBlue );
    fidArrowZ.SetFillColor( kBlue );
    fidArrowZ.Draw( );
    if( isAlphaCalib == true ) {
        CreateDrawText( 0.6, 0.65, Form("Mean: %0.2lf [cm]", zMean ), 0.045 );
        CreateDrawText( 0.6, 0.58, Form("Sigma: %0.2lf [cm]", zRes ), 0.045 );
    }

    cvs.SaveAs( Form( "%s/absZ.png", outputDir.c_str( ) ) );
    cvs.SaveAs( Form( "%s/absZ.pdf", outputDir.c_str( ) ) );

    // pHistEnergy->GetXaxis()->SetTitle( "#SigmaADC count" );
    pHistEnergy->GetXaxis()->SetTitle( "Energy [keV]" );
    pHistEnergy->GetYaxis()->SetTitle( "Events / 5 keV" );
    pHistEnergy->GetYaxis()->SetRangeUser( 0.0, pHistEnergy->GetMaximum( ) * 1.5 );

    pHistEnergy->SetLineColor( kBlack );
    pHistEnergy->SetLineWidth( 2 );
    pHistEnergy->SetMarkerColor( kBlack );
    pHistEnergy->SetMarkerStyle( 20 );
    pHistEnergyCut->SetLineColor( kRed );
    pHistEnergyCut->SetLineWidth( 2 );
    pHistEnergyCut->SetMarkerColor( kRed );
    pHistEnergyCut->SetMarkerStyle( 20 );

    pHistEnergy->Draw( "e" );
    pHistEnergyCut->Draw( "esame" );

    TF1 fitGauss( "fitGauss", "gaus", 70, 140 );
    fitGauss.SetLineColor( kRed );
    pHistEnergyCut->Fit(&fitGauss,"LMI","",70, 140);
    double peakVal = fitGauss.GetParameter( 1 );
    double resoVal = fitGauss.GetParameter( 2 );
    // CreateDrawText( 0.48, 0.85, Form( "Default cal. factor = %0.6lf", calFactor ), 0.035 );
    // CreateDrawText( 0.48, 0.8, Form( "New cal. factor = %0.6lf", calFactor * SIM_ENERGY_DEPOSIT / peakVal ), 0.035 );
    double detReso = 100.0 * sqrt(resoVal*resoVal - 8.33*8.33) / peakVal; // percent
    
    TLegend* pLegEn = new TLegend( 0.6, 0.55, 0.9, 0.7 );
    pLegEn->SetFillStyle( 0 );
    pLegEn->SetBorderSize( 0 );
    pLegEn->SetTextFont( 42 );
    pLegEn->AddEntry( pHistEnergy, "Before cut", "lep" );
    pLegEn->AddEntry( pHistEnergyCut, "After cut", "lep" );
    pLegEn->Draw( );
    
    CreateDrawText( 0.15, 0.88, "^{241}Am run", 0.04);
    CreateDrawText( 0.15, 0.81, Form( "Gas gain: %0.0lf", 2800.0 ), 0.04 );
    // CreateDrawText( 0.15, 0.76, Form( "Energy resolution: %0.1lf %%", 20.0 ), 0.04 );
    CreateDrawText( 0.15, 0.76, Form( "#sigma_{E} / E = %0.0lf%% @128 keV", detReso ), 0.04 );
    cvs.SaveAs( Form( "%s/energy.pdf", outputDir.c_str( ) ) );
    cvs.SaveAs( Form( "%s/energy.png", outputDir.c_str( ) ) );

    pHistEnergyPulse->Draw();
    TF1 fitGaussPulse( "fitGaussPulse", "gaus", 90, 160 );
    pHistEnergyPulse->Fit( &fitGaussPulse, "LMI", "", 90, 160 );
    cvs.SaveAs( Form( "%s/energyPulse.pdf", outputDir.c_str( ) ) );
    
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 1.00 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.80, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable( NRGBs, stops, red, green, blue, NCont );
    gStyle->SetNumberContours( NCont );
    gPad->SetRightMargin( 0.2 );
    cvs.SetGridx( 1 ); cvs.SetGridy( 1 );

    pHistHitmapXY->SetXTitle( "X [cm]" );
    pHistHitmapXY->SetYTitle( "Y [cm]" );
    pHistHitmapXY->SetZTitle( "Events" );
    pHistHitmapXY->Draw( "colz" );
    cvs.SaveAs( Form( "%s/hitmapXY.png", outputDir.c_str( ) ) );

    pHistHitmapXYCut->SetXTitle( "X [cm]" );
    pHistHitmapXYCut->SetYTitle( "Y [cm]" );
    pHistHitmapXYCut->SetZTitle( "Events" );
    pHistHitmapXYCut->Draw( "colz" );
    cvs.SaveAs( Form( "%s/hitmapXYCut.png", outputDir.c_str( ) ) );

    pHistHitmapXYLenCut->SetXTitle( "X [cm]" );
    pHistHitmapXYLenCut->SetYTitle( "Y [cm]" );
    pHistHitmapXYLenCut->SetZTitle( "Events" );
    pHistHitmapXYLenCut->Draw( "colz" );
    cvs.SaveAs( Form( "%s/hitmapXYLenCut.png", outputDir.c_str( ) ) );

    pHistHitmapXZ->SetXTitle( "X [cm]" );
    pHistHitmapXZ->SetYTitle( "Z [cm]" );
    pHistHitmapXZ->SetZTitle( "Events" );
    pHistHitmapXZ->Draw( "colz" );
    cvs.SaveAs( Form( "%s/hitmapXZ.png", outputDir.c_str( ) ) );

    pHistHitmapYZ->SetXTitle( "Y [cm]" );
    pHistHitmapYZ->SetYTitle( "Z [cm]" );
    pHistHitmapYZ->SetZTitle( "Events" );
    pHistHitmapYZ->Draw( "colz" );
    cvs.SaveAs( Form( "%s/hitmapYZ.png", outputDir.c_str( ) ) );

    pHistHitmapXMaxMin->SetXTitle( "X max [cm]" );
    pHistHitmapXMaxMin->SetYTitle( "X min [cm]" );
    pHistHitmapXMaxMin->SetZTitle( "Events" );
    pHistHitmapXMaxMin->Draw( "colz" );
    cvs.SaveAs( Form( "%s/hitmapXMaxMin.png", outputDir.c_str( ) ) );

    pHistHitmapYMaxMin->SetXTitle( "Y max [cm]" );
    pHistHitmapYMaxMin->SetYTitle( "Y min [cm]" );
    pHistHitmapYMaxMin->SetZTitle( "Events" );
    pHistHitmapYMaxMin->Draw( "colz" );
    cvs.SaveAs( Form( "%s/hitmapYMaxMin.png", outputDir.c_str( ) ) );

    pHistHitmapXMaxMinCut->SetXTitle( "X max [cm]" );
    pHistHitmapXMaxMinCut->SetYTitle( "X min [cm]" );
    pHistHitmapXMaxMinCut->SetZTitle( "Events" );
    pHistHitmapXMaxMinCut->Draw( "colz" );
    cvs.SaveAs( Form( "%s/hitmapXMaxMinCut.png", outputDir.c_str( ) ) );

    pHistHitmapYMaxMinCut->SetXTitle( "Y max [cm]" );
    pHistHitmapYMaxMinCut->SetYTitle( "Y min [cm]" );
    pHistHitmapYMaxMinCut->SetZTitle( "Events" );
    pHistHitmapYMaxMinCut->Draw( "colz" );
    cvs.SaveAs( Form( "%s/hitmapYMaxMinCut.png", outputDir.c_str( ) ) );


    pHistEnergyToT->SetXTitle( "Energy [keV]" );
    pHistEnergyToT->SetYTitle( "#SigmaToT [us]" );
    pHistEnergyToT->SetZTitle( "Events" );
    pHistEnergyToT->Draw( "colz" );
    cvs.SaveAs( Form( "%s/energyToT.png", outputDir.c_str( ) ) );

    pHistEnergyToTEn->SetXTitle( "Energy [keV]" );
    pHistEnergyToTEn->SetYTitle( "#SigmaToT / Energy [us/keV]" );
    pHistEnergyToTEn->SetZTitle( "Events" );
    pHistEnergyToTEn->Draw( "colz" );
    cvs.SaveAs( Form( "%s/energyToTEn.png", outputDir.c_str( ) ) );

    pHistEnergyDt->SetXTitle( "Energy [keV]" );
    pHistEnergyDt->SetYTitle( "dt [us]" );
    pHistEnergyDt->SetZTitle( "Events" );
    pHistEnergyDt->Draw( "colz" );
    cvs.SaveAs( Form( "%s/energyDt.png", outputDir.c_str( ) ) );

    pHistEnergyLengthXZ->SetXTitle( "Energy [keV]" );
    pHistEnergyLengthXZ->SetYTitle( "Length [cm]" );
    pHistEnergyLengthXZ->SetZTitle( "Events / 10 keV / 0.05 cm" );
    pHistEnergyLengthXZ->GetZaxis( )->SetTitleOffset(1.5);
    pHistEnergyLengthXZ->Draw( "colz" );
    // g_SRIM_F->Draw( "same" );
    // g_SRIM_He->Draw( "same" );
    // g_SRIM_H->Draw( "same" );

    TF1 cutFunc( "cutFunc", "[0] + [1] * x", 0.0, 400.0 );
    cutFunc.SetParameter( 0, 0.4 );
    cutFunc.SetParameter( 1, 0.004 );
    cutFunc.SetLineColor( kRed );
    cutFunc.SetLineWidth( 4 );
    
    TLine cutLine( 0.0, 1.3, 400.0, 1.3 );
    cutLine.SetLineWidth( 2 );
    cutLine.SetLineColor( kRed );

    if( isAlphaCalib == true ) {
        CreateDrawText( 0.55, 0.2, "^{241}Am run", 0.05 );
        cutLine.Draw( );
    }
    else {
        CreateDrawText( 0.6, 0.2, "^{252}Cf run", 0.06 );
        cutFunc.Draw( "same" );
    }

    cvs.SaveAs( Form( "%s/energyLengthXZ.pdf", outputDir.c_str( ) ) );
    cvs.SaveAs( Form( "%s/energyLengthXZ.png", outputDir.c_str( ) ) );

    pHistEnergyLengthYZ->SetXTitle( "Energy [keV]" );
    pHistEnergyLengthYZ->SetYTitle( "Length [cm]" );
    pHistEnergyLengthYZ->SetZTitle( "Events" );
    pHistEnergyLengthYZ->Draw( "colz" );
    g_SRIM_F->Draw( "same" );
    // g_SRIM_He->Draw( "same" );
    // g_SRIM_H->Draw( "same" );
    cvs.SaveAs( Form( "%s/energyLengthYZ.png", outputDir.c_str( ) ) );

    pHistEnergyLengthXY->SetXTitle( "Energy [keV]" );
    pHistEnergyLengthXY->SetYTitle( "Length [cm]" );
    pHistEnergyLengthXY->SetZTitle( "Events" );
    pHistEnergyLengthXY->Draw( "colz" );
    g_SRIM_F->Draw( "same" );
    // g_SRIM_He->Draw( "same" );
    // g_SRIM_H->Draw( "same" );
    cvs.SaveAs( Form( "%s/energyLengthXY.png", outputDir.c_str( ) ) );

    pHistEnergyLengthXZCut->SetXTitle( "Energy [keV]" );
    pHistEnergyLengthXZCut->SetYTitle( "Length [cm]" );
    pHistEnergyLengthXZCut->SetZTitle( "Events" );
    pHistEnergyLengthXZCut->Draw( "colz" );
    g_SRIM_F->Draw( "same" );
    // g_SRIM_He->Draw( "same" );
    // g_SRIM_H->Draw( "same" );
    cvs.SaveAs( Form( "%s/energyLengthXZCut.png", outputDir.c_str( ) ) );

    pHistEnergyLengthYZCut->SetXTitle( "Energy [keV]" );
    pHistEnergyLengthYZCut->SetYTitle( "Length [cm]" );
    pHistEnergyLengthYZCut->SetZTitle( "Events" );
    pHistEnergyLengthYZCut->Draw( "colz" );
    g_SRIM_F->Draw( "same" );
    // g_SRIM_He->Draw( "same" );
    // g_SRIM_H->Draw( "same" );
    cvs.SaveAs( Form( "%s/energyLengthYZCut.png", outputDir.c_str( ) ) );

    pHistEnergyToTCut->SetXTitle( "Energy [keV]" );
    pHistEnergyToTCut->SetYTitle( "#SigmaToT [us]" );
    pHistEnergyToTCut->SetZTitle( "Events" );
    pHistEnergyToTCut->Draw( "colz" );
    cvs.SaveAs( Form( "%s/energyToTCut.png", outputDir.c_str( ) ) );

    pHistEnergyToTEnCut->SetXTitle( "Energy [keV]" );
    pHistEnergyToTEnCut->SetYTitle( "#SigmaToT / Energy [us/keV]" );
    pHistEnergyToTEnCut->SetZTitle( "Events" );
    pHistEnergyToTEnCut->Draw( "colz" );
    cvs.SaveAs( Form( "%s/energyToTEnCut.png", outputDir.c_str( ) ) );


    double eff50Mino  = (double)selEvt50Mino  / (double)selEvt50 ;
    double eff100Mino = (double)selEvt100Mino / (double)selEvt100;
    double eff200Mino = (double)selEvt200Mino / (double)selEvt200;
    double eff400Mino = (double)selEvt400Mino / (double)selEvt400;
    double effAllMino = (double)selEvtAllMino / (double)selEvtAll;

    std::cout << "The number of events (30  < E < 50  keV)\t:\t" << selEvt50  << "\t" << selEvt50Mino  << "\t" << eff50Mino  << "\t" << sqrt( eff50Mino  * (1.0 - eff50Mino)  / (double)selEvt50  ) << std::endl;
    std::cout << "The number of events (50  < E < 100 keV)\t:\t" << selEvt100 << "\t" << selEvt100Mino << "\t" << eff100Mino << "\t" << sqrt( eff100Mino * (1.0 - eff100Mino) / (double)selEvt100 ) << std::endl;
    std::cout << "The number of events (100 < E < 200 keV)\t:\t" << selEvt200 << "\t" << selEvt200Mino << "\t" << eff200Mino << "\t" << sqrt( eff200Mino * (1.0 - eff200Mino) / (double)selEvt200 ) << std::endl;
    std::cout << "The number of events (200 < E < 400 keV)\t:\t" << selEvt400 << "\t" << selEvt400Mino << "\t" << eff400Mino << "\t" << sqrt( eff400Mino * (1.0 - eff400Mino) / (double)selEvt400 ) << std::endl;
    std::cout << "The number of events (30 < E keV)\t:\t"        << selEvtAll << "\t" << selEvtAllMino << "\t" << effAllMino << "\t" << sqrt( effAllMino * (1.0 - effAllMino) / (double)selEvtAll ) << std::endl;


    gPad->SetRightMargin( 0.1 );

    pHistHitmapAllSelXY->SetXTitle( "#it{x} [cm]" );
    pHistHitmapAllSelXY->SetYTitle( "#it{y} [cm]" );
    pHistHitmapAllSelXY->SetZTitle( "Events" );
    pHistHitmapAllSelXY->SetMarkerStyle(24);
    pHistHitmapAllSelXY->SetMarkerColor(kRed);
    pHistHitmapAllSelXY->Draw( "" );
    // pHistHitmapAllSelXY->Draw( "colz" );
    CreateDrawText( 0.15, 0.88, "^{252}Cf run, X-Y projection hit map", 0.05);
    // CreateDrawText( 0.45, 0.88, "X-Y projection hit map", 0.05);
    // CreateDrawText( 0.15, 0.83, "X-Y projection hit map", 0.04);
    TBox fidAreaXY(0.0, -1.28, 1.28, 1.28 );
    fidAreaXY.SetFillStyle(0);
    fidAreaXY.SetLineColor(kBlue);
    fidAreaXY.SetLineWidth(2);
    fidAreaXY.Draw();
    cvs.SaveAs( Form( "%s/hitmapAllSelXY.png", outputDir.c_str( ) ) );
    cvs.SaveAs( Form( "%s/hitmapAllSelXY.pdf", outputDir.c_str( ) ) );

    pHistHitmapAllSelXZ->SetXTitle( "#it{x} [cm]" );
    pHistHitmapAllSelXZ->SetYTitle( "#it{z} [cm]" );
    pHistHitmapAllSelXZ->SetZTitle( "Events" );
    pHistHitmapAllSelXZ->SetMarkerStyle(24);
    pHistHitmapAllSelXZ->SetMarkerColor(kRed);
    pHistHitmapAllSelXZ->Draw( "" );
    // pHistHitmapAllSelXZ->Draw( "colz" );
    CreateDrawText( 0.15, 0.88, "^{252}Cf run, X-Z projection hit map", 0.05);
    // CreateDrawText( 0.45, 0.88, "X-Z projection hit map", 0.05);
    // CreateDrawText( 0.15, 0.83, "X-Z projection hit map", 0.04);
    TBox fidAreaXZ(0.0, 0.0, 1.28, 14.4 );
    fidAreaXZ.SetFillStyle(0);
    fidAreaXZ.SetLineColor(kBlue);
    fidAreaXZ.SetLineWidth(2);
    fidAreaXZ.Draw();
    cvs.SaveAs( Form( "%s/hitmapAllSelXZ.png", outputDir.c_str( ) ) );
    cvs.SaveAs( Form( "%s/hitmapAllSelXZ.pdf", outputDir.c_str( ) ) );

    pHistHitmapAllSelYZ->SetXTitle( "#it{y} [cm]" );
    pHistHitmapAllSelYZ->SetYTitle( "#it{z} [cm]" );
    pHistHitmapAllSelYZ->SetZTitle( "Events" );
    pHistHitmapAllSelYZ->SetMarkerStyle(24);
    pHistHitmapAllSelYZ->SetMarkerColor(kRed);
    pHistHitmapAllSelYZ->Draw( "" );
    // pHistHitmapAllSelYZ->Draw( "colz" );
    CreateDrawText( 0.15, 0.88, "^{252}Cf run, Y-Z projection hit map", 0.05);
    // CreateDrawText( 0.45, 0.88, "Y-Z projection hit map", 0.05);
    // CreateDrawText( 0.15, 0.83, "Y-Z projection hit map", 0.04);
    TBox fidAreaYZ(-1.28, 0.0, 1.28, 14.4 );
    fidAreaYZ.SetFillStyle(0);
    fidAreaYZ.SetLineColor(kBlue);
    fidAreaYZ.SetLineWidth(2);
    fidAreaYZ.Draw();
    cvs.SaveAs( Form( "%s/hitmapAllSelYZ.png", outputDir.c_str( ) ) );
    cvs.SaveAs( Form( "%s/hitmapAllSelYZ.pdf", outputDir.c_str( ) ) );


    pHistRiseMeanSigma->SetXTitle( "Mean risetime [#mus]" );
    pHistRiseMeanSigma->SetYTitle( "Sigma risetime [#mus]" );
    pHistRiseMeanSigma->Draw("colz");
    cvs.SaveAs( Form( "%s/rise.png", outputDir.c_str( ) ) );
    cvs.SaveAs( Form( "%s/rise.pdf", outputDir.c_str( ) ) );

    pHistRiseSigmaY->SetXTitle( "Y [cm]" );
    pHistRiseSigmaY->SetYTitle( "Sigma risetime [#mus]" );
    pHistRiseSigmaY->Draw("colz");
    cvs.SaveAs( Form( "%s/riseSigmaY.png", outputDir.c_str( ) ) );
    cvs.SaveAs( Form( "%s/riseSigmaY.pdf", outputDir.c_str( ) ) );

    TCanvas cvs2( "cvs2", "cvs2", 800, 600 );
    const double eneBinning[6] = { 0.0, 30.0, 50.0, 100.0, 200.0, 400.0 };
    TH1F histNum( "histMinoEffNum", "histMinoEffNum", 5, eneBinning );
    TH1F histDen( "histMinoEffDen", "histMinoEffDen", 5, eneBinning );
    
    for(int i = 0; i < selEvt50Mino;  ++i ) histNum.Fill( 40.0  );
    for(int i = 0; i < selEvt100Mino; ++i ) histNum.Fill( 75.0  );
    for(int i = 0; i < selEvt200Mino; ++i ) histNum.Fill( 150.0 );
    for(int i = 0; i < selEvt400Mino; ++i ) histNum.Fill( 300.0 );

    for(int i = 0; i < selEvt50;  ++i ) histDen.Fill( 40.0  );
    for(int i = 0; i < selEvt100; ++i ) histDen.Fill( 75.0  );
    for(int i = 0; i < selEvt200; ++i ) histDen.Fill( 150.0 );
    for(int i = 0; i < selEvt400; ++i ) histDen.Fill( 300.0 );

    histNum.Draw( "histe" );
    cvs2.SaveAs( Form( "%s/eneEffNum.png", outputDir.c_str( ) ) );
    histDen.Draw( "histe" );
    cvs2.SaveAs( Form( "%s/eneEffDen.png", outputDir.c_str( ) ) );

    TEfficiency eff( histNum, histDen );
    eff.SetStatisticOption( TEfficiency::kFCP );
    eff.SetMarkerStyle( 20 );
    eff.SetMarkerColor( kRed );
    eff.SetLineWidth( 2.0 );
    eff.SetLineColor( kRed );
    
    TH1* pHistFrame = cvs2.DrawFrame( 0.0, 0.0, 400.0, 1.5 );
    pHistFrame->GetXaxis()->SetTitle( "Energy [keV]" );
    pHistFrame->GetXaxis()->SetTitleSize( 0.05 );
    pHistFrame->GetXaxis()->SetTitleOffset( 1.0 );
    pHistFrame->GetXaxis()->SetLabelSize( 0.05 );
    pHistFrame->GetYaxis()->SetTitle( "2-peak detection efficiency" );
    pHistFrame->GetYaxis()->SetTitleSize( 0.05 );
    pHistFrame->GetYaxis()->SetTitleOffset( 1.0 );
    pHistFrame->GetYaxis()->SetLabelSize( 0.05 );
    gPad->Update( );
    eff.Draw( "same" );

    CreateDrawText( 0.16, 0.85, "^{252}Cf run", 0.07 );
    CreateDrawText( 0.16, 0.75, Form( "Total: %0.0lf #pm %0.0lf %%", effAllMino * 100.0, sqrt( effAllMino * (1.0 - effAllMino) / (double)selEvtAll ) * 100.0 ), 0.07 );

    TLine effLine( 0.0, 1.0, 400.0, 1.0 );
    effLine.SetLineWidth( 2 );
    effLine.SetLineColor( kBlack );
    effLine.SetLineStyle( 3 );
    effLine.Draw( );

    cvs2.SaveAs( Form( "%s/eneEff.png", outputDir.c_str( ) ) );
    cvs2.SaveAs( Form( "%s/eneEff.pdf", outputDir.c_str( ) ) );


    pHistAllDt->GetXaxis()->SetTitle( "#Deltat [#mus]" );
    pHistAllDt->GetYaxis()->SetTitle( "#Strips" );
    pHistAllDt->Draw( "" );
    cvs2.SaveAs( Form( "%s/allDt.png", outputDir.c_str( ) ) );
    cvs2.SaveAs( Form( "%s/allDt.pdf", outputDir.c_str( ) ) );


    TCanvas cvs3D("cvs3D", "cvs3D", 800, 800 );
    cvs3D.SetGridx( 1 ); cvs3D.SetGridy( 1 );
    pHistHitmapXYZ->SetFillColor( kRed );
    pHistHitmapXYZ->GetXaxis()->SetRangeUser( -1.5, 1.5 );
    pHistHitmapXYZ->GetXaxis()->SetTitle( "#it{x} [cm]" );
    pHistHitmapXYZ->GetXaxis()->SetTitleOffset( 2.0 );
    // pHistHitmapXYZ->GetXaxis()->SetLavelSize( 0.9 );
    pHistHitmapXYZ->GetXaxis()->SetRangeUser( -0.1, 1.2 );
    pHistHitmapXYZ->GetYaxis()->SetTitle( "#it{y} [cm]" );
    pHistHitmapXYZ->GetYaxis()->SetTitleOffset( 2.5 );
    pHistHitmapXYZ->GetYaxis()->SetLabelOffset( 0.012 );
    pHistHitmapXYZ->GetZaxis()->SetTitle( "#it{z} [cm]" );
    pHistHitmapXYZ->GetZaxis()->SetTitleOffset( 1.2 );
    // pHistHitmapXYZ->GetZaxis()->SetLavelSize( 0.9 );
    pHistHitmapXYZ->Draw( "box" );
    pHistHitmapXYZ->GetXaxis()->SetRangeUser( -0.75, 2.25 );
    pHistHitmapXYZ->SetMarkerStyle( 20 );
    pHistHitmapXYZ->SetMarkerColorAlpha( kRed, 0.65 );
    pHistHitmapXYZ->Draw( "" );

    CreateDrawText( 0.05, 0.94, "^{252}Cf run", 0.05);
    CreateDrawText( 0.05, 0.88, "3D hit map", 0.05);

    cvs3D.SaveAs( Form( "%s/hitmapXYZ.png", outputDir.c_str( ) ) );
    cvs3D.SaveAs( Form( "%s/hitmapXYZ.pdf", outputDir.c_str( ) ) );


    return;
}
