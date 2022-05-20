#include "NAPStyle.h"
#include "NIConfig.hh"
#include "NIConfig.cc"
#include "NAUtil.h"
#include "NAUtil.cc"

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
    TTree* pTree = dynamic_cast< TTree* >( file.Get( "anal_tree" ) );
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

    double calFactor = pConfig->cal_factor;

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

    TH1D* pHistEnergy = new TH1D( "histEnergy", "histEnergy", 100, 0.0, 400.0 );
    TH1D* pHistEnergyCut = new TH1D( "histEnergyCut", "histEnergyCut", 100, 0.0, 400.0 );

    TH2D* pHistHitmapXY = new TH2D( "histHitmapXY", "histHitmapXY", 100, -1.5, 1.5, 100, -1.5,  1.5 );
    TH2D* pHistHitmapXZ = new TH2D( "histHitmapXZ", "histHitmapXZ", 100, -1.5, 1.5, 100,  0.0, 12.0 );
    TH2D* pHistHitmapYZ = new TH2D( "histHitmapYZ", "histHitmapYZ", 100, -1.5, 1.5, 100,  0.0, 12.0 );

    TH2D* pHistHitmapXMaxMin = new TH2D( "histHitmapXMaxMin", "histHitmapXMaxMin", 100, -1.5, 1.5, 100, -1.5,  1.5 );
    TH2D* pHistHitmapYMaxMin = new TH2D( "histHitmapYMaxMin", "histHitmapYMaxMin", 100, -1.5, 1.5, 100, -1.5,  1.5 );

    TH2D* pHistEnergyToT    = new TH2D( "histEnergyToT", "histEnergyToT", 100, 0.0, 400.0, 100, 0.0, 2000.0 );
    TH2D* pHistEnergyToTEn  = new TH2D( "histEnergyToTEn", "histEnergyToTEn", 100, 0.0, 400.0, 100, 0.0, 10.0 );
    TH2D* pHistEnergyLengthXZ = new TH2D( "histEnergyLengthXZ", "histEnergyLengthXZ", 100, 0.0, 400.0, 100, 0.0, 2.0 );
    TH2D* pHistEnergyLengthYZ = new TH2D( "histEnergyLengthYZ", "histEnergyLengthYZ", 100, 0.0, 400.0, 100, 0.0, 2.0 );
    TH2D* pHistEnergyLengthXY = new TH2D( "histEnergyLengthXY", "histEnergyLengthXY", 100, 0.0, 400.0, 100, 0.0, 2.0 );

    TH2D* pHistHitmapXY_cut = new TH2D( "histHitmapXY_cut", "histHitmapXY_cut", 100, -1.5, 1.5, 100, -1.5,  1.5 );

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
    for( int evt = 0; evt < totEvt; ++evt ) {
        NAUtil::PrintProgressBar( evt, totEvt );
        pTree->GetEntry( evt );
        
        if( pVec_yz_y->size( ) <= 0 ) continue;

        pHistHitmapXY->Fill( ave_x, ave_y );
        pHistHitmapXZ->Fill( ave_x, ave_z );
        pHistHitmapYZ->Fill( ave_y, ave_z );

        double sumToT = 0.0;
        int nHit = pVec_a_hg_mainrise_time->size( );
        if( pVec_a_hg_mainfall_time->size( ) == nHit ) {
            for( int i = 0; i < nHit; ++i ) {
                sumToT += ( pVec_a_hg_mainfall_time->at( i ) - pVec_a_hg_mainrise_time->at( i ) );
            }
            pHistEnergyToT->Fill( a_hg_sum_charge * calFactor, sumToT );
            pHistEnergyToTEn->Fill( a_hg_sum_charge * calFactor, sumToT/(a_hg_sum_charge * calFactor) );
        }
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
        // if( pVec_yz_y->size( ) <= 0 ) {
        //     lengthXZ = 0.0;
        //     lengthYZ = 0.0;
        // }

        // if( lengthXZ > 2.0 ) pHistHitmapXY_cut->Fill( ave_x, ave_y );

        pHistHitmapXMaxMin->Fill( xz_x_max, xz_x_min );
        pHistHitmapYMaxMin->Fill( yz_y_max, yz_y_min );

        // if( fabs( ave_x ) < 0.8 && fabs( ave_y ) < 0.8 && a_hg_sum_charge * calFactor > 10.0 &&
        //     yz_y_min > -1.0 && yz_y_max < 1.0 && xz_x_min > -1.0 && xz_x_max < 1.0 ) {
        // if( fabs( ave_x ) < 0.8 && fabs( ave_y ) < 0.8 && a_hg_sum_charge * calFactor > 10.0 && xz_x_max < 1.0 && xz_x_min > -1.0 && yz_y_min > -1.0 && yz_y_max < 1.0 ) {
        // if( fabs( ave_x ) < 0.8 && fabs( ave_y ) < 0.8 && a_hg_sum_charge * calFactor > 100.0 && xz_x_max > 0.8 && xz_x_min < -0.8 ) {
        // if( lengthXZ < 1.5 && a_hg_sum_charge * calFactor > 100.0 )
        if( xz_x_max > 1 && xz_x_min < 0 )
            std::cout << Form("%d\t%lf\t%lf\t%lf\t%lf\t%lf",ev, a_hg_sum_charge * calFactor, lengthXZ, lengthXY, lengthX, lengthY) << std::endl;

        pHistEnergyLengthXZ->Fill( a_hg_sum_charge * calFactor, lengthXZ );
        pHistEnergyLengthXY->Fill( a_hg_sum_charge * calFactor, lengthXY );
        pHistEnergyLengthYZ->Fill( a_hg_sum_charge * calFactor, lengthYZ );

        pHistEnergy->Fill( a_hg_sum_charge * calFactor );
        // if( lengthXZ > 2.0 || lengthYZ > 2.0 ) pHistEnergyCut->Fill( a_hg_sum_charge * calFactor );
        // if( ( lengthXZ > 1.0 || lengthYZ > 1.0 ) && ( lengthXZ < 1.3 || lengthYZ < 1.3 ) ) pHistEnergyCut->Fill( a_hg_sum_charge * calFactor );
        // std::cout << "xmin= " << xz_x_min << ", xmax= " << xz_x_max << ", ymin= " << yz_y_min << ", ymax= " << yz_y_max << ", zmin= " << xz_z_min << ", zmax= " << xz_z_max << std::endl;
    }
    
    TCanvas cvs( "cvs", "cvs", 800, 800 );

    gPad->SetRightMargin( 0.1 );

    pHistEnergy->GetXaxis()->SetTitle( "#SigmaADC count" );
    pHistEnergy->GetYaxis()->SetTitle( "Events" );
    pHistEnergyCut->SetLineColor( kRed );

    pHistEnergy->Draw( "" );
    pHistEnergyCut->Draw( "same" );

    TF1 fitGauss( "fitGauss", "gaus", 50, 150 );
    pHistEnergyCut->Fit(&fitGauss,"LMI","",50, 150);
    double peakVal = fitGauss.GetParameter( 1 );
    CreateDrawText( 0.48, 0.85, Form( "Default cal. factor = %0.6lf", calFactor ), 0.035 );
    CreateDrawText( 0.48, 0.8, Form( "New cal. factor = %0.6lf", calFactor * SIM_ENERGY_DEPOSIT / peakVal ), 0.035 );
    cvs.SaveAs( Form( "%s/energy.png", outputDir.c_str( ) ) );

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

    pHistHitmapXY_cut->SetXTitle( "X [cm]" );
    pHistHitmapXY_cut->SetYTitle( "Y [cm]" );
    pHistHitmapXY_cut->SetZTitle( "Events" );
    pHistHitmapXY_cut->Draw( "colz" );
    cvs.SaveAs( Form( "%s/hitmapXY_cut.png", outputDir.c_str( ) ) );

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

    pHistEnergyLengthXZ->SetXTitle( "Energy [keV]" );
    pHistEnergyLengthXZ->SetYTitle( "Length [cm]" );
    pHistEnergyLengthXZ->SetZTitle( "Events" );
    pHistEnergyLengthXZ->Draw( "colz" );
    g_SRIM_F->Draw( "same" );
    g_SRIM_He->Draw( "same" );
    g_SRIM_H->Draw( "same" );
    cvs.SaveAs( Form( "%s/energyLengthXZ.png", outputDir.c_str( ) ) );

    pHistEnergyLengthYZ->SetXTitle( "Energy [keV]" );
    pHistEnergyLengthYZ->SetYTitle( "Length [cm]" );
    pHistEnergyLengthYZ->SetZTitle( "Events" );
    pHistEnergyLengthYZ->Draw( "colz" );
    g_SRIM_F->Draw( "same" );
    g_SRIM_He->Draw( "same" );
    g_SRIM_H->Draw( "same" );
    cvs.SaveAs( Form( "%s/energyLengthYZ.png", outputDir.c_str( ) ) );

    pHistEnergyLengthXY->SetXTitle( "Energy [keV]" );
    pHistEnergyLengthXY->SetYTitle( "Length [cm]" );
    pHistEnergyLengthXY->SetZTitle( "Events" );
    pHistEnergyLengthXY->Draw( "colz" );
    g_SRIM_F->Draw( "same" );
    g_SRIM_He->Draw( "same" );
    g_SRIM_H->Draw( "same" );
    cvs.SaveAs( Form( "%s/energyLengthXY.png", outputDir.c_str( ) ) );

    return;
}
