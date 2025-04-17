#ifndef NAPStyle_hh
#define NAPStyle_hh


#include <TROOT.h>

TStyle* NAPStyle( )
{
    TStyle *pNAPStyle = new TStyle("NAP","NAP style");

    // use plain black on white colors
    Int_t icol=0; // WHITE
    pNAPStyle->SetFrameBorderMode(icol);
    pNAPStyle->SetFrameFillColor(icol);
    pNAPStyle->SetCanvasBorderMode(icol);
    pNAPStyle->SetCanvasColor(icol);
    pNAPStyle->SetPadBorderMode(icol);
    pNAPStyle->SetPadColor(icol);
    pNAPStyle->SetStatColor(icol);
    //pNAPStyle->SetFillColor(icol); // don't use: white fill color for *all* objects

    // set the paper & margin sizes
    pNAPStyle->SetPaperSize(20,26);

    // set margin sizes
    pNAPStyle->SetPadTopMargin(0.05);
    pNAPStyle->SetPadRightMargin(0.05);
    pNAPStyle->SetPadBottomMargin(0.12);
    pNAPStyle->SetPadLeftMargin(0.12);

    // set title offsets (for axis label)
    pNAPStyle->SetTitleXOffset(1.3);
    pNAPStyle->SetTitleYOffset(1.3);

    // use large fonts
    //Int_t font=72; // Helvetica italics
    Int_t font=42; // Helvetica
    Double_t tsize=0.04;
    pNAPStyle->SetTextFont(font);

    pNAPStyle->SetTextSize(tsize);
    pNAPStyle->SetLabelFont(font,"x");
    pNAPStyle->SetTitleFont(font,"x");
    pNAPStyle->SetLabelFont(font,"y");
    pNAPStyle->SetTitleFont(font,"y");
    pNAPStyle->SetLabelFont(font,"z");
    pNAPStyle->SetTitleFont(font,"z");

    pNAPStyle->SetLabelSize(tsize,"x");
    pNAPStyle->SetTitleSize(tsize,"x");
    pNAPStyle->SetLabelSize(tsize,"y");
    pNAPStyle->SetTitleSize(tsize,"y");
    pNAPStyle->SetLabelSize(tsize,"z");
    pNAPStyle->SetTitleSize(tsize,"z");

    // use bold lines and markers
    // pNAPStyle->SetMarkerStyle(20);
    pNAPStyle->SetMarkerSize(1.2);
    pNAPStyle->SetHistLineWidth(1);
    // pNAPStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

    // get rid of X error bars
    //pNAPStyle->SetErrorX(0.001);
    // get rid of error bar caps
    pNAPStyle->SetEndErrorSize(0.);

    // do not display any of the standard histogram decorations
    pNAPStyle->SetOptTitle(0);
    //pNAPStyle->SetOptStat(1111);
    pNAPStyle->SetOptStat(0);
    //pNAPStyle->SetOptFit(1111);
    pNAPStyle->SetOptFit(0);

    // put tick marks on top and RHS of plots
    pNAPStyle->SetPadTickX(1);
    pNAPStyle->SetPadTickY(1);

    return pNAPStyle;
}

void SetNAPStyle( )
{
    static TStyle* pNAPStyle = nullptr;
    std::cout << "\nApplying NAP style settings...\n" << std::endl ;
    if ( pNAPStyle==0 ) pNAPStyle = NAPStyle();
    gROOT->SetStyle("NAP");
    gROOT->ForceStyle();

    return;
}

#endif // NAPStyle_hh
