// rootlogon for NIuTPC analysis
{
    std::string decoderPath = gSystem->Getenv("NI_DECODER_DIR");
    std::string includePath = decoderPath + "/source/include";
    std::string includePath = decoderPath + "/source/src";
    gInterpreter->AddIncludePath( includePath.c_str( ) );
    std::cout << "" << gSystem->GetIncludePath( ) << std::endl;
}
