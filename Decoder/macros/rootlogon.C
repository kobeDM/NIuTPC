// rootlogon for NIuTPC analysis
{
    std::string decoderPath = gSystem->Getenv("NI_DECODER_DIR");
    std::string incPath     = decoderPath + "/source/include";
    std::string srcPath     = decoderPath + "/source/src";
    gInterpreter->AddIncludePath( incPath.c_str( ) );
    gInterpreter->AddIncludePath( srcPath.c_str( ) );
    std::cout << "" << gSystem->GetIncludePath( ) << std::endl;
}
