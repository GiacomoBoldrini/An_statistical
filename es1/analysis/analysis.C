void analysis(){

    std::fstream uni_genFile;
    uni_genFile.open("../data/mc_hm_bw.txt");

    TCanvas * c1 = new TCanvas("uni_gen","uni_gen", 700,700);

    double min = -10;
    double max = 10.;
    int bins = 100;

    TH1D * uni_Histo = new TH1D("uni_gen","uni_gen", bins,min,max);
    
    double entry;

    while(uni_genFile >> entry){
        uni_Histo -> Fill(entry);
    }

    uni_Histo -> Draw();

    uni_Histo->SetFillColorAlpha(kBlue, 0.7);
    c1->Print("mc_hw_bw.pdf","pdf");
}