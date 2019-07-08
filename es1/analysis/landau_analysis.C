double landau_fun(double x, double mu, double sigma){

    return TMath::Landau(x,mu,sigma);
}

void line(){

    std::cout << "\n ======================================================= \n" << std::endl;
}


 void landau_analysis(string file_name = "../data/landau.txt"){

     //string file_name = "../data/landau.txt"; //data under study
     string save_path = "../analysis/landau.pdf"; //canvas print path

     std::cout << "File: " << file_name << std::endl;

     std::fstream landau_file;
     landau_file.open(file_name);

     //Histograms parameters
     double min = 0;
     double max = 10;
     double bins = 100;

     //Canvas
     TCanvas* c1 = new TCanvas("BW","BW",700,700);

     //Histo
     TH1D* landau_histo = new TH1D("landau_histo","landau_histo",bins,min,max);

     double entry;

     while(landau_file >> entry){

         landau_histo -> Fill(entry);
     }

     double norm = landau_histo->GetEntries()*(max-min)/bins;
     double over = landau_histo->GetBin(bins+1);
     double under = landau_histo->GetBin(0);

     //std::cout << "Norm: " << norm << std::endl;
     //std::cout << "Overflow: " << over << std::endl;
     //std::cout << "Underflow: " << under << std::endl;

     //TF1
     TF1 * landau_pdf = new TF1("landauPDF","[2]*landau_fun(x,[0],[1])",min,max);

     landau_pdf->SetParameters(2,0.5);
     landau_pdf->SetParName(0,"#mu");
     landau_pdf->SetParName(1,"#sigma");
     //bw_pdf->SetParName(2,"Norm");

     gStyle->SetOptStat(1110011);
     gStyle->SetOptFit(1111);

     landau_histo->Fit("landauPDF");

     landau_histo->Draw();   
          
     c1->Print(save_path.c_str(),"pdf"); 

     //Hip test
     double chi_sq = landau_pdf->GetChisquare();
     double dof = landau_pdf->GetNDF();
     double p = landau_pdf->GetProb();
     
     line();
     std::cout << "Chi Square tilde: " << chi_sq / dof << " Probability: " << p << std::endl;
     line();   

         
 }
