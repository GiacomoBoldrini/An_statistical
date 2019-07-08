double landau_fun(double x, double mu, double sigma){

    return TMath::Landau(x,mu,sigma);
}

void line(){

    std::cout << "\n ======================================================= \n" << std::endl;
}


 void analysis(string file_name = "prova.txt"){

     //string file_name = "../data/landau.txt"; //data under study
     string save_path = ""; //canvas print path

     std::cout << "File: " << file_name << std::endl;

     std::fstream landau_file;
     landau_file.open(file_name);

     //Histograms parameters
     double min = 0;
     double max = 100;
     double bins = 1000;

     //Canvas
     TCanvas* c1 = new TCanvas("BW","BW",1000,800);

     //Histo
     TH1D* landau_histo = new TH1D("diff pdf","diff pdf",bins,min,max);

     double entry;

     while(landau_file >> entry){
         std::cout << entry << std::endl;
         landau_histo -> Fill(entry);
     }  

     landau_histo->Draw();   
          
     c1->Print("prova.root","root"); 

     //Hip test
     /*/double chi_sq = landau_pdf->GetChisquare();
     double dof = landau_pdf->GetNDF();
     double p = landau_pdf->GetProb();
     
     line();
     std::cout << "Chi Square tilde: " << chi_sq / dof << " Probability: " << p << std::endl;
     line();*/

         
 }
