
int TMVATrain(){

    TMVA::Tools::Instance();

    //TTree reading
    
    TFile * dataset = TFile::Open("dataset.root");

    TTree* signal = (TTree*) dataset ->Get("Signal");
    TTree* background = (TTree*) dataset ->Get("Background");

    TString outputfileName = ("TMVATrain.root");
    TFile* outputfile = new TFile(outputfileName, "RECREATE");

    //Factory definition

    TMVA::Factory* factory = new TMVA::Factory
    (
        "TMVAClassification",
        outputfile,
        "!V:!Silent:Color:DrawProgressBar:Transformations=I;P;G:AnalysisType=Classification"

    );

    //Adding input variables

    TMVA::DataLoader * dataloader = new TMVA::DataLoader("NN_classification");

    dataloader->AddVariable("x",'F');
    dataloader->AddVariable("y",'F');

    dataloader->AddSignalTree(signal, 1.); //tree, weights
    dataloader->AddBackgroundTree(background, 1.);

    TString sample_train = "8000";
    TCut cut_sig = "";
    TCut cut_bkg = cut_sig;


    dataloader->PrepareTrainingAndTestTree
    (
        cut_sig,
        cut_bkg,
        "nTrain_signal=8000:nTrain_Background=8000:SplitMode=Random:NormMode=EqualNumEvents:!V"
    );
 

      factory->BookMethod 
    (
      dataloader, 
      TMVA::Types::kMLP, 
      "MLP", 
      "!H:!V:NeuronType=tanh:VarTransform=Gauss:NCycles=200:HiddenLayers=5:TestRate=15:TrainingMethod=BP:SamplingTesting=True:ConvergenceImprove=1e-5:ConvergenceTests=75" 
    );

      factory->BookMethod
    ( 
       dataloader,
       TMVA::Types::kSVM, 
       "SVM"
    );
    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();   
    
  factory->Write();
  //ROC->Write();
  dataloader->Write();
  
  outputfile->Close () ;

  delete factory ;
  delete dataloader ;
  delete signal;
  delete background ;
  delete outputfile ;

  if (!gROOT->IsBatch()) TMVA::TMVAGui( outputfileName );
    
    return 0; 
    
}

