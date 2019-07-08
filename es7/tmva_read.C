int TMVARead(){

    TMVA::Tools::Instance();

    //Tree Reading

    TFile* dataset = TFile::Open("dataset.root");

    TTree* signal = (TTree*) dataset->Get("Signal");
    TTree* background = (TTree*) dataset->Get("Background");

    TFile inputClass("TMVATrain.root");

    TDirectoryFile* NN_classification = (TDirectoryFile*) inputClass.Get("NN_classification");

    TTree* trainTree = (TTree*)NN_classification->Get("TrainTree");
    TTree* testTree = (TTree*)NN_classification->Get("TestTree");

    TMVA::Factory* factoru = (TMVA::Factory*) NN_classification->Get("factory");
    TMVA::DataLoader* dataloader = (TMVA::DataLoader*)NN_classification->Get("dataloader");

    TString outputfilename = "TMVARead.root";

    TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");


    //Defining variables
    float x,y;


    reader->AddVariable("x", &x);
    reader->AddVariable("y", &y);

    reader->BookMVA("NN", "./NN_classification/weights/TMVAClassification_MLP.weights.xml");

    signal->SetBranchAddress("x", &x);
    signal->SetBranchAddress("y", &y);

    background->SetBranchAddress("x", &x);
    background->SetBranchAddress("y", &y);



return 0;
}