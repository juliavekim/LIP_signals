void TestGaussian() {

   using namespace RooFit;

   TCanvas *C = new TCanvas( "c", "c",  800, 800);
   //C->SetLogx();
   //C->SetLogy();


   // gaus sig
   TH1F h1("h1", "histogram from a gaussian", 100, 0, 10);
   //h1.FillRandom("gaus", 10000);
   TRandom3 rndgen;
   for (int i = 0; i < 100000; i++) h1.Fill(rndgen.Gaus(5, 1));
    
   // constant bkg
   TH1F h2("h2", "histogram from uniform", 1000, 0, 10);
   h2.FillRandom("pol0", 100000);


   RooRealVar E("E", "E", 0, 10);
   RooDataHist hgaus("hgaus", "hgaus", E, Import(h1));
   RooDataHist hfbg("hfbg", "hfbg", E, Import(h2));

   RooHistPdf bkg("bg_pdf", "bg_pdf", E, hfbg, 1);
   RooHistPdf sig("sig_pdf", "sig_pdf", E, hgaus, 1);
  
   // plotting PDFs
   //RooPlot* frame1 = E.frame(Title("Gaussian Signal & Flat Background PDFs"));
   //sig.plotOn(frame1, LineColor(kRed));
   //bkg.plotOn(frame1);
   //frame1->Draw();

    // build model with initial singal/background ratio 20:80
    RooRealVar norm_s("norm_s","N_{s}",0,-1000,1000);
    RooRealVar norm_b("norm_b","N_{b}",10000,0,20000);
    
    // model 1: (sig+bkg) for generating dataset
    const RooArgList components(sig, bkg);
    const RooArgList coeffs(norm_s, norm_b);    
    RooAddPdf model("model","f_{s+b}",components,coeffs);


    Int_t N_seeds = 10000;
    Int_t Starting_seed = 0;


   // open a new root file to save the results.
   TFile *f = new TFile("TestGaussian.root", "recreate");

   TH1F norm_s_hat("norm_s_hat", "norm_s_hat", 100, -500, 500);
   TH1F norm_b_hat("norm_b_hat", "norm_b_hat", 100, 9500, 10500);
   for (int seed=Starting_seed; seed< Starting_seed+N_seeds; seed++) {
     cout << "Processing seed number " << seed << endl;
     RooRandom::randomGenerator()->SetSeed(seed);

     norm_s.setVal(0.);
     norm_b.setVal(10000.);

     // Getting our total number of events from Poisson
     //Int_t t_events = random.Poisson(norm_b.getVal());
     
     // MC Toys
     RooDataSet *newdata = model.generate(E);
     
     // create a likelihood function.
     RooAbsReal* nll = model.createNLL(*newdata, Extended(), NumCPU(8));
     RooMinuit(*nll).migrad();

     cout<< "sum: " << norm_s.getVal() + norm_b.getVal() << endl;
     //norm_s_hat.Fill(norm_s.getVal());
     //norm_b_hat.Fill(norm_b.getVal());
     delete newdata;
     delete nll;
  }
  
  f->Write();
  f->Close();
 
}
                 
