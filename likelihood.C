{
   using namespace RooFit;

   TCanvas *C = new TCanvas( "c", "c",  800, 800);
   //C->SetLogx();
   //C->SetLogy();
   
   // open a new root file to save the results.
   TFile *f = new TFile("limit.root", "recreate");
   TH1F *limit_hist = new TH1F("limit_hist","limit_hist",1000, 0, 1000); 

   //TFile *f = new TFile("limit.root", "update");
   //TH1F *limit_hist = (TH1F*)f->Get("limit_hist");  

   // set up bin edges for histogram
   // same bin edges used in CUTE bgexplorer 
  const Int_t NBINS = 3100;
   Double_t edges[NBINS + 1] = {};
   for (int i = 0; i < 10; i++) {
       edges[i] = 0.01*i;
   }
   for (int i = 10; i < 110; i++) {
       edges[i] = 0.1 + 0.1*(i-10);
   }
   for (int i = 110; i < 3102; i++) {
       edges[i] = 11 + i-110;
   }
   


   TH1 *bg_hist = new TH1D("bg_PDF", "bg_PDF", NBINS, edges);
   TH1 *sig_hist = new TH1D("sig_PDF", "sig_PDF", NBINS, edges);

   // read in bkg events
   fstream file0;
   file0.open("mock_G124Singles_bg.txt", ios::in);

   double bg_value;
   while(1)
   {
       file0 >> bg_value;
       bg_hist->Fill(bg_value);
       if(file0.eof()) break;
   }
   file0.close();

   // read in sig events
   fstream file1;
   file1.open("mock_LIP.txt", ios::in);

   double sig_value;
   while(1)
   {
       file1 >> sig_value;
       sig_hist->Fill(sig_value);
       if(file1.eof()) break;
   }
   file1.close();

   // declaring observables and choose analysis window
   RooRealVar E("E", "E", 0.1, 10);
   RooDataHist dh0("dh0", "dh0", E, Import(*bg_hist));
   RooHistPdf bkg("bg_pdf", "bg_pdf", E, dh0, 1);
   RooDataHist dh1("dh1", "dh1", E, Import(*sig_hist));
   RooHistPdf sig("sig_pdf", "sig_pdf", E, dh1, 1);

   //RooPlot* frame1 = E.frame(Title("bkg pdf"));
   //bkg.plotOn(frame1);
   //frame1->Draw();
   //RooPlot* frame1 = E.frame(Title("sig pdf"));
   //sig.plotOn(frame1);
   //frame1->Draw();


    // build model
    RooRealVar norm_s("norm_s","N_{s}",0,0,10000);
    RooRealVar norm_b("norm_b","N_{b}",1,0,10000);
    const RooArgList components(sig, bkg);
    const RooArgList coeffs(norm_s, norm_b);
    RooAddPdf model("model","f_{s+b}",components,coeffs);

    // generate a toy MC sample from model
    Int_t N_seeds = 50;
    Int_t Starting_seed = 0;
    // loop through different random seeds to generate multiple datasets.
    for (int seed=Starting_seed; seed< Starting_seed+N_seeds; seed++) {
      cout << "Processing seed number " << seed << endl;
      // reset norm_s and norm_b values
      norm_s.setVal(0);
      norm_b.setVal(1);

      RooRandom::randomGenerator()->SetSeed(seed);
      RooDataSet *newdata = model.generate(E, 10000);

      // create a likelihood function.
      RooAbsReal* nll = model.createNLL(*newdata, NumCPU(8));
      // minimize the likelihood.
      RooMinuit(*nll).migrad();
      // plot likelihood scan
      //RooPlot *frame = norm_s.frame();
      //nll->plotOn(frame, ShiftToZero());
      // plot profile likelihood scan
      RooAbsReal* pll_norm_s = nll->createProfile(norm_s);
      //pll_norm_s->plotOn(frame, LineColor(kRed));    
      //frame->Draw();
   
      TF1 *pll = pll_norm_s -> asTF(norm_s);
      // loop through norm_s to find the 95% CL point
      // according to chi-squared distribution
    
      Double_t q_test = 2* (pll->Eval(1));
      Double_t q_test_last = 2* (pll->Eval(0)); 
      for (int i= 1; i < 10000; i++) {
          if ((q_test > 1.64*1.64) && (q_test_last < 1.64*1.64)) {
              cout << "q_test: " << q_test << endl;
              cout << "limit: " << i << endl; 
              limit_hist->Fill(i);
              break;
          }
          q_test_last = q_test;
          q_test = 2* (pll->Eval(i));
      }
   }   
   
   //limit_hist->Draw();
   f->Write();
   f->Close();   
}
