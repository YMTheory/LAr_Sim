

void ellipse_fit()
{
    TFile* ff = new TFile("../outputs/coutour_Ba.root", "read");
    TCanvas* cc = (TCanvas*)ff->Get("c1");
    TGraph* best   = (TGraph*)cc->GetPrimitive("best");
    TGraph* sigma1 = (TGraph*)cc->GetPrimitive("sigma1");
    TGraph* sigma5 = (TGraph*)cc->GetPrimitive("sigma5");

}
