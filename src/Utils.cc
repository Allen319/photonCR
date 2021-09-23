#include <Utils.h>

#include <TF1.h>
#include <TH1.h>

namespace utils {

  TH1F * SmoothenHisto(TH1 *h1) {
    TF1 * function = new TF1("func", "[0] + [1] * x", 0, 100);
    h1->Fit("func");
    for (unsigned int i = 0; i < h1->GetNbinsX(); i++) {
      h1->SetBinContent(i, function->Eval(h1->GetBinCenter(i)));
    }
    delete function;
    TH1F *histo = (TH1F*) h1->Clone();
    return histo;
  }
  
};

