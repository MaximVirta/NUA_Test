#include <TMath.h>
double e1;
double e2;

void binErrors(UInt_t nevt1 = 5000, UInt_t nevt2 = 50000) {
	TFile *fin1 = new TFile(Form("outputs/output_%u_evt.root", nevt1), "READ");
	TFile *fin2 = new TFile(Form("outputs/output_%u_evt.root", nevt2), "READ"); 

	TH1D *h1 = (TH1D*) fin1->Get("hPhiCorrected");
	TH1D *h2 = (TH1D*) fin2->Get("hPhiCorrected");

	for (UInt_t i = 1; i<h1->GetNbinsX()+1; i++) {
		e1 = h1->GetBinError(i);
		e2 = h2->GetBinError(i);

		// printf("%.3f\n", e1/e2*e1/e2);
	}
	e1 = h1->GetMeanError();
	e2 = h2->GetMeanError();
	printf("Increasing statistics from %d to %d, reduces error by: %.3f\nSqrt(nevt2/nevt1)=%.3f\n", nevt1, nevt2, e1/e2, TMath::Sqrt(1.0*nevt2/nevt1));
	fin1->Close();
	fin2->Close();
}