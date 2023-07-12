#include "TLegend.h"
#include "TH1.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TString.h"
#include <TStopwatch.h>
#include <TComplex.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <TSystem.h> 
#include "src/JBaseEventHeader.h"
#include "src/JBaseTrack.h"
#include "src/JHistos.h"
#include "src/JPDF.h"


const Double_t kEtaMax = 0.8;
bool b2holes = 1; //1=two holes, 0=1 hole
Int_t GetCentralityBin( Double_t dice);
void NormalizeSample(TH1D *hist);

int main(int argc, char **argv)
{
TROOT root("flow","run mc");
if ( argc<3 ) {
	cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	cout<<"+  "<<argv[0]<<" <fRootFile> <Nevt> <random seed>  <bNUE>"<<endl;
	cout << endl << endl;
	exit(1);
}
// Arguments 
char *sFile = argv[1];
Int_t Nevt= atoi(argv[2]);
Int_t random_seed = atoi(argv[3]);
Int_t bNUE = atoi(argv[4]); //if bNUE = false = 0, uses the uniform acceptance for sampling
// Declare variables
// cout<< strCentrality[0]<<endl;


JFlowInput *jinput = new JFlowInput();
jinput->LoadAliceData();
JPDF *jpdf = new JPDF(jinput);
jpdf->CreatePDF();

JHistos *jhisto = new JHistos();
jhisto->CreateFlowHistos();

//-----------------------------Generating pdfs--------------------------------------
//Making acceptance function fNUE
TString NUEFormula = "[0]*(1-(x > 1.65)*(x < 2.2)*0.5";
if(b2holes)NUEFormula+="-(x > 0.3)*(x < 0.4)*0.7";//background with two gaps in phi
NUEFormula+=")";
// NUEFormula="[0]";
TF1 *fNUE = new TF1("fNUE",NUEFormula,0.0,2.0*TMath::Pi());
fNUE->SetParameter(0,1.0);

//-------Random number needed for sampling
TRandom3 *prng = new TRandom3(random_seed);
gRandom->SetSeed(2*random_seed); //used for GetRandom()
//---------------------------End of generating pdfs---------------------------------

int ieout = Nevt/10;
if (ieout<1) ieout=1;

TStopwatch timer;
timer.Start();

if (bNUE==1) {
	TFile *fRoot = new TFile(sFile, "RECREATE");
	fRoot->cd();
	for (UInt_t iEvent=0; iEvent<Nevt; iEvent++) {
		if(iEvent % ieout == 0) { cout << iEvent << "\t" << int(float(iEvent)/Nevt*100) << "%" << endl ;}
		
		double cent = prng->Uniform(0.0,60.0);
		Int_t ic = GetCentralityBin(cent);
		if(ic < 0)
			continue;
		 
		jhisto->hCentSample->Fill(ic);
		UInt_t Nch=inputNch[ic];
		
		jpdf->GeneratePDF(prng,ic);
		TF1 *fourier = jpdf->GetPDF();

		for (UInt_t t=0; t<Nch; t++){
			double phi = fourier->GetRandom();
			if(prng->Uniform(0,1) > fNUE->Eval(phi,0,0))//For a 1-d function give y=0 and z=0 
				continue;
			jhisto->hSample->Fill(phi);
		}
	}

	jhisto->hSample->Write("hNUA");

	fRoot->Close();
	delete fRoot;
} else if (bNUE==0) {
	TFile *fRoot = new TFile(sFile, "READ");
	TFile *fOutRoot = new TFile(Form("output_%d_evt.root", Nevt), "RECREATE");
	TH1D *hNUA = (TH1D*) fRoot->Get("hNUA");
	TH1D *hCorrected = new TH1D("hPhiCorrected","hPhiCorrected",200, 0.0, 2.0*TMath::Pi());
	for (UInt_t iEvent=0; iEvent<Nevt; iEvent++) {
		if(iEvent % ieout == 0) { cout << iEvent << "\t" << int(float(iEvent)/Nevt*100) << "%" << endl ;}
		
		double cent = prng->Uniform(0.0,60.0);
		Int_t ic = GetCentralityBin(cent);
		if(ic < 0)
			continue;

		UInt_t Nch=inputNch[ic];
		
		jpdf->GeneratePDF(prng,ic);
		TF1 *fourier = jpdf->GetPDF();

		for (UInt_t t=0; t<Nch; t++){
			double phi = fourier->GetRandom();
			if(prng->Uniform(0,1) > fNUE->Eval(phi,0,0))//For a 1-d function give y=0 and z=0 
				continue;
			double weight = hNUA->GetBinContent(hNUA->FindBin(phi));
			hCorrected->Fill(phi,1./weight);
			jhisto->hSample->Fill(phi);
		}
	}

	fOutRoot->cd();
	hCorrected->Write();
	jhisto->hSample->Write();
	fOutRoot->Close();
	delete fOutRoot;
	fRoot->Close();
	delete fRoot;
}


delete prng;
delete fNUE;
delete jhisto;
cout <<"Successfully finished."<<endl;
timer.Print();
return 0;
}

Int_t GetCentralityBin(Double_t cent) {
for(UInt_t i = 0; i < NC; ++i)
	if(cent < centBins[i])
		return i;
return -1;
}