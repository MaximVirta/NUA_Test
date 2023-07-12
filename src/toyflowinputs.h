#ifndef JTreeGenerator_H
#define JTreeGenerator_H

const int NC = 7; // # of centralities
const double centBins[NC] = {5.0,10.0,20.0,30.0,40.0,50.0,60.0};
const TString strCentrality[NC] = {"0-5\%","5-10\%","10-20\%","20-30\%","30-40\%","40-50\%","50-60\%"};
const float inputNch[NC] ={1943, 1603, 1196, 802, 521, 322, 185};
const int NH = 3; // # of harmonic orders
//rows: input vn from n=2 to n=7
//cols: vn for multiplicities from central to peripheral
// mean value of vn, needs to be combined with the fluctuations
const double inputVn[NH][NC] = {
	{0.,0.,0.,0.,0.,0.,0.},
	{0.0279896,0.0449772,0.0645474,0.0839821,0.0959965,0.10125,0.0988201},
	{0.0206594,0.0240024,0.0269388,0.0298739,0.0316964,0.0318984,0.0303094}
    /*{0.044325430932513066,0.13409520177544662,0.11868377868652139},
    {0.020319043993158488,0.06147004658105161,0.05440535759432017},
    {0.009314371910439144,0.028178238867966162,0.02493974296844666},
    {0.0042697640752778244,0.01291707408507019,0.011432528097143131},
    {0.001957285518962443,0.005921264409071094,0.0052407396041462785},
    {0.0008972314477330478,0.0027143431996458607,0.0024023865382260113},*/
};

#endif