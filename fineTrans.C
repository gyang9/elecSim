#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TTree.h>
#include <TVector3.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TPolyMarker3D.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include "TSystem.h"
#include <TRandom3.h>
#include "TG4Event.h"

#include "kalman-test.cpp"

void fineTransfer( int inputF, int ident){

Double_t E,startX,startY,startZ,stopX,stopY,stopZ,startT,stopT;
Int_t eventN;
std::string det;
Int_t iii=0;

//constants for energy calibration
const double CBIRKS = 0.00208; // mm/MeV
const double EdepToPhotConv_FGD = 70.8; // CLHEP::MeV; // contains both collection in fiber and edep->gamma conversion 
const double DistMPPCscint_FGD = 41; //*CLHEP::mm;
const double LongCompFrac_FGD = 0.816;
const double LongAtt_FGD = 11926.; //*CLHEP::mm;
const double ShortAtt_FGD = 312.; //*CLHEP::mm;
const double DecayLength_FGD = 0.0858; // CLHEP::mm;
const double Lbar_FGD = 1864.3; //* CLHEP::mm;
const double TransTimeInFiber = 1./28. *10.; //mm/ns
// SuperFGD constants
const double MPPCEff_SuperFGD = 0.38;

// Approximate collection factors from PDG2016 (Detectors and accelerators section)
const double CollFactor_SingleClad = 0.06;
//const double CollFactor_DoubleClad = 0.054; // from Licciardi's thesis  
const double CollFactor_DoubleClad = 0.10;

const double Pedestal = 0;//145;  // pedeltal of ADC counts
const double Gain = 10;  // Gain ADC counts of high gain channel
const double LowGain  = 1;  // Gain ADC counts of low gain channel
const double ElecNoise = 1.7;  // sigma of high gain electronics noise
const double LowElecNoise = 1.2;  // sigma of low gain electronics noise
const double PixelGainVari = 0.031;  // gain variation among pixels

double a=0.;        // long attenuation component fraction
double d=0.;        // distance MPPC-scint outside the bar
double LongAtt=0.;  // long attenuation length
double ShortAtt=0.; // short attenuation length
double Ldecay=0.;   // decay length
double Lbar=0.;     // bar length
  
double hitLocation[3000][3]={},hitPE[3000][3]={},hitT[3000][3]={};
double adc_tmp[3000][3]={},loadc_tmp[3000][3]={},Q[3000][3]={},loQ[3000][3]={},adc[3000][3]={};
double loadc[3000][3]={};
Int_t prim[3000],PDG[3000];
double ener[3000];
double trueMom,trueLen;
double true3Mom[3]={};
double true4Mom[100][4]={};
double trueCos[100]={};
double Enu;
int Mode;
int nupdg;
Int_t if3DST[3000]={};
Int_t ifTPC[3000]={};
double vtxPoint[4] = {};
Int_t NHits = 0;
Int_t hitPE_m[3000][3]={};

TFile* outFile = TFile::Open(Form("/dune/app/users/gyang/elecSim/full3DST.neutrino.eleSim.file%d.patch%d.root", inputF, ident),"RECREATE");
TTree* c = new TTree("EDepSimTree","EDepSimTree");
c->Branch("event",&eventN,"event/I");
c->Branch("hitLocation",&hitLocation,"hitLocation[3000][3]/D");
c->Branch("hitPE_mean",&hitPE,"hitPE_mean[3000][3]/D");
c->Branch("hitPE_measure",&hitPE_m,"hitPE_measure[3000][3]/I");
c->Branch("hitT",&hitT,"hitT[3000][3]/D");
c->Branch("hitADC",&adc,"hitADC[3000][3]/D");
c->Branch("hitQ",&Q,"hitQ[3000][3]/D");
c->Branch("hitPrim",&prim,"hitPrim[3000]/I");
c->Branch("hitPDG",&PDG,"hitPDG[3000]/I");
c->Branch("hitE",&ener,"hitE[3000]/D");
c->Branch("trueCos",&trueCos,"trueCos[100]/D");
c->Branch("true4Mom",&true4Mom,"true4Mom[100][4]/D");
c->Branch("if3DST",&if3DST,"if3DST[3000]/I");
c->Branch("ifTPC",&ifTPC,"ifTPC[3000]/I");
c->Branch("vtxPoint",&vtxPoint,"vtxPoint[4]/D");
c->Branch("Mode",&Mode,"Mode/I");
c->Branch("Enu",&Enu,"Enu/D");
c->Branch("nuPDG",& nupdg,"nuPDG/I");
//c->Branch("trueLen",&trueLen,"trueLen/D");
//c->Branch("trueMom",&trueMom,"trueMom/D");
//c->Branch("hitLowQ",&loQ,"hitLowQ[3]/D");
//c->Branch("hitLowADC",&loadc,"hitLowADC[3]/D");

TFile g(Form("/pnfs/dune/persistent/users/gyang/3DST/edep/fullGeo/standardGeo10/PROD106/full3DST.neutrino.%d.edepsim.root",inputF));
TTree* events = (TTree*) g.Get("EDepSimEvents");

TG4Event* event=NULL;
events->SetBranchAddress("Event",&event);
Int_t nevent = 2000; //events->GetEntries();

Int_t StdHepPdgb[10000];
Int_t StdHepStatusb[10000];
Int_t StdHepN;
Double_t vtxb[4]={};
Double_t ivtxb[10000][4]={};
Double_t imomb[10000][4]={};
Bool_t flagb=false;
Bool_t wflagb=false;
Int_t rightSignb=0,wrongSignb=0;
double nuenergy[10000]={};
int interactionMode[10000]={};
int nuPDG[10000]={};
double vtxin[10000][3]={};
TRandom3* random1 = new TRandom3();
TRandom3* random2 = new TRandom3();
random2->SetSeed(6666);

TFile f1b(Form("/pnfs/dune/persistent/users/gyang/3DST/genie/fullGeo/standardGeo10/PROD101/full3DST.neutrino.%d.rootracker.root",inputF));
TTree* h1b = (TTree*) f1b.Get("gRooTracker");

h1b->SetBranchAddress("StdHepPdg", &StdHepPdgb);
h1b->SetBranchAddress("StdHepStatus", &StdHepStatusb);
h1b->SetBranchAddress("EvtVtx", &vtxb);
h1b->SetBranchAddress("StdHepX4", &ivtxb);
h1b->SetBranchAddress("StdHepP4", &imomb);
h1b->SetBranchAddress("StdHepN", &StdHepN);

Int_t nentriesb = (Int_t) h1b->GetEntries();

for(Int_t ib=nevent*ident; ib<nevent*(ident+1) ; ib++) {

  h1b->GetEntry(ib);

  vtxin[ib][0] = vtxb[0];
  vtxin[ib][1] = vtxb[1];
  vtxin[ib][2] = vtxb[2];
  vtxin[ib][3] = vtxb[3];
  
  //std::cout<<"neutrino energy "<<imomb[0][3]<<" vertex "<<vtxin[ib][0]<< std::endl;
  for(Int_t iib=0;iib<3;iib++){
    if(StdHepPdgb[iib]==14){
      nuPDG[ib] = 14;
      nuenergy[ib] = imomb[iib][3]; break;
    }
    if(StdHepPdgb[iib]==-14){
      nuPDG[ib] = -14;
      nuenergy[ib] = imomb[iib][3]; break;
    }
    if(StdHepPdgb[iib]==12){
      nuPDG[ib] = 12;
      nuenergy[ib] = imomb[iib][3]; break;
    }
    if(StdHepPdgb[iib]==-12){
      nuPDG[ib] = -12;
      nuenergy[ib] = imomb[iib][3]; break;
    }
  }
  int cPi = 0; int zPi = 0;
  for(Int_t iib=0; iib<15; iib++){
    if(TMath::Abs(StdHepPdgb[iib])==211){
      cPi ++;
    }
    if(TMath::Abs(StdHepPdgb[iib])==111){
      zPi ++;
    }
  }
  if(cPi==0 || zPi==0) interactionMode[ib] = 1;
  else if(cPi>0 || zPi==0) interactionMode[ib] = 2;
  else if(cPi==0 || zPi>0) interactionMode[ib] = 3;
  else if(cPi>0 || zPi>0) interactionMode[ib] = 4;  
}
/////////////////////////////////////////////////////////////////////////////////////

for(Int_t ii=nevent*ident;ii<nevent*(ident+1);ii++){

  events->GetEntry(ii);
  std::cout<<"event number "<<ii<<"----------------- number of prim. particle "<<event->Primaries[0].Particles.size()<<std::endl;

  double randomNumber = random1->Gaus(0,0.03);
  eventN = ii;
  Enu = nuenergy[ii];
  Mode = interactionMode[ii];
  nupdg = nuPDG[ii];
  vtxPoint[0] = vtxin[ii][0];
  vtxPoint[1] = vtxin[ii][1];
  vtxPoint[2] = vtxin[ii][2];
  vtxPoint[3] = vtxin[ii][3];
  //std::cout<<"re-check neutrino energy "<<Enu<<" and vertex "<<vtxin[ii][0]<<" "<<vtxPoint[0]<<std::endl;

  NHits = 0;

  for (Int_t ccloop1 = 0;ccloop1< 3000 ; ccloop1++){
    prim[ccloop1] = -1;
    PDG[ccloop1] = -1;
    ener[ccloop1] = -1;
    if3DST[ccloop1] = -1;
    ifTPC[ccloop1] = -1;
    if (ccloop1 < 30){
      trueCos[ccloop1]=-1;  
    }
    for(Int_t ccloop2 = 0; ccloop2<3; ccloop2++){
      hitLocation[ccloop1][ccloop2]=-1;
      hitPE[ccloop1][ccloop2]=-1;
      hitPE_m[ccloop1][ccloop2]=-1;
      hitT[ccloop1][ccloop2]=-1;
      adc[ccloop1][ccloop2]=-1;
      Q[ccloop1][ccloop2]=-1;
    }
  }

  for (Int_t ccloop1 = 0;ccloop1< 100 ; ccloop1++){
    for(Int_t ccloop2 = 0; ccloop2<4; ccloop2++){
      true4Mom[ccloop1][ccloop2]=-1;
    }
  }

for(auto sd : event->SegmentDetectors)
{
  for(Int_t i=0;i<sd.second.size();i++){

  /*  T2K upgrade numbers
  --- Cut levels --------
  < highlandCorrections.mom_corrections.sigma_x = 0.800 >   // minimum accum level to save the event
  < highlandCorrections.mom_corrections.B = 0.2 >   // minimum accum level to save the event


  --- Event Weights ------------

  Enable/disable configurations with a single systematic (when EnableSingleWeightSystConfigurations = 1)
  and enable systematics in the "all_syst" configuration (when EnableAllSystematics = 1)

  < highlandCorrections.mom_corrections.x0 = 115.103 > 
  */
    //std::cout<<sd.first<<std::endl;
    TString det = sd.first;
    if( det.Contains("TPC") ){
      /*
      ////////////////////////////////////////////////////////////////////////////
      // this function is used in superFGD	    
      double pt_iv = anaUtils::ComputeInversePT(*cross)*units::MeV;
      double sigma = Sigma(cross,partID);

      Float_t u_t = sqrt(1-pow(cross->DirectionStart[0],2));
      double pt_iv_smeared = gRandom->Gaus(pt_iv,sigma);
      double p_smeared     = (1/pt_iv_smeared)*(1/u_t);
  
      cross->MomentumError=sigma;
      return p_smeared;

      ////////////////////////////////////////////////////////////////////////////////
      // this function is used in superFGD
      double LYZ = cross->DeltaLYZ/units::m;
      double sin_lambda = sqrt(1-pow(cross->DirectionStart[2],2));
      double cos_lambda = fabs(cross->DirectionStart[2]);
  
      int N=0;
      //if (cross->DirectionStart[1] < cross->DirectionStart[2])
      N = abs(cross->DeltaLYZ/(9.7*units::mm)) > 72 ? 72:abs(cross->DeltaLYZ/(9.7*units::mm));
      //else
      //N = abs(cross->DeltaLYZ*cross->DirectionStart[1]/(6.9*units::mm));

      if (LYZ < 1e-6)
        return 0;

      //  if     (pdg==211)  mass = 139.570; // pi+-
      //  else if(pdg==13)   mass = 105.658; // muon
      //  else if(pdg==2212) mass = 938.272; // proton
      //  else if(pdg==11)   mass = 0.511;   // electron
      double mass = anaUtils::GetParticleMass(partID);

      double bg      = cross->Momentum/mass;
      double beta   = bg / sqrt(1. + bg * bg);
      double sigma_x = ND::params().GetParameterD("highlandCorrections.mom_corrections.sigma_x");
      double B       = ND::params().GetParameterD("highlandCorrections.mom_corrections.B");
      double X0      = ND::params().GetParameterD("highlandCorrections.mom_corrections.x0");

      double p_iv = anaUtils::ComputeInversePT(*cross)*units::GeV;

      double sigma_inv_p_point = 1e-3/(300*B*LYZ*LYZ)*sigma_x*sqrt(720/(N+4));
      double sigma_inv_p_MS    = 1e-3*0.016*p_iv/(300*B*LYZ*cos_lambda)*sqrt(LYZ/X0);
      double sigma_inv_p       = sqrt(TMath::Power(sigma_inv_p_point,2)+TMath::Power(sigma_inv_p_MS,2));
      //std::cout << "sigma comp.: " << LYZ << " " << cos_lambda << " " << N << " " << 1/p_iv << " " << sigma_inv_p_point << " " << sigma_inv_p_MS << " " << sigma_inv_p << std::endl;
   
      return sigma_inv_p;
      /////////////////////////////////////////////////////////////////////////////////
      */
      if(sd.second[i].TrackLength>0 && NHits<3000){

        ifTPC[NHits] = 1;

	// point resolution for T2K TPC 0.7 um : https://arxiv.org/pdf/1012.0865.pdf
        double xlocation = random2->Gaus((sd.second[i].Stop.X()+sd.second[i].Start.X())/2.,0.7);
        double ylocation = random2->Gaus((sd.second[i].Stop.Y()+sd.second[i].Start.Y())/2.,0.7);
        double zlocation = random2->Gaus((sd.second[i].Stop.Z()+sd.second[i].Start.Z())/2.-5000. ,0.7);

        prim[NHits] = sd.second[i].PrimaryId;
        PDG[NHits] = event->Primaries[0].Particles[prim[NHits]].PDGCode;

        trueMom = event->Primaries[0].Particles[prim[NHits]].Momentum.Energy();
        trueLen = 0;
        true3Mom[0]=event->Primaries[0].Particles[prim[NHits]].Momentum.Px();
        true3Mom[1]=event->Primaries[0].Particles[prim[NHits]].Momentum.Py();
        true3Mom[2]=event->Primaries[0].Particles[prim[NHits]].Momentum.Pz();

        int primTemp = prim[NHits];
        if(primTemp < 30){
          true4Mom[primTemp][0] = event->Primaries[0].Particles[prim[NHits]].Momentum.Px();
          true4Mom[primTemp][1] = event->Primaries[0].Particles[prim[NHits]].Momentum.Py();
          true4Mom[primTemp][2] = event->Primaries[0].Particles[prim[NHits]].Momentum.Pz();
          true4Mom[primTemp][3] = event->Primaries[0].Particles[prim[NHits]].Momentum.Energy();
          trueCos[primTemp] = event->Primaries[0].Particles[prim[NHits]].Momentum.CosTheta();
        }

        hitLocation[NHits][0]=xlocation;
        hitLocation[NHits][1]=ylocation;
        hitLocation[NHits][2]=zlocation;

        ener[NHits] = sd.second[i].EnergyDeposit;
        NHits++;
      }
    }
    //////////////////////////////////////////////////////////////////////////////////

    else if( det.Contains("Cube") ){
      if(sd.second[i].TrackLength>0 && NHits<3000){

        if3DST[NHits] = 1;

        int aveX= ((sd.second[i].Stop.X()+sd.second[i].Start.X())/2.)/10 ;
        int aveY= ((sd.second[i].Stop.Y()+sd.second[i].Start.Y())/2.)/10 ;
        int aveZ= ((sd.second[i].Stop.Z()+sd.second[i].Start.Z())/2.)/10 ;

        double xlocation = aveX*10. + 5 + 0.;
        double ylocation = aveY*10. + 5 + 0.;
        double zlocation = aveZ*10. + 5 + 500. - 5500;

        prim[NHits] = sd.second[i].PrimaryId;
        PDG[NHits] = event->Primaries[0].Particles[prim[NHits]].PDGCode; 

        trueMom = event->Primaries[0].Particles[prim[NHits]].Momentum.Energy();
        trueLen = 0;
	true3Mom[0]=event->Primaries[0].Particles[prim[NHits]].Momentum.Px();
	true3Mom[1]=event->Primaries[0].Particles[prim[NHits]].Momentum.Py();
	true3Mom[2]=event->Primaries[0].Particles[prim[NHits]].Momentum.Pz();

	int primTemp = prim[NHits];
	if(primTemp < 30){
	  true4Mom[primTemp][0] = event->Primaries[0].Particles[prim[NHits]].Momentum.Px();
	  true4Mom[primTemp][1] = event->Primaries[0].Particles[prim[NHits]].Momentum.Py();
	  true4Mom[primTemp][2] = event->Primaries[0].Particles[prim[NHits]].Momentum.Pz();
	  true4Mom[primTemp][3] = event->Primaries[0].Particles[prim[NHits]].Momentum.Energy();
	  trueCos[primTemp] = event->Primaries[0].Particles[prim[NHits]].Momentum.CosTheta();
	}

	hitLocation[NHits][0]=xlocation;
	hitLocation[NHits][1]=ylocation;
	hitLocation[NHits][2]=zlocation;

	ener[NHits] = sd.second[i].EnergyDeposit;

	Double_t dedx = sd.second[i].EnergyDeposit/sd.second[i].TrackLength;
	Double_t edep= sd.second[i].EnergyDeposit/(1. + CBIRKS*dedx);

        // Account for the 3 fibers in the same scintillator cube
        double collfact = CollFactor_DoubleClad;
        double fact_fib1 = collfact;
        double fact_fib2 = (1-fact_fib1)*collfact;
        double fact_fib3 = (1-fact_fib2)*collfact;
        double CollFactAve = (fact_fib1+fact_fib2+fact_fib3)/3.;
        double NormShadowLight = CollFactAve / collfact; // fraction 
        double Nphot = edep * EdepToPhotConv_FGD * NormShadowLight;

        a = LongCompFrac_FGD;
        d = DistMPPCscint_FGD;
        LongAtt = LongAtt_FGD;
        ShortAtt = ShortAtt_FGD;
        Ldecay= DecayLength_FGD;
        Lbar = Lbar_FGD;

  	double xx = 2400 - xlocation;
  	double yy = 2400 - ylocation;
  	double zz = 2000 - zlocation;  
  	double NphotXY = Nphot * ( a*TMath::Exp((-zz-d)/LongAtt) + (1-a)*TMath::Exp((-zz-d)/ShortAtt) ) * (1/3.);
  	double NphotXZ = Nphot * ( a*TMath::Exp((-yy-d)/LongAtt) + (1-a)*TMath::Exp((-yy-d)/ShortAtt) ) * (1/3.);
  	double NphotYZ = Nphot * ( a*TMath::Exp((-xx-d)/LongAtt) + (1-a)*TMath::Exp((-xx-d)/ShortAtt) ) * (1/3.);

  	double TimeDelayXY =  sd.second[i].Start.T()+TransTimeInFiber * zz;
  	double TimeDelayXZ =  sd.second[i].Start.T()+TransTimeInFiber * yy;
  	double TimeDelayYZ =  sd.second[i].Start.T()+TransTimeInFiber * xx;

  	double peXY = NphotXY * MPPCEff_SuperFGD;
  	double peXZ = NphotXZ * MPPCEff_SuperFGD;
  	double peYZ = NphotYZ * MPPCEff_SuperFGD;

  	hitT[NHits][0]=TimeDelayXY;
  	hitT[NHits][1]=TimeDelayXZ;
  	hitT[NHits][2]=TimeDelayYZ;
  	hitPE[NHits][0]=peXY;
  	hitPE[NHits][1]=peXZ;
  	hitPE[NHits][2]=peYZ;

  	random1->SetSeed(NHits*10);
  	if(peXY>0)
    	  hitPE_m[NHits][0]= (int)random1->Poisson(peXY) ;
  	random1->SetSeed(NHits*11);
  	if(peXZ>0)
    	  hitPE_m[NHits][1]= (int)random1->Poisson(peXZ) ;
  	random1->SetSeed(NHits*12);
  	if(peYZ>0)
    	  hitPE_m[NHits][2]= (int)random1->Poisson(peYZ) ;

  	for(Int_t dim =0;dim<3;dim++){
  	  //PE to ADC
  	  adc_tmp[NHits][dim] = Pedestal + (hitPE[NHits][dim])*Gain;
  	  loadc_tmp[NHits][dim] = Pedestal + (hitPE[NHits][dim])*LowGain*14.29/13.55;

  	  //Electronics noise
  	  adc_tmp[NHits][dim] = random1->Gaus(adc_tmp[NHits][dim],ElecNoise);
  	  loadc_tmp[NHits][dim] = random1->Gaus(loadc_tmp[NHits][dim],LowElecNoise);

  	  //ADC to Charge
  	  //Q=(adc_tmp+53)/217;
  	  //loQ=(loadc_tmp+82)/26;
  	  Q[NHits][dim]=(adc_tmp[NHits][dim])/135.5;
  	  loQ[NHits][dim]=(loadc_tmp[NHits][dim])/14.29;

  	  //Non linearlity of high gain ADC
          if(Q[NHits][dim]<0.65) adc[NHits][dim]=135.5*Q[NHits][dim];
  	  else if(Q[NHits][dim]<3.2)  adc[NHits][dim]=217*Q[NHits][dim]-53;
  	  else if(Q[NHits][dim]<4.2)  adc[NHits][dim]=158.6*Q[NHits][dim]+133.9;
  	  else if(Q[NHits][dim]<14)  adc[NHits][dim]=5.1*Q[NHits][dim]+778.6;
  	  else  adc[NHits][dim]=850;

  	  //Non linearlity of low gain ADC
  	  if(loQ[NHits][dim]<7)  loadc[NHits][dim]=14.29*loQ[NHits][dim];
  	  else if(loQ[NHits][dim]<27)  loadc[NHits][dim]=26*loQ[NHits][dim]-82;
  	  else if(loQ[NHits][dim]<35.5)  loadc[NHits][dim]=21.12*loQ[NHits][dim]+48.24;
  	  else if(loQ[NHits][dim]<178.4)  loadc[NHits][dim]=0.7*loQ[NHits][dim]+775.1;
 	  else  loadc[NHits][dim]=900;
  	}	

  	NHits++;
      }
    }
  }
}
c->Fill();

}
outFile->Write();
outFile->Close();
}


int main(int argc, char* argv[]){

  //for(int ipatch=0;ipatch<1;ipatch++)
  //  fineTransfer(atoi(argv[1]), ipatch);

  std::cout<<"electronics simulation done. "<<std::endl;

  // Do kalman filter
  //
  TH2D* h2[20];
  TH2D* h22[20];
  for(int i=0;i<20;i++){ 
    h2[i] = new TH2D("","",2000,-1000,1000,2400,-1200,1200);
    h22[i] = new TH2D("","",2000,-1000,1000,2400,-1200,1200);
  }

  TH1D* signAll = new TH1D("","",50,0,5000);
  TH1D* signSel = new TH1D("","",50,0,5000);
  std::cout<<"ready to run Kalman filter .."<<std::endl;
  int n = 5; // Number of states
  int m = 3; // Number of measurements

  double dt = 1.0/30; // Time step

  Eigen::MatrixXd A(n, n); // System dynamics matrix
  Eigen::MatrixXd C(m, n); // Output matrix
  Eigen::MatrixXd Q(n, n); // Process noise covariance
  Eigen::MatrixXd R(m, m); // Measurement noise covariance
  Eigen::MatrixXd P(n, n); // Estimate error covariance

  // Discrete LTI projectile motion, measuring position only
  A << 1, 0, 0, 0, 0,
       0, 1, 0, 0, 0,
       0, 0, 1, 0, 0,
       0, 0, 0, 1, 0,
       0, 0, 0, 0, 1;
  C << 1, 0, 0, 0, 0,
       0, 1, 0, 0, 0,
       0, 0, 1, 0, 0;

  // Reasonable covariance matrices
  Q << 1, .05, .05, .0, .0,
       .05, 1, .05, .0, .0,
       .05, .05, 1, .0, .05,
       .0, .0, .0, 1, .05,
       .0, .0, .0, .0, 1,
  R << 1, 0, 0,
       0, 1, 0,
       0, 0, 1;
  P << 100, 1, 1, 1, 1,
       1, 100, 1, 1, 1,
       1, 1, 100, 1, 1,
       1, 1, 1, 100, 1,
       1, 1, 1, 1, 100;

  std::cout << "A: \n" << A << std::endl;
  std::cout << "C: \n" << C << std::endl;
  std::cout << "Q: \n" << Q << std::endl;
  std::cout << "R: \n" << R << std::endl;
  std::cout << "P: \n" << P << std::endl;

  // Construct the filter
  KalmanFilterAlgorithm kf(0, A, C, Q, R, P);

  int usePDG = 13;
  TString useDet = "3DSTTPC";
  // set up measurements
  std::cout<<"reading in the measurement "<<std::endl;
  std::map<double, std::vector<std::vector<double>>> _measurements = kf.readInMeasure(atoi(argv[1]), usePDG, useDet);
  std::vector<std::vector<double>> outR;
  int mCount = 0;
  int iEvt = 0;
  //std::vector<double> measurements;
  std::cout<<"looping through the events "<<std::endl;

  for(std::map<double,std::vector<std::vector<double>>>::iterator Measurements = _measurements.begin(); Measurements != _measurements.end(); ++ Measurements){
  //for(std::vector<std::vector<double>>::iterator iloop = _measurements.begin(); iloop != _measurements.end(); ++ iloop){

    //for(std::vector<double>::iterator jloop = iloop->begin(); jloop != iloop->end(); ++ jloop){
      //measurements.push_back(jloop);
      //std::cout<<jloop<<std::endl;
    //}
    
    double trueE = Measurements->first;
    std::vector<std::vector<double>> measurements = Measurements->second;
    std::cout<<"size of the event hit list: "<<measurements.size()<<std::endl; 
    if(measurements.size()>3){
      // Best guess of initial states
      Eigen::VectorXd beforeUpdate(n);
      Eigen::VectorXd x0(n);
      beforeUpdate<< 0, 0, 1, 0, 0;
      x0 << 0, 0, 1, 0, 0; //-9.81;
      kf.init(0, x0);

      // Feed measurements into filter, output estimated states
      double t              = 0;
      double signVote       = 0;
      double meanFirstOrder = 0.;
      double meanSecondOrder= 0.;
      Eigen::VectorXd y(m);
      std::cout << "t = " << t << ", " << "x_hat[0]: " << kf.state().transpose() << std::endl;
      double latest = 0.;
      for(int i = 0; i < measurements.size(); i++) {
      //for (std::vector<double*>::iterator i = measurements->begin(); i != measurements->end(); i++) {
        t += dt;
	//std::cout<<"doing propagate , beforeUpdate and kf state "<<beforeUpdate.transpose() <<" "<< kf.state().transpose()<<std::endl;
	Eigen::MatrixXd AA = kf.propagate(beforeUpdate, kf.state());
	//std::cout<<"doing projection"<<std::endl;
	Eigen::MatrixXd CC = kf.projection(kf.state(), 0.01);
        y << (measurements.at(i))[0], (measurements.at(i))[1], (measurements.at(i))[2];
	//std::cout<<"measurements "<<(measurements->at(i))[0]<<" "<< (measurements->at(i))[1]<<" "<< (measurements->at(i))[2]<<std::endl;
	//std::cout<<"measurements "<<i[0]<<" "<<i[1]<<" "<<i[2]<<std::endl;
	beforeUpdate = kf.state();
	std::cout<<"updating "<<std::endl;
	kf.update(y, t, AA, CC);
        //std::cout << "t = " << t << ", " << "y[" << i << "] = " << y.transpose()
        //    << ", x_hat[" << i << "] = " << kf.state().transpose() << std::endl;
        std::cout<<"event "<<iEvt<<"; point "<<i<<"; measurements "<<y.transpose()<<"; predictions "<<kf.state().transpose()<<std::endl;
        //std::cout<<"test "<<kf.state().transpose()(0)<<"  |  "<<kf.state().transpose()(1)<<std::endl;
        signVote += kf.state().transpose()(2);
        meanFirstOrder  +=  kf.state().transpose()(1);
        meanSecondOrder +=  kf.state().transpose()(2);
	if (iEvt<20){
	  h2[iEvt] -> Fill((measurements.at(i))[0],(measurements.at(i))[1], 1);
	  Eigen::VectorXd tempVec = CC * kf.state();
	  h22[iEvt]-> Fill(tempVec[0], tempVec[1], 1);
	}	  
        //outR[iEvt].push_back(kf.state().transpose()(2));
	latest = kf.state().transpose()(2);
      }
      iEvt ++;

      meanFirstOrder /= measurements.size();
      meanSecondOrder /= measurements.size();
      int m_muonSign2 = 0;
      signAll -> Fill(trueE);
      if (latest < 0){
	signSel -> Fill(trueE);
        m_muonSign2 = -1;
        mCount++;
      }  
      else{
        m_muonSign2 = 1;
      }
    }
  }
  TFile* outFile = TFile::Open("outfile_0.4.root","RECREATE");
  for(Int_t ii=0; ii< 20; ii++){
    h2[ii]->Write(Form("measure_track_%d",ii));
    h22[ii]->Write(Form("state_track_%d",ii));
  }
  signAll->Write("signAll");
  signSel->Write("signSel");
  outFile->Close();
  std::cout<<"total number of events "<<iEvt<<" ; number of muons with negative sign : "<<mCount<<std::endl;
}

