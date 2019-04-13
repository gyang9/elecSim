{
//#include "TG4Event.h"

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
  
double hitLocation[1000][3]={},hitPE[1000][3]={},hitT[1000][3]={};
double adc_tmp[1000][3]={},loadc_tmp[1000][3]={},Q[1000][3]={},loQ[1000][3]={},adc[1000][3]={};
double loadc[1000][3]={};
Int_t prim[1000],PDG[1000];
double ener[1000];
double trueMom,trueLen;
double true3Mom[3]={};
double true4Mom[30][4]={};
double trueCos[30]={};
double Enu;
double Mode;
Int_t if3DST[1000]={};
Int_t intPoint =0;
Int_t NHits = 0;
Int_t hitPE_m[2000][3]={};

TFile* outFile = TFile::Open("/dune/data/users/gyang/elecSim/outputTest.root","RECREATE");
TTree* c = new TTree("EDepSimTree","EDepSimTree");
c->Branch("event",&eventN,"event/I");
c->Branch("hitLocation",&hitLocation,"hitLocation[1000][3]/D");
c->Branch("hitPE_mean",&hitPE,"hitPE_mean[1000][3]/D");
c->Branch("hitPE_measure",&hitPE_m,"hitPE_measure[1000][3]/I");
c->Branch("hitT",&hitT,"hitT[1000][3]/D");
c->Branch("hitADC",&adc,"hitADC[1000][3]/D");
c->Branch("hitQ",&Q,"hitQ[1000][3]/D");
c->Branch("hitPrim",&prim,"hitPrim[1000]/I");
c->Branch("hitPDG",&PDG,"hitPDG[1000]/I");
c->Branch("hitE",&ener,"hitE[1000]/D");
c->Branch("trueCos",&trueCos,"trueCos[30]/D");
c->Branch("true4Mom",&true4Mom,"true4Mom[30][4]/D");
c->Branch("if3DST",&if3DST,"if3DST[1000]/I");
//c->Branch("intPoint",&intPoint,"intPoint/I");
//c->Branch("Mode",&Mode,"Mode/D");
//c->Branch("Enu",&Enu,"Enu/D");
//c->Branch("trueLen",&trueLen,"trueLen/D");
//c->Branch("trueMom",&trueMom,"trueMom/D");
//c->Branch("hitLowQ",&loQ,"hitLowQ[3]/D");
//c->Branch("hitLowADC",&loadc,"hitLowADC[3]/D");

//TFile g("../3DST_nu+e_222.root");
//TFile g("../3DST_222_particleGun1GeVMuon.root");
TFile g(Form("/home/guang/work/elecSim/PROD101/full3DST.neutrino.899.edepsim.root"));
TTree* events = (TTree*) g.Get("EDepSimEvents");

//TGeoManager* geo = new TGeoManager;
//geo->Import("../Juan.gdml");

TG4Event* event=NULL;
events->SetBranchAddress("Event",&event);

Int_t nevent = 200;//events->GetEntries();


Int_t StdHepPdgb[1000];
Int_t StdHepStatusb[1000];
Double_t vtxb[4]={};
Double_t ivtxb[1000][4]={};
Double_t imomb[1000][4]={};
Bool_t flagb=false;
Bool_t wflagb=false;
Int_t rightSignb=0,wrongSignb=0;
double nuenergy[10000]={};
int interactionMode[10000]={};
TRandom3* random1 = new TRandom3();


for(Int_t ii=nevent*atoi(gApplication->Argv(7));ii<nevent*(atoi(gApplication->Argv(7))+1);ii++){

events->GetEntry(ii);
std::cout<<"event number "<<ii<<"----------------- number of prim. particle "<<event->Primaries[0].Particles.size()<<std::endl;

double randomNumber = random1->Gaus(0,0.03);
eventN = ii;
Enu = nuenergy[ii];
Mode = interactionMode[ii];

NHits = 0;

for (Int_t ccloop1 = 0;ccloop1< 1000 ; ccloop1++){
  prim[ccloop1] = 0;
  PDG[ccloop1] = 0;
  ener[ccloop1] = 0;
  if3DST[ccloop1] = 0;
  if (ccloop1 < 30){
    trueCos[ccloop1]=0;  
  }
  for(Int_t ccloop2 = 0; ccloop2<3; ccloop2++){
    hitLocation[ccloop1][ccloop2]=0;
    hitPE[ccloop1][ccloop2]=0;
    hitPE_m[ccloop1][ccloop2]=0;
    hitT[ccloop1][ccloop2]=0;
    adc[ccloop1][ccloop2]=0;
    Q[ccloop1][ccloop2]=0;
  }
}

for (Int_t ccloop1 = 0;ccloop1< 30 ; ccloop1++){
  for(Int_t ccloop2 = 0; ccloop2<4; ccloop2++){
    true4Mom[ccloop1][ccloop2]=0;
  }
}

for(auto sd : event->SegmentDetectors)
{
for(Int_t i=0;i<sd.second.size();i++){
//for (std::vector<TG4HitSegment>::iterator h = sd.second.begin(); h != sd.second.end();++h){

if(sd.first != "3DScint"){

if(i==0){
intPoint = 1;
}
else {intPoint = 0;}

// 3% resolution for ECAL 
ener[NHits] = ((std::vector<TG4HitSegment>)(sd.second))[i].EnergyDeposit * (1 + randomNumber);
//ener[NHits] = h.EnergyDeposit * (1 + randomNumber);

prim[NHits] = ((TG4HitSegment)(sd.second[i])).PrimaryId;
//prim[NHits] = h.PrimaryId;
PDG[NHits] = event->Primaries[0].Particles[prim[NHits]].PDGCode;

//int aveX= (1200 + 1199)/10 ;
//int aveY= (1200 + 1199)/10 ;
//int aveZ= (1000 + 999)/10 ;

// with 3DST+ECAL.gdml
double xlocation = -1000; //aveX*10. + 5 + 0.;
double ylocation = -1000; //aveY*10. + 5 + 55000.;
double zlocation = -1000; //aveZ*10. + 5 - 594000.;

hitLocation[NHits][0]=xlocation;
hitLocation[NHits][1]=ylocation;
hitLocation[NHits][2]=zlocation;

//std::cout<<"in the detector other than 3DScint.. PDG: "<<PDG<<std::endl;
//c->Fill();
if3DST[NHits] = 0;

NHits++;

}

if(sd.first == "3DScint"){
if(sd.second[i].TrackLength>0){
//if(h.TrackLength>0)

if(i==0){
intPoint = 1;
}
else {intPoint = 0;}

if3DST[NHits] = 1;

int aveX= ((sd.second[i].Stop.X()+sd.second[i].Start.X())/2.)/10 ;
int aveY= ((sd.second[i].Stop.Y()+sd.second[i].Start.Y())/2.)/10 ;
int aveZ= ((sd.second[i].Stop.Z()+sd.second[i].Start.Z())/2.)/10 ;

// with 3DST_beam_222.gdml
//double xlocation = aveX*10. + 5 + 0.;
//double ylocation = aveY*10. + 5 + 55000.;
//double zlocation = aveZ*10. + 5 - 591450.;

// with 3DST+ECAL.gdml
double xlocation = aveX*10. + 5 + 0.;
double ylocation = aveY*10. + 5 + 0.;
double zlocation = aveZ*10. + 5 + 500.;

//std::cout<<sd.second[i].Start.X()<<" "<<sd.second[i].Start.Y()<<" "<<sd.second[i].Start.Z()<<std::endl;
//std::cout<<detNumber<<" "<< (sd.second[i].Stop.Z()+sd.second[i].Start.Z())/2.<<" "<<xlocation<<" "<<ylocation<<" "<<zlocation<<std::endl;

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
    //cout << "NormShadowLight = " << NormShadowLight << endl;   
    double Nphot = edep * EdepToPhotConv_FGD * NormShadowLight;
    //std::cout<<"before dividing Nphot: "<<Nphot<<std::endl;

  a = LongCompFrac_FGD;
  d = DistMPPCscint_FGD;
  LongAtt = LongAtt_FGD;
  ShortAtt = ShortAtt_FGD;
  Ldecay= DecayLength_FGD;
  Lbar = Lbar_FGD;

  double xx = 2400 - xlocation;
  double yy = 2400 - ylocation;
  double zz = 2000 - zlocation;  
  //std::cout<<"xx "<<xx<<std::endl;
  double NphotXY = Nphot * ( a*TMath::Exp((-zz-d)/LongAtt) + (1-a)*TMath::Exp((-zz-d)/ShortAtt) ) * (1/3.);
  double NphotXZ = Nphot * ( a*TMath::Exp((-yy-d)/LongAtt) + (1-a)*TMath::Exp((-yy-d)/ShortAtt) ) * (1/3.);
  double NphotYZ = Nphot * ( a*TMath::Exp((-xx-d)/LongAtt) + (1-a)*TMath::Exp((-xx-d)/ShortAtt) ) * (1/3.);

  //std::cout<<"a xx d longAtt ShortAtt xNphot "<<a<<" "<<xx<<" "<<d<<" "<<LongAtt<<" "<<ShortAtt<<" "<<xNphot<<" "<<TMath::Exp((-xx-d)/LongAtt)<<std::endl;
  //std::cout<<"short att "<<TMath::Exp((-xx-d)/ShortAtt)<<" "<<(1-a)*TMath::Exp((-xx-d)/ShortAtt)<<" "<<a*exp((-zz-d)/LongAtt)<<std::endl;
  double TimeDelayXY =  sd.second[i].Start.T()+TransTimeInFiber * zz;
  double TimeDelayXZ =  sd.second[i].Start.T()+TransTimeInFiber * yy;
  double TimeDelayYZ =  sd.second[i].Start.T()+TransTimeInFiber * xx;

  double peXY = NphotXY * MPPCEff_SuperFGD;
  double peXZ = NphotXZ * MPPCEff_SuperFGD;
  double peYZ = NphotYZ * MPPCEff_SuperFGD;

  //std::cout<<xNphot<<" "<<xpe<<std::endl;

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

  //std::cout<<"time AND hitPE: "<<hitT[NHits][0]<<" "<<hitT[NHits][1]<<" "<<hitT[NHits][2]<<" "<<hitPE_m[NHits][0]<<" "<<hitPE_m[NHits][1]<<" "<<hitPE_m[NHits][2]<<std::endl;

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

//  std::cout<<"location, PE, Q, ADC along x: "<<hitLocation[0]<<" "<<hitPE[0]<<" "<<Q[0]<<" "<<adc[0]<<" "<<std::endl; 
//  std::cout<<"location, PE, Q, ADC along y: "<<hitLocation[1]<<" "<<hitPE[1]<<" "<<Q[1]<<" "<<adc[1]<<" "<<std::endl;
//  std::cout<<"location, PE, Q, ADC along z: "<<hitLocation[2]<<" "<<hitPE[2]<<" "<<Q[2]<<" "<<adc[2]<<" "<<std::endl;
/*
det=sd.first;
startX=sd.second[i].Start.X();
startY=sd.second[i].Start.Y();
startZ=sd.second[i].Start.Z();
startT=sd.second[i].Start.T();
stopX=sd.second[i].Stop.X();
stopY=sd.second[i].Stop.Y();
stopZ=sd.second[i].Stop.Z();
stopT=sd.second[i].Stop.T();
*/
//c->Fill();
//events->Fill();
}
}
}
}
c->Fill();

}
//c->AddFriend(events);
outFile->Write();
outFile->Close();
}
