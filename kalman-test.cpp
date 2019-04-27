/**
 *  @file   ExampleAlgorithms/KalmanFilterAlgorithm.cc
 *
 *  @brief  Implementation of the list preparation algorithm class.
 *
 *  $Log: $
 */

#include "kalman-test.h"
#include <fstream>

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<std::vector<double>> KalmanFilterAlgorithm::readInMeasure(TString inFile, int pdg, TString det){

  TFile infile(inFile);
  TTree* t = (TTree*)infile.Get("EDepSimTree");

  double hitLocation[2000][3]={};
  double hitPE_mean[2000][3]={};
  double hitPE_measure[2000][3]={};
  double hitT[2000][3]={};
  double hitPrim[2000]={};
  int hitPDG[2000]={};
  double hitE[2000]={};
  double true4Mom[30][4]={};
  int if3DST[2000]={};
  int ifTPC[2000]={};
  int use3DST=0;
  int useTPC=0;
  std::vector<std::vector<double>> output;

  t->SetBranchAddress("hitLocation",&hitLocation);
  t->SetBranchAddress("hitE",&hitE);
  t->SetBranchAddress("hitPDG",&hitPDG);
  t->SetBranchAddress("if3DST",&if3DST);
  t->SetBranchAddress("ifTPC",&ifTPC);

  if(det.Contains("TPC")){
    useTPC = 1;
  }
  if(det.Contains("3DST")){
    use3DST = 1;
  }

  for(Int_t i=0;i<t->GetEntries();i++){

    t->GetEntry(i);
    for(Int_t j=0;j<2000;j++){
      if(if3DST[j] == use3DST && if3DST[j] == 1){
        if(hitPDG[j] == pdg)
          output[i].push_back(hitLocation[j][1]);    
      }
      else if(ifTPC[j] == useTPC && ifTPC[j] == 1){
        if(hitPDG[j] == pdg)
          output[i].push_back(hitLocation[j][1]);
      }  
    }
  }
  return output;
}


KalmanFilterAlgorithm::KalmanFilterAlgorithm(
    double kdt,
    const Eigen::MatrixXd& kA,
    const Eigen::MatrixXd& kC,
    const Eigen::MatrixXd& kQ,
    const Eigen::MatrixXd& kR,
    const Eigen::MatrixXd& kP)
  : kA(kA), kC(kC), kQ(kQ), kR(kR), kP0(kP),
    m(kC.rows()), n(kA.rows()), kdt(kdt), kinitialized(false),
    kI(n, n), kx_hat(n), kx_hat_new(n)
{
  kI.setIdentity();
}


void KalmanFilterAlgorithm::init(double t0, const Eigen::VectorXd& x0) {
  kx_hat = x0;
  kP = kP0;
  this->kt0 = t0;
  kt = t0;
  kinitialized = true;
}


void KalmanFilterAlgorithm::update(const Eigen::VectorXd& y) {

  if(!kinitialized)
    throw std::runtime_error("Filter is not initialized!");

  kx_hat_new = kA * kx_hat;
  kP = kA*kP*kA.transpose() + kQ;
  kK = kP*kC.transpose()*(kC*kP*kC.transpose() + kR).inverse();
  kx_hat_new += kK * (y - kC*kx_hat_new);
  kP = (kI - kK*kC)*kP;
  kx_hat = kx_hat_new;

  kt += kdt;
}

void KalmanFilterAlgorithm::update(const Eigen::VectorXd& y, double dt, const Eigen::MatrixXd A) {

  this->kA = A;
  this->kdt = dt;
  update(y);
}

