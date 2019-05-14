/**
 *  @file   ExampleAlgorithms/KalmanFilterAlgorithm.cc
 *
 *  @brief  Implementation of the list preparation algorithm class.
 *
 *  $Log: $
 */

#include "kalman-test.h"
#include <TF1.h>
#include <fstream>

//------------------------------------------------------------------------------------------------------------------------------------------

std::map<double, std::vector<std::vector<double>>> KalmanFilterAlgorithm::readInMeasure(int inFile, int pdg, TString det){

  TFile infile(Form("/home/guang/work/elecSim/full3DST.neutrino.eleSim.file%d.patch0.root",inFile));
  TTree* t = (TTree*)infile.Get("EDepSimTree");

  std::cout<<"inside the readInMeasure() "<<std::endl;
  double hitLocation[3000][3]={};
  double hitPE_mean[3000][3]={};
  double hitPE_measure[3000][3]={};
  double hitT[3000][3]={};
  double hitPrim[3000]={};
  int hitPDG[3000]={};
  double hitE[3000]={};
  double true4Mom[100][4]={};
  int if3DST[3000]={};
  int ifTPC[3000]={};
  double vtxPoint[4]={};
  int use3DST=0;
  int useTPC=0;
  double trueE = 0;
  std::map<double, std::vector<std::vector<double>>> output;

  t->SetBranchAddress("hitLocation",&hitLocation);
  t->SetBranchAddress("hitE",&hitE);
  t->SetBranchAddress("hitPDG",&hitPDG);
  t->SetBranchAddress("if3DST",&if3DST);
  t->SetBranchAddress("ifTPC",&ifTPC);
  t->SetBranchAddress("true4Mom",&true4Mom);
  t->SetBranchAddress("hitPrim",&hitPrim);
  t->SetBranchAddress("vtxPoint",&vtxPoint);

  if(det.Contains("TPC")){
    useTPC = 1;
  }
  if(det.Contains("3DST")){
    use3DST = 1;
  }

  std::cout<<"require pdg "<<pdg<<std::endl;
  //std::vector<double*> temp;
  //double tem[3];
  for(Int_t i=0;i<t->GetEntries();i++){

    t->GetEntry(i);

    for(Int_t j=0;j<3000;j++){
      if(hitPDG[j] == pdg && hitPrim[j] > 0 ){
        int temploca = hitPrim[j];
	if(true4Mom[temploca][3]){
	  trueE = true4Mom[temploca][3];
	  break;
	}	
      }
    }
    std::cout<<"true energy is "<<trueE<<" and vtx "<<vtxPoint[0]<<" "<<vtxPoint[1]<<" "<<vtxPoint[2]<<std::endl;

    std::vector<std::vector<double>> temp;
    std::vector<double> tem;
    if (trueE>300 && TMath::Abs(vtxPoint[0])<0.8 && TMath::Abs(vtxPoint[1])<0.8 && vtxPoint[2]>4.2 && vtxPoint[2]<5.8){
      for(Int_t j=0;j<3000;j++){
        if(if3DST[j] == use3DST && use3DST == 1){
          if(hitPDG[j] == pdg && j%3 == 1){
	    // the coordinate in KF is different, so use this:	
	    tem.push_back( hitLocation[j][1]);
	    tem.push_back( hitLocation[j][2]);
	    tem.push_back( hitLocation[j][0]);
	    //std::cout<< tem[0]<<std::endl;
            temp.push_back(tem);    
	  }
        }
        else if(ifTPC[j] == useTPC && useTPC == 1){
          if(hitPDG[j] == pdg && j%3 == 1){
            // the coordinate in KF is different, so use this:    
            tem.push_back( hitLocation[j][1]);
            tem.push_back( hitLocation[j][2]);
            tem.push_back( hitLocation[j][0]);
            temp.push_back(tem);
	  }
        }
        tem.erase (tem.begin(), tem.end());  
      }
      //for(Int_t itest = 0; itest< temp.size(); itest++)
      //  std::cout<<temp.size()<<" "<<(temp.at(itest))[0]<<" "<<std::endl; //<<(temp.at(itest))[1]<<" "<<(temp.at(itest))[2]<<std::endl;
      std::cout<<"------------------------------------"<<std::endl;
      //output.emplace(trueE, temp);
      output.insert(make_pair(trueE, temp));
      //output.push_back(temp);
    }
    trueE = -1;
    temp.erase (temp.begin(),temp.end());
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


//http://www-jlc.kek.jp/subg/offl/kaltest/doc/ReferenceManual.pdf
Eigen::MatrixXd KalmanFilterAlgorithm::projection(double drho, double theta0, double kappa, double dz, double tanlambda, double theta){
  Eigen::MatrixXd AA(5,3);
  double alpha = 1./(speedOfLight * Bfield);
  AA<< TMath::Cos(theta0), TMath::Sin(theta0), 0,
       -(drho + alpha/kappa)*TMath::Sin(theta0) + (alpha/kappa)*TMath::Sin(theta0+theta),
       (drho + alpha/kappa)*TMath::Cos(theta0) - (alpha/kappa)*TMath::Cos(theta0+theta), 0,
       -(alpha/(kappa*kappa))* (TMath::Cos(theta0) - TMath::Cos(theta0-theta)),
       -(alpha/(kappa*kappa))*(TMath::Sin(theta0) - TMath::Sin(theta0-theta)), (alpha/(kappa*kappa))*theta*tanlambda,
       0, 0, 1,
       0, 0, -(alpha/kappa)*theta;
  return AA.transpose();
}

Eigen::MatrixXd KalmanFilterAlgorithm::projection(Eigen::VectorXd v, double theta){
  return projection(v(0), v(1), v(2), v(3), v(4), theta);
}


Eigen::MatrixXd KalmanFilterAlgorithm::propagate(double drho, double theta0, double kappa, double dz, double tanlambda, 
		double drho_prime, double theta0_prime, double kappa_prime, double dz_prime, double tanlambda_prime){
  Eigen::MatrixXd AA(5,5);
  double alpha = 1./(speedOfLight * Bfield);
  double aa = -TMath::Power((drho_prime + alpha/kappa),-1) * TMath::Sin(theta0_prime - theta0);
  double ab = (drho + alpha/kappa)*TMath::Power((drho_prime + alpha/kappa),-1) * TMath::Cos(theta0_prime - theta0);
  double ac = alpha/(kappa*kappa) * TMath::Power((drho_prime + alpha/kappa),-1) * TMath::Sin(theta0_prime - theta0);
  double ad = 0;
  double ae = 0;

  //std::cout<<"finished one propagate "<<std::endl;
  double ba = TMath::Cos(theta0_prime - theta0);
  double bb = (drho_prime + alpha/kappa) * TMath::Sin(theta0_prime - theta0);
  double bc = alpha/(kappa*kappa) * (1- TMath::Cos(theta0_prime - theta0));
  double bd = 0;
  double be = 0;

  //std::cout<<"finished two propagate "<<std::endl;
  double ca = 0;
  double cb = 0;
  double cc = 1;
  double cd = 0;
  double ce = 0;

  double ea = 0;
  double eb = 0;
  double ec = 0;
  double ed = 0;
  double ee = 1;
 
  double da = (alpha/kappa) * TMath::Power(drho_prime + alpha/kappa,-1) * tanlambda * TMath::Sin(theta0_prime - theta0);
  double db = (alpha/kappa) * tanlambda * (1 - (drho+alpha/kappa) * TMath::Power(drho_prime + alpha/kappa,-1) * TMath::Cos(theta0_prime - theta0));
  double dc = alpha/(kappa*kappa) * tanlambda * (theta0_prime - theta0 - (alpha/kappa) * TMath::Power(drho_prime + alpha/kappa,-1) * TMath::Sin(theta0_prime - theta0));
  double dd = 1; 
  double de = -(alpha/kappa) * (theta0_prime - theta0);

  //std::cout<<"finished all propagate "<<std::endl;
  AA << ba, bb, bc, bd, be,
        aa, ab, ac, ad, ae,
	ca, cb, cc, cd, ce,
        da, db, dc, dd, de,
	ea, eb, ec, ed, ee;
  //std::cout<<AA<<std::endl;
  return AA;
}

Eigen::MatrixXd KalmanFilterAlgorithm::propagate(Eigen::VectorXd v, Eigen::VectorXd v_prime){
  return propagate (v(0), v(1), v(2), v(3), v(4), v_prime(0), v_prime(1), v_prime(2), v_prime(3), v_prime(4) );
}	


// https://ww2.odu.edu/~skuhn/NucPhys/KalmanFilter.pdf
// http://www-jlc.kek.jp/subg/offl/lib/docs/helix_manip/node3.html
// kappa = Q/Pt ; rho = alpha/kappa; alpha = 1/cB
// with Q being the charge, $\rho$ being the signed radius of the helix, and $\alpha \equiv 1/c B$ being a magnetic-field-dependent constant, dz is the distrance of the helix from the pivotal point in the z direction, and $\tan\lambda$ is the dip angle. 
// $d_\rho$ is the distance of the helix from the pivotal point in the xy plane, $\phi_0$ is the azimuthal angle to specifies the pivotal point with respect to the helix center, $\kappa$ is the signed reciprocal transverse momentum

Eigen::MatrixXd KalmanFilterAlgorithm::currentPropagator(double drho, double theta0, double kappa, double dz, double tanlambda, Eigen::VectorXd& x0, double theta){
  
  double alpha = 1./ (speedOfLight * Bfield);

  double Xc = x0(0) + (drho + alpha/kappa)* TMath::Cos(theta0);
  double Yc = x0(1) + (drho + alpha/kappa)* TMath::Sin(theta0);	  
  
  double x_prime = x0_prime(x0(0), drho, theta0, alpha, kappa, theta);
  double y_prime = y0_prime(x0(1), drho, theta0, alpha, kappa, theta);
  double z_prime = z0_prime(x0(2), dz, alpha, kappa, tanlambda, theta);

  double par[5] = {};
  double xx[5] = {};

  par[0] = drho;
  par[1] = theta0;
  par[2] = kappa;
  par[3] = dz;
  par[4] = tanlambda;

  xx[0] = x0(0);
  xx[1] = x0(1);
  xx[2] = alpha;
  xx[3] = theta;
  xx[4] = theta0_prime(xx, par);
  double tt = xx[4];

  TF1* f1 = new TF1("f1", theta0_prime(xx, par), -100000, 100000, 5);
  Double_t derivative1[5]={};
  f1->GradientPar(xx, derivative1);

  TF1* f2 = new TF1("f2", drho_prime(xx, par), -100000, 100000, 5);
  Double_t derivative2[5]={};
  f2->GradientPar(xx, derivative2);

  TF1* f3 = new TF1("f3", kappa_prime(xx, par), -100000, 100000, 5);
  Double_t derivative3[5]={};
  f3->GradientPar(xx, derivative3);

  xx[0] = x0(2);
  xx[1] = alpha;
  xx[2] = tt;
  xx[3] = theta;

  TF1* f4 = new TF1("f4", dz_prime(xx, par), -100000, 100000, 5);
  Double_t derivative4[5]={};
  f4->GradientPar(xx, derivative4);

  TF1* f5 = new TF1("f5", tanlambda_prime(xx, par), -100000, 100000, 5);
  Double_t derivative5[5]={};
  f5->GradientPar(xx, derivative5);  

  Eigen::MatrixXd AA(5, 5); 
  AA<< derivative1[0], derivative1[1], derivative1[2], derivative1[3], derivative1[4],
       derivative2[0], derivative2[1], derivative2[2], derivative2[3], derivative2[4],
       derivative3[0], derivative3[1], derivative3[2], derivative3[3], derivative3[4],
       derivative4[0], derivative4[1], derivative4[2], derivative4[3], derivative4[4],
       derivative5[0], derivative5[1], derivative5[2], derivative5[3], derivative5[4];
  return AA;
}

Double_t KalmanFilterAlgorithm::x0_prime(double x0, double drho, double theta0, double alpha, double kappa, double theta){
  return x0 + drho*TMath::Cos(theta0) + (alpha/kappa)*(TMath::Cos(theta0)- TMath::Cos(theta0+theta));
}

Double_t KalmanFilterAlgorithm::y0_prime(double x0, double drho, double theta0, double alpha, double kappa, double theta){
  return x0 + drho*TMath::Sin(theta0) + (alpha/kappa)*(TMath::Sin(theta0)- TMath::Sin(theta0+theta));
}	

Double_t KalmanFilterAlgorithm::z0_prime(double x0, double dz, double alpha, double kappa, double tanlambda, double theta){
  return x0 + dz - (alpha/kappa) * tanlambda * theta;
}

Double_t KalmanFilterAlgorithm::drho_prime(double* xx, double* par){
  double drho = par[0];
  double theta0 = par[1];
  double kappa = par[2];
  double dz = par[3];
  double tanlambda = par[4];

  double x0 = xx[0];
  double y0 = xx[1];
  double alpha = xx[2];
  double theta = xx[3];
  double theta0_prime = xx[4];

  double Xc = x0 + (drho + alpha/kappa)* TMath::Cos(theta0);
  double Yc = y0 + (drho + alpha/kappa)* TMath::Sin(theta0);
  double x0_prime = x0 + drho*TMath::Cos(theta0) + (alpha/kappa)*(TMath::Cos(theta0)- TMath::Cos(theta0+theta));
  double y0_prime = y0 + drho*TMath::Sin(theta0) + (alpha/kappa)*(TMath::Sin(theta0)- TMath::Sin(theta0+theta));  
  return (Xc - x0_prime)* TMath::Cos(theta0_prime) + (Yc - y0_prime) * TMath::Sin(theta0_prime) - kappa/alpha;
}

Double_t KalmanFilterAlgorithm::theta0_prime(double* xx, double* par){

  double drho = par[0];
  double theta0 = par[1];
  double kappa = par[2];
  double dz = par[3];
  double tanlambda = par[4];

  double x0 = xx[0];
  double y0 = xx[1];
  double alpha = xx[2];
  double theta = xx[3];

  double Xc = x0 + (drho + alpha/kappa)* TMath::Cos(theta0);
  double Yc = y0 + (drho + alpha/kappa)* TMath::Sin(theta0);
  double x0_prime = x0 + drho*TMath::Cos(theta0) + (alpha/kappa)*(TMath::Cos(theta0)- TMath::Cos(theta0+theta));
  double y0_prime = y0 + drho*TMath::Sin(theta0) + (alpha/kappa)*(TMath::Sin(theta0)- TMath::Sin(theta0+theta));

  return 1./ (TMath::Tan((y0_prime - Yc)/(x0_prime - Xc)));
}

Double_t KalmanFilterAlgorithm::kappa_prime(double* xx, double* par){
  double kappa = par[2];
  return kappa;
}

Double_t KalmanFilterAlgorithm::dz_prime (double* xx, double* par){
  double drho = par[0];
  double theta0 = par[1];
  double kappa = par[2];
  double dz = par[3];
  double tanlambda = par[4];

  double z0 = xx[0];
  double alpha = xx[1];
  double theta0_prime = xx[2];
  double theta = xx[3];    

  double z0_prime = z0 + dz - (alpha/kappa) * tanlambda * theta;

  return z0 - z0_prime + dz - (alpha/kappa) * (theta0_prime - theta0) * tanlambda;
}

Double_t KalmanFilterAlgorithm::tanlambda_prime (double* xx, double* par){
  double tanlambda = par[4];
  return tanlambda;
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
  
  // predict
  kx_hat_new = kA * kx_hat;
  kP = kA*kP*kA.transpose() + kQ;

  // update
  kK = kP*kC.transpose()*(kC*kP*kC.transpose() + kR).inverse();
  kx_hat_new += kK * (y - kC*kx_hat_new);
  kP = (kI - kK*kC)*kP;
  kx_hat = kx_hat_new;

  kt += kdt;
}

void KalmanFilterAlgorithm::update(const Eigen::VectorXd& y, double dt, const Eigen::MatrixXd A, const Eigen::MatrixXd C) {

  this->kA = A;
  this->kdt = dt;
  this->kC = C;
  update(y);
}

