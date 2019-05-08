#include <vector>
#include <Eigen/Dense>

class KalmanFilterAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    KalmanFilterAlgorithm();
  /**
  * Create a Kalman filter with the specified matrices.
  *   A - System dynamics matrix
  *   C - Output matrix
  *   Q - Process noise covariance
  *   R - Measurement noise covariance
  *   P - Estimate error covariance
  */
   
    KalmanFilterAlgorithm(
        double kdt,
        const Eigen::MatrixXd& kA,
        const Eigen::MatrixXd& kC,
        const Eigen::MatrixXd& kQ,
        const Eigen::MatrixXd& kR,
        const Eigen::MatrixXd& kP
    );

  // read in the electronics simulated files
  std::map<double,std::vector<std::vector<double>>> readInMeasure(int t, int pdg, TString det);

  // Matrices for computation
  Eigen::MatrixXd kA, kC, kQ, kR, kP, kK, kP0;

  // System dimensions
  int m, n;

  // Initial and current time
  double kt0, kt;

  // Discrete time step
  double kdt;

  // Is the filter initialized?
  bool kinitialized;

  // n-size identity
  Eigen::MatrixXd kI;

  // Estimated states
  Eigen::VectorXd kx_hat, kx_hat_new;


  /**
  * Create a blank estimator.
  */
  //KalmanFilterAlgorithm();

  /**
  * Initialize the filter with initial states as zero.
  */
  //void init();

  /**
  * Initialize the filter with a guess for initial states.
  */
  void init(double t0, const Eigen::VectorXd& x0);

  /**
  * Update the estimated state based on measured values. The
  * time step is assumed to remain constant.
  */
  void update(const Eigen::VectorXd& y);

  /**
  * Update the estimated state based on measured values,
  * using the given time step and dynamics matrix.
  */
  void update(const Eigen::VectorXd& y, double dt, const Eigen::MatrixXd A, const Eigen::MatrixXd C);

  /**
  * Return the current state and time.
  */
  Eigen::VectorXd state() { return kx_hat; };
  double time() { return kt; };

  Eigen::MatrixXd currentPropagator(double drho, double theta0, double kappa, double dz, double tanlambda, Eigen::VectorXd& x0, double theta);

  Double_t x0_prime(double x0, double drho, double theta0, double alpha, double kappa, double theta);

  Double_t y0_prime(double x0, double drho, double theta0, double alpha, double kappa, double theta);

  Double_t z0_prime(double x0, double dz, double alpha, double kappa, double tanlambda, double theta);

  Double_t drho_prime(double* xx, double* par);

  Double_t theta0_prime(double* xx, double* par);

  Double_t kappa_prime(double* xx, double* par);

  Double_t dz_prime (double* xx, double* par);

  Double_t tanlambda_prime (double* xx, double* par);

  Eigen::MatrixXd projection(double drho, double theta0, double kappa, double dz, double tanlambda, double theta);

  Eigen::MatrixXd propagate(double drho, double theta0, double kappa, double dz, double tanlambda,
                double drho_prime, double theta0_prime, double kappa_prime, double dz_prime, double tanlambda_prime);

  Eigen::MatrixXd projection(Eigen::VectorXd v, double theta);

  Eigen::MatrixXd propagate(Eigen::VectorXd v1, Eigen::VectorXd v2);

  double depositEnergy ;
  double trueEnergy;
  //std::vector<double*> measurements;

  std::string         m_inputCaloHitListName;             ///< The input calo hit list name
  std::string         m_outputCaloHitListName;           ///< The output calo hit list name for TPC_VIEW_U hits
  std::string         m_currentCaloHitListReplacement;    ///< The name of the calo hit list to replace the current list (optional)
  std::string         m_outputPfoListName;

};



