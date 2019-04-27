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


  std::vector<std::vector<double>> readInMeasure(TString t, int pdg, TString det);

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
  void update(const Eigen::VectorXd& y, double dt, const Eigen::MatrixXd A);

  /**
  * Return the current state and time.
  */
  Eigen::VectorXd state() { return kx_hat; };
  double time() { return kt; };

  double depositEnergy ;
  double trueEnergy;
  std::vector<double> measurements;

  std::string         m_inputCaloHitListName;             ///< The input calo hit list name
  std::string         m_outputCaloHitListName;           ///< The output calo hit list name for TPC_VIEW_U hits
  std::string         m_currentCaloHitListReplacement;    ///< The name of the calo hit list to replace the current list (optional)
  std::string         m_outputPfoListName;

};



