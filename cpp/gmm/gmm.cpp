#include <iostream>
#include <stdio.h>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cmath>
#include <string>
#include <climits>
#include <algorithm>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <fstream>
#include <algorithm>

using std::min;
using std::max;
using std::numeric_limits;
using std::vector;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::istream;
using std::ofstream;
using std::ostream;
using std::sort;
using std::pair;
using std::make_pair;

struct 
CovarianceMatrix{
  vector<vector<double> > m, m_inv;
  size_t dim;
  double det;
  size_t det_sgn;

  CovarianceMatrix (size_t num_dims){
    dim = num_dims;
    m = vector <vector <double> > (dim, vector <double> (dim, 0.0));
    m_inv = vector <vector <double> > (dim, vector <double> (dim, 0.0));
  }

  // Debug print
  void print(bool inv = false, ostream &os = cout) const {
    for(size_t i = 0; i < dim; ++i){
      for(size_t j = 0; j < dim; ++j)
        if(!inv)
          os << m[i][j]  << "\t";
        else
          os << m_inv[i][j] << "\t";

      os << "\n";
    }
  }

  // Resets covariance to identity
  void reset(){
    for(size_t i = 0; i < dim; ++i)
      for(size_t j = i; j < dim; ++j)
        if(i == j)
          m[i][i] = m_inv[i][i] = 1;
        else
          m[i][j] = m[j][i] = m_inv[i][j] = m_inv[j][i] = 0;
    det_sgn = 1;
    det = 1.0;
  }
  // Finds the inverse of the covariance matrix and its determinant.
  // to be used for density computation and called after a new instance
  // of the covariance matrix has been computed
  bool invert_det(){
    bool invertible = true;

    gsl_matrix *gm = gsl_matrix_alloc(dim,dim);
    for(size_t i = 0; i < dim; ++i)
      for(size_t j = 0; j < dim; ++j)
        gsl_matrix_set(gm, i,j,m[i][j]);

    gsl_permutation *p = gsl_permutation_alloc(dim);
    int s;

    // LU decomposition
    gsl_linalg_LU_decomp(gm, p, &s);
    
    //Log of determinant
    det = gsl_linalg_LU_lndet (gm);
    det_sgn = gsl_linalg_LU_sgndet(gm, s);
       
    // Test if det = 0
    // exp(-310) is roughly underflow for double precision
    if(isnan(det) || det < -300.0){
      reset();
      invertible = false;
    } else {
      // Inverse
      gsl_matrix *inv = gsl_matrix_alloc(dim,dim);
      gsl_linalg_LU_invert(gm, p, inv);

      for(size_t i = 0; i < dim; ++i)
        for(size_t j = i; j < dim; ++j)
          m_inv[j][i] = m_inv[i][j] = gsl_matrix_get(inv,i,j);

      gsl_matrix_free(inv);
    }
    
    gsl_permutation_free(p);
    gsl_matrix_free(gm);
    return invertible;
  }
};

double
dist (const vector <double> &a, const vector <double> &b){
  double ans = 0,z;
  assert(a.size() == b.size());
  for(size_t i = 0; i < a.size(); ++i){
    z = a[i] - b[i];
    ans += z*z;
  }
  return ans;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//// #    #       #     # #######    #    #     #  #####
//// #   #        ##   ## #         # #   ##    # #     #
//// #  #         # # # # #        #   #  # #   # #
//// ###    ##### #  #  # #####   #     # #  #  #  #####
//// #  #         #     # #       ####### #   # #       #
//// #   #        #     # #       #     # #    ## #     #
//// #    #       #     # ####### #     # #     #  #####
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

// Sets the k-th cluster centroid as random points
void
kmeans_init_centroids (vector <vector <double> > &means, 
                       const vector <vector <double> > &in){

  size_t cur_cluster = 0, cur_elems = 0;
  const size_t num_samples = in.size();
  const size_t num_clusters = means.size();
  const size_t num_dims = in[0].size();

  for(size_t k = 0; k < num_clusters; ++k)
    for(size_t j = 0; j < num_dims; ++j)
      means[k][j] = in[k][j];
}

// Assigns elements to a centroid
// returns true if any of them changed
// if false we can stop k-means
bool
kmeans_assign (const vector <vector <double> > &in,
               const vector <vector <double>> &means,
               vector <size_t> &which_cluster){
  bool changed = false;
  for(size_t i = 0; i < in.size(); i++){
    double min_dist = dist(in[i], means[0]),
           cur_dist;
    size_t which = 0;
    for(size_t k = 1; k < means.size(); k++){
      cur_dist = dist(in[i], means[k]);
      if(cur_dist< min_dist){
        min_dist = cur_dist;
        which = k;
      }
    }

    if(which != which_cluster[i])
      changed = true;
    which_cluster[i] = which;
  }

  return changed;
}

// Updates centroids based on assignment
void
kmeans_update(const vector <vector <double> > &in,
              vector <vector <double> > &means,
              const vector <size_t> &which_cluster){

  const size_t num_dims = in[0].size(),
               num_clusters = means.size(),
               num_samples = in.size();
            
  vector <size_t> num_elements = vector<size_t>(num_clusters,0);

  // Set the means back to zero
  for(size_t k = 0; k < num_clusters;++k)
    means[k] = vector<double>(num_dims, 0.0);

  for(size_t i = 0; i < num_samples; ++i){
    for(size_t j = 0; j < num_dims; ++j){
      means[which_cluster[i]][j] += in[i][j];
    }
    num_elements[which_cluster[i]]++;
  }

  for(size_t k = 0; k < num_clusters; ++k)
    for(size_t j = 0; j < num_dims; ++j)
      means[k][j] = means[k][j]/double(num_elements[k]);
}

////////////////////////////
////////////////////////////
///  #####  #     # #     # 
/// #     # ##   ## ##   ## 
/// #       # # # # # # # # 
/// #  #### #  #  # #  #  # 
/// #     # #     # #     # 
/// #     # #     # #     # 
///  #####  #     # #     # 
////////////////////////////
////////////////////////////

// MVN density (log)
double 
dnorm (const vector <double> &x,
                const vector <double> &mean,
                const CovarianceMatrix &cov){

  size_t num_dims = mean.size();
  double exp_val = 0.0;

  // L = (x - mu)^t * Sigma^{-1} 
  vector <double> left_mult = vector <double> (num_dims, 0.0);
  for(size_t i = 0; i < num_dims; ++i)
    for(size_t j = 0; j < num_dims; ++j)
      left_mult[i] += (x[j] - mean[j])*cov.m_inv[j][i];

  // R = L * (x - mu)
  for(size_t j = 0; j < num_dims; ++j)
    exp_val += left_mult[j]*(x[j] - mean[j]);
  
  // If matrix is singular or not psd, probability = 0
  if(exp_val < 0.0)
    return numeric_limits<double>::min();

  return (-0.5 * (cov.det + double(num_dims)*log(2*M_PI) + exp_val));
}
                
// Calculates the log sum of exponentials by making sure
// all numbers exponentiated are negative to avoid overflow
// More exactly, calculates log(\sum_i \pi_k N_k (x_i))
double
log_sum_of_exps(const vector <double> &pi,
                const vector <double> &log_density){
  double ans = 0.0;
  size_t which_max = 0;
  for(size_t k = 1; k < log_density.size(); ++k)
    if(log_density[k] > log_density[which_max])
      which_max = k;

  for(size_t k = 0; k < log_density.size(); ++k)
    ans += pi[k]*exp(log_density[k] - log_density[which_max]);

  return (log_density[which_max] + log(ans));
}

// Computes the responsibilities but with some numerical tricks
// to avoid underflow
void
e_step(vector <vector <double> > &gamma,
       const vector <vector <double> > &in,
       const vector <double> &pi,
       const vector <vector <double> > &means,
       const vector <CovarianceMatrix> &cov ){

  const size_t num_clusters = pi.size();
  const size_t num_samples = in.size();

  // Summation of small exps
  for(size_t i = 0; i < num_samples; ++i){
    vector <double> log_density = vector<double>(num_clusters,0.0);
    for(size_t k = 0; k < num_clusters; ++k)
      log_density[k] = dnorm(in[i], means[k], cov[k]);

    double sum = log_sum_of_exps(pi, log_density);
    for(size_t k = 0; k < num_clusters; ++k){
      if(pi[k] == 0.0) {
        gamma[i][k] = 0.0;
      } else {
        gamma[i][k]  = log(pi[k]) + log_density[k];
        gamma[i][k] -= sum;
        gamma[i][k]  = exp(gamma[i][k]);
      }
      
      assert(!isnan(gamma[i][k]));
    }
  }
}

// Maximum likelihood params
void
m_step(const vector< vector <double> > &gamma,
       const vector< vector <double> > &in,
       vector <double> &pi,
       vector <vector <double> > &means,
       vector <CovarianceMatrix> &cov){
  const size_t num_clusters = pi.size();
  const size_t num_samples = in.size();
  const size_t num_dims = in[0].size();

  for(size_t k = 0; k < num_clusters; ++k){

    // Pi
    pi[k] = 0.0;
    for(size_t i = 0; i < num_samples; ++i)
      pi[k] += gamma[i][k];

    // Unused cluster, send it back to the center
    if(pi[k] <= 0.0){
      for(size_t j = 0; j < num_dims; ++j)
        means[k][j] = 0.0;
      cov[k].reset();

    // Some elements have been assigned, estiamte mean and cov
    } else {

      // Means
      for(size_t j = 0; j < num_dims; ++j){
        means[k][j] = 0.0;
        for(size_t i = 0; i < num_samples; ++i)
          means[k][j] += gamma[i][k]*in[i][j];       
        means[k][j] = means[k][j] / pi[k];
      }

      // Sigma
      for(size_t j = 0; j < num_dims; ++j){
        for(size_t jp = j; jp < num_dims; ++jp){
          cov[k].m[j][jp] = cov[k].m[jp][j] = 0.0;

          for(size_t i = 0; i < num_samples; ++i){
            double z = gamma[i][k]*(in[i][j] - means[k][j])*(in[i][jp] - means[k][jp]);
            cov[k].m[j][jp] += z;
            if(jp != j)
              cov[k].m[jp][j] += z;
          }
          cov[k].m[j][jp] = cov[k].m[jp][j] = cov[k].m[j][jp]/(pi[k]);
        }
      }  
      //Dont forget this!!!
      cov[k].invert_det();

      // pi goes to [0,1]
      pi[k] = pi[k]/double(num_samples);
    }
  }
}

// Total log likelihood of data given EM parameters
double
log_likelihood(const vector <vector <double> > &in,
               const vector <vector <double> > &means,
               const vector <CovarianceMatrix> &cov,
               const vector <double> &pi){
  double loglik = 0.0;
  const size_t num_clusters = pi.size();
  const size_t num_samples = in.size();

 for(size_t i = 0; i < num_samples; ++i){
   double sum = 0;
   vector<double> log_density = vector<double>(num_clusters, 0.0);
   for(size_t k = 0; k < num_clusters; ++k)
     log_density[k] = dnorm(in[i], means[k], cov[k]);  

   loglik += log_sum_of_exps(pi, log_density);
 }

 return loglik;
}

// Bayesian Information Criterion
double
information_criterion(const vector <vector <double> > &in,
    const vector <vector <double> > &means,
    const vector <CovarianceMatrix> &cov,
    const vector <double> &pi,
    bool bayesian = true){

  size_t num_params = 0;
  size_t num_samples = in.size(),
         num_clusters = pi.size(),
         num_dims = in[0].size();

  num_params += num_clusters * num_dims; //means
  num_params += num_clusters - 1; //pi. Minus one cause they sum to one
  num_params += num_clusters*num_dims*(num_dims + 1) >> 1; //cov

  double ans = double(num_params);
  if(bayesian)
    ans *= log(double(num_samples));
  else 
    ans *= 2;

  return ans - 2*log_likelihood(in,means,cov,pi);
}

void 
print_debug(const vector <vector <double> > &gamma, 
            const vector <double> &pi,
            const vector <vector <double> > &means,
            const vector <CovarianceMatrix> &cov){
  const size_t num_clusters = pi.size();
  const size_t num_dims = means[0].size();

  for(size_t k = 0; k < num_clusters; ++k){
    cout << "=========CLUSTER " << k << "========";
    cout << "\n---pi:\t" << pi[k];
    cout << "\n---mean:\t";
    for(size_t j = 0; j < num_dims; ++j)
      cout << means[k][j] << "\t";

    cout << "--cov det:\t" << cov[k].det << "\n";
    //cout << "\n---cov\n";
    //cov[k].print();
    //cout << "\n---cov-inv\n";
    //cov[k].print(true);
  }

  cout << "========GAMMA==========\n";
  for(size_t i = 0; i < gamma.size(); ++i){
    for(size_t k = 0; k < num_clusters; ++k){
      if(gamma[i][k] > 1e-5)
        cout << gamma[i][k] << "\t";
      else
        cout << "0\t";
    }
    cout << "\n";
  }
}

////////////////////////////////
////////////////////////////////
/// #     #    #    ### #     #
/// ##   ##   # #    #  ##    #
/// # # # #  #   #   #  # #   #
/// #  #  # #     #  #  #  #  #
/// #     # #######  #  #   # #
/// #     # #     #  #  #    ##
/// #     # #     # ### #     #
////////////////////////////////
////////////////////////////////

// Returns the clustering BIC
double
cluster (const vector <vector <double> > &in, 
         const size_t num_clusters,
         const vector <size_t> &kmeans_init){
  const size_t num_samples = in.size(),
               num_dims = in[0].size();

  vector <vector <double> > means = vector <vector <double> > (num_clusters, vector <double> (num_dims, 0.0));

  /////////////////////////////////////////////////
  // GMM algorithm to soften the k-means assignment
  /////////////////////////////////////////////////
  
  // mixing probabilities
  vector <double> pi = vector <double> (num_clusters, 0.0);

  // Responsibilites
  vector <vector <double> > gamma = vector <vector <double> > (num_samples, vector <double> (num_clusters, 0.0));

  // Covariance matrices
  vector <CovarianceMatrix> cov = vector <CovarianceMatrix> (num_clusters, CovarianceMatrix(num_dims));

  // Read the cluster assignment from k means
  for(size_t i = 0; i < num_samples; ++i)
    gamma[i][kmeans_init[i]] = 1.0;

  // Initial guess of parameters given k means
  m_step(gamma,in,pi,means,cov);

  double prev = log_likelihood(in,means,cov,pi),cur;
  // Some steps to test if loglik is increasing
  
  double eps = 1e-2;
  bool stop = false;
  for(size_t nsteps = 0; !stop; ++nsteps){
    e_step(gamma,in,pi,means,cov);
    m_step(gamma,in,pi,means,cov);
    cur = log_likelihood(in,means,cov,pi);
    if(fabs(prev - cur) < eps)
      stop = true;
    prev = cur;
  }
  
  //print_debug(gamma,pi,means,cov);
  double ans = information_criterion(in,means,cov,pi, false);
  return ans;
}

int 
main(){
  const size_t max_cl = 50;
  size_t num_samples, num_dims;
  cin >> num_samples >> num_dims;

  vector <vector <double> > in  = vector < vector <double> > (num_samples, vector <double> (num_dims, 0.0));
  for(size_t i = 0; i < num_samples; i++)
    for(size_t j = 0; j < num_dims; j++)
      cin >> in[i][j];

  vector <vector <size_t> > kmeans_init  = vector <vector <size_t> > (max_cl, vector <size_t> (num_samples, 0));
  for(size_t k = 0; k < max_cl; ++k)
    for(size_t i = 0; i < num_samples; ++i)
      cin >> kmeans_init[k][i];

  vector<double> ics = vector<double> (max_cl, 0.0);

#ifdef _OPENMP
  #pragma omp parallel for 
#endif
 for(size_t cl = 1; cl <= max_cl; ++cl){
    cerr << "Calculating k = " << cl << "...\n";
    ics[cl - 1] = cluster(in,cl, kmeans_init[cl - 1]);
  }

  size_t which_min = 0;
  for(size_t i = 0; i < max_cl; ++i)
    cout << i+1 << "\t" << ics[i] << "\n";
  
  for(size_t i = 1; i < max_cl; ++i)
    if(ics[i] < ics[which_min])
      which_min = i;

  cout << "Number of clusters: " << which_min + 1 << endl;
}

