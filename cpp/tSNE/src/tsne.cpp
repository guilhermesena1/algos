/* tsne: what this program does
 *
 * Copyright (C) 2019 Guilherme De Sena Brandine
 *                    Andrew D. Smith
 *
 * Authors: Guilherme De Sena Brandine
 *          Andrew D. Smith
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include <iostream>
#include <numeric>
#include <cmath>
#include <random>
#include <vector>
#include <cassert>
#include <fstream>
#include <cstdlib>

using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;

// ============================================================
// ================== PRINT AND DEBUG =========================
// ============================================================

// Read matrix from file
// Input file has two numbers r and c in the first row (rows and columns)
// the next r rows have c numbers
static void
read_mtx(const string &filename, vector<vector<double> > &v) {
  std::ifstream is(filename);
  size_t n_cells, n_dims;
  is >> n_cells >> n_dims;
  v = vector<vector<double> >(n_cells, vector<double>(n_dims, 0.0));
  for (size_t i = 0; i < n_cells; ++i)
    for (size_t j = 0; j < n_dims; ++j)
      is >> v[i][j];
}

// Preallocates all matrices
static void
allocate(const size_t n_cells, const size_t low_dim,
         const double eta_init,
         vector<vector<double> > &p,
         vector<vector<double> > &q,
         vector<double> &eta,
         vector<double> &deltabar,
         vector<double> &grad,
         vector<double> &denom) {

  // Vector
  eta = vector<double>(low_dim, eta_init);
  deltabar = vector<double>(low_dim,  0.0);
  grad = vector<double> (low_dim,0.0);
  denom = vector<double>(n_cells,0.0);

  // Matrices
  p = vector<vector<double> >(n_cells, vector<double> (n_cells,0.0));
  q = vector<vector<double> >(n_cells, vector<double> (n_cells,0.0));
}

// ============================================================
// ======== AUXILIARY DISTANCE/PROB FUNCTIONS =================
// ============================================================

// Computes the square Euclidean distance between two vectors
static double
sq_dist(const vector<double> &a, const vector<double> &b) {
  size_t sa = a.size(), sb = b.size();
  assert(sa == sb);
  double ans = 0.0, _tmp;
  for(size_t i = 0; i < sa; ++i){
    _tmp = a[i] - b[i];
    ans += (_tmp)*(_tmp);
  }
  return ans;
}

// ============================================================
// ======================= ENTROPY AND KL =====================
// ============================================================

static double
get_average_nn(const vector<vector<double> > &p){
  double sum_nn = 0.0, tmp;
  for(size_t i = 0; i < p.size(); ++i){
    tmp = 0.0;
    for(size_t j = 0; j < p[0].size(); ++j){
      if(i != j){
        if(p[i][j] > 0) //0log(0) = 0
          tmp -= p[i][j] * log2(p[i][j]);
      }
    }
    sum_nn += pow(2, tmp);
  }

  return sum_nn/p.size();
}

// KL divergence between two matrices
static double
kl(const vector<vector<double> > &p, const vector<vector<double> > &q) {
  double ans = 0;
  size_t sz = p.size();
  for (size_t i = 0 ; i < sz; ++i) {
    for (size_t j = 0; j < sz; ++j) {
      if (i != j) {
        if (p[i][j] > 0)
          ans += p[i][j] * (log(p[i][j]) - log(q[i][j]));
      }
    }
  }
  return ans;
}


// Computes the conditional probability that a high-dimension point was generated
// from the distribution of one of the other points with uniform prior
static void
cond_prob (vector<vector<double> > &p,
           const vector<vector<double> > &cells,
           const double sigma,
           const bool symmetrize = true){

  size_t n_cells = cells.size();
  const double denom = -2*sigma*sigma;


  // Calculates row sums and exp euclidean distance
  vector<double> rowsums = vector<double>(n_cells, 0.0);
  for(size_t i = 0; i < n_cells - 1; ++i){
    for(size_t j = i+1; j < n_cells; ++j){
      double d = exp(sq_dist(cells[i], cells[j]) / denom);
      rowsums[i] += d;
      rowsums[j] += d;
      p[i][j] = p[j][i] = d;
    }
  }

  // Now divide rows by row sums
  for(size_t i = 0; i < n_cells; ++i){
    for(size_t j = 0; j < n_cells; ++j){
      if(rowsums[i]>0)
        p[i][j] = exp(log(p[i][j]) - log(rowsums[i]));
    }
  }

  // Only symmetrize when not searching for sigma
  if(symmetrize){
    for(size_t i = 0; i < n_cells - 1; ++i){
      for(size_t j = i+1; j < n_cells; ++j){
        p[i][j] = p[j][i] = (p[i][j] + p[j][i])/(2*n_cells);
      }
    }
  }
}

// Finds sigma given perplexity
static double
find_sigma(vector<vector<double> > &p,
           const double perplexity,
           const vector<vector<double> > &v){
  double sigma_min = 0.00001;
  double sigma_max = 100;
  double mid;
  // double ent;
  double average_nn;

  while(sigma_max - sigma_min > 0.0000001){
    mid = (sigma_min + sigma_max)/2.0;
    cond_prob(p, v, mid, false);

    average_nn = get_average_nn(p);

    // Entropy too high -> decrease sigma
    if(average_nn > perplexity){
      sigma_max = mid;
    } else {
      sigma_min = mid;
    }
  }

  return mid;
}
//
// ============================================================
// ============================ Y =============================
// ============================================================

// Random initial guess for Y
static void
y_init(const size_t n_cells, const size_t dim, const double sd,
       vector<vector<double> > &y) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> d(0, sd);

  y = vector<vector<double> >(n_cells, vector<double>(dim));
  for (size_t i = 0; i < n_cells; ++i)
    for (size_t j = 0; j < dim; ++j)
      y[i][j] = d(gen);
}

// Conditional distribution of a set of t-distributed points
static void  y_q (vector<vector<double> > &q,
                  const vector<vector<double> > &y){

  size_t n_cells = y.size();
  double sum = 0.0;

  for(size_t i = 0; i < n_cells - 1; ++i){
    for(size_t j = i+1; j < n_cells; ++j){
      q[i][j] = q[j][i] = 1/(1 + sq_dist(y[i],y[j]));
      sum += 2*q[i][j];
    }
  }

  for(size_t i = 0; i < n_cells - 1; ++i){
    for(size_t j = i+1; j < n_cells; ++j){
      if(i != j){
        q[i][j] = q[j][i] = q[i][j]/sum;
       }
    }
  }
}

// ============================================================
// =================== GRADIENT DESCENT =======================
// ============================================================

// Equation 5 from the tSNE paper
static void
gradient (vector<double> &grad,
          vector<double> &denom,
          const vector<vector<double> > &p,
          const vector<vector<double> > &q,
          const vector<vector<double> > &y,
          size_t i){
  size_t low_dim = y[0].size(),
    n_cells = p.size();

  for(size_t j = 0; j < n_cells; j++)
    denom[j] = (1 + sq_dist(y[i],y[j]));

  for(size_t k = 0; k < low_dim; ++k){
    grad[k] = 0.0;
    for(size_t j = 0; j < n_cells; ++j){
      grad[k] += (p[i][j] - q[i][j])*(y[i][k] - y[j][k])/denom[j];
    }
  }
}

//Delta-bar-delta method described by Jacobs, 1988
static void
update_eta (double &eta, double &deltabar, const double &grad){
  const double kappa = 100.0, phi = 0.5, theta = 0.5;
  double increment = 0.0;

  if(grad * deltabar > 0.0) {
    increment = kappa;
  } else if (grad * deltabar < 0.0){
    increment = -phi*eta;
  }

  eta += increment;
  deltabar = (1 - theta)*grad + theta*deltabar;
}

// Computes next y based on gradient descent
static void 
next_y(vector<vector<double> > &y,
       const vector<vector<double> > &p,
       const vector<vector<double> > &q,
       vector<double> &eta,
       vector<double> &deltabar,
       vector<double> &grad,
       vector<double> &denom){

  size_t n_cells = y.size(), n_dims = y[0].size();
  for(size_t i = 0; i < n_cells; ++i){
    gradient(grad,denom,p,q,y,i);
    for(size_t j = 0; j < n_dims; ++j){
      update_eta(eta[j], deltabar[j], grad[j]);
      y[i][j] -= eta[j]*grad[j];
    }
  }
}

// ===========================================================
// ================= MAIN FOR TESTING ========================
// ===========================================================

int
main(int argc, const char **argv) {
  try {

    double perplexity = 10.0;
    double sigma = 0.0;
    double low_sigma = 1.0;
    double eta_init = 100.0;

    size_t low_dim = 2;

    /* FILES */
    string outfile;
    bool VERBOSE = false;

    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse(strip_path(argv[0]), "", "<matrix>");
    opt_parse.add_opt("output", 'o', "output file (default: stdout)",
                      false , outfile);
    opt_parse.add_opt("perp", 'p', "perplexity", false , perplexity);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false , VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string matrix_file(leftover_args.front());
    /**********************************************************************/

    // Matrices
    vector<vector<double> > v;
    // Fills the v vector with input matrix
    read_mtx (matrix_file, v);
    const size_t n_cells = v.size();
    if (VERBOSE)
      cerr << "n_cells=" << n_cells << endl;

    // Vectors
    vector<double> eta, deltabar, grad, denom;
    vector<vector<double> > p,q,y;
    // Preallocates all matrices
    allocate(n_cells, low_dim, eta_init,
             p, q,
             // y,
             eta, deltabar, grad, denom);

    // Binary search to convert perplexity to sigma
    if (VERBOSE)
      cerr << "Finding sigma for perplexity = " << perplexity << "...\n";
    sigma = find_sigma(p,perplexity,v);
    if (VERBOSE)
      cerr << "Perplexity = " << perplexity << " -> Sigma = " << sigma << "\n";

    // Calculates conditional probabilities from high dimension data
    cond_prob(p, v, sigma);

    // Initial random solution for y
    // TODO: implement an educated guess
    y_init(n_cells, low_dim, low_sigma, y);

    size_t n_iter = 0;
    double prev_kl = 1e9; // MAGIC!
    double cur_kl = 1e5; // MAGIC!
    for (; cur_kl < prev_kl; n_iter++) {
      // Calculates q matrix
      y_q(q,y);

      // Calculates kl divergence between p and q
      const double tmp_kl = kl(p,q);
      if (n_iter % 100 == 0) {
        prev_kl = cur_kl;
        cur_kl = tmp_kl;
        if (VERBOSE)
          cerr << "iteration=" << n_iter << endl
               << "kl_divergence= " << cur_kl << endl;
      }
      next_y(y, p, q, eta, deltabar, grad, denom);
    }

    if (VERBOSE)
      cerr << "iteration=" << n_iter << endl
           << "kl_divergence= " << cur_kl << endl;

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    for (size_t i = 0; i < y.size(); ++i) {
      // out << row_names[i];
      copy(begin(y[i]), end(y[i]),
           std::ostream_iterator<double>(out, "\t"));
      out << endl;
    }
  }
  catch (std::exception &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
