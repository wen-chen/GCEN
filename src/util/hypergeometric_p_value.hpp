#ifndef __HYPERGEOMETRIC_H_
#define __HYPERGEOMETRIC_H_


#include <iostream>
#include <cmath>
#include <cstdlib>


inline void initLogFacs(double* logFacs, const int n) {
  logFacs[0] = 0;
  for(int i = 1; i < n + 1; ++i) {
    logFacs[i] = logFacs[i - 1] + log((double)i); // only n times of log() calls
  }
}


inline double logHypergeometricProb(double* logFacs, const int a, const int b, const int c, const int d) {
  return logFacs[a + b] + logFacs[c + d] + logFacs[a + c] + logFacs[b + d] - logFacs[a] - logFacs[b] - logFacs[c] - logFacs[d] - logFacs[a + b + c + d];
}


double fisher_exact(int a, int b, int c, int d) {
  int n = a + b + c + d;
  double* logFacs = new double[n + 1]; // *** dynamically allocate memory logFacs[0..n] ***
  initLogFacs(logFacs, n); // *** initialize logFacs array ***
  double logpCutoff = logHypergeometricProb(logFacs, a, b, c, d); // *** logFacs added
  double pFraction = 0;
  for(int x = 0; x <= n; ++x) {
    if ( a + b - x >= 0 && a + c - x >= 0 && d - a + x >=0 ) {
      double l = logHypergeometricProb(logFacs, x, a + b - x, a + c - x, d - a + x);
      if ( l <= logpCutoff ) {
        pFraction += exp(l - logpCutoff);
      }
    }
  }
  double logpValue = logpCutoff + log(pFraction);
  delete [] logFacs;
  return exp(logpValue);
}


double calc_p_hypergeometric(int pop_n, int pop_count, int study_n, int study_count) {
  int a = study_count;
  int b = study_n - study_count;
  int c = pop_count - study_count;
  int d = pop_n - pop_count - b;
  return fisher_exact(a, b, c, d);
}

#endif