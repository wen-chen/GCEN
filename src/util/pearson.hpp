#ifndef __PEARSON_H_
#define __PEARSON_H_

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

inline double pearson(const std::vector<double>& x,
                      const std::vector<double>& y) {
  constexpr double TINY = std::numeric_limits<double>::min();
  int n = x.size();

  double ax = 0.0, ay = 0;
  for (int i = 0; i < n; i++) {
    ax += x[i];
    ay += y[i];
  }
  ax = ax / n;
  ay = ay / n;

  double xt = 0.0, yt = 0.0, sxx = 0.0, syy = 0.0, sxy = 0.0;
  for (int i = 0; i < n; i++) {
    xt = x[i] - ax;
    yt = y[i] - ay;
    sxx += xt * xt;
    syy += yt * yt;
    sxy += xt * yt;
  }

  double cor = sxy / std::sqrt(sxx * syy + TINY);
  return cor;
}

/*
Regularized Incomplete Beta Function comes from
https://github.com/codeplea/incbeta Lewis Van Winkle
*/
inline double incbeta(const double a, const double b, const double x) {
  const double TINY = 1.0e-30;
  const double STOP = 1.0e-8;

  if (x < 0.0 || x > 1.0) return 1.0 / 0.0;

  /*The continued fraction converges nicely for x < (a+1)/(a+b+2)*/
  if (x > (a + 1.0) / (a + b + 2.0)) {
    return (1.0 -
            incbeta(b, a, 1.0 - x)); /*Use the fact that beta is symmetrical.*/
  }

  /*Find the first part before the continued fraction.*/
  const double lbeta_ab = lgamma(a) + lgamma(b) - lgamma(a + b);
  const double front = exp(log(x) * a + log(1.0 - x) * b - lbeta_ab) / a;

  /*Use Lentz's algorithm to evaluate the continued fraction.*/
  double f = 1.0, c = 1.0, d = 0.0;

  int i, m;
  for (i = 0; i <= 200; ++i) {
    m = i / 2;

    double numerator;
    if (i == 0) {
      numerator = 1.0; /*First numerator is 1.0.*/
    } else if (i % 2 == 0) {
      numerator = (m * (b - m) * x) /
                  ((a + 2.0 * m - 1.0) * (a + 2.0 * m)); /*Even term.*/
    } else {
      numerator = -((a + m) * (a + b + m) * x) /
                  ((a + 2.0 * m) * (a + 2.0 * m + 1)); /*Odd term.*/
    }

    /*Do an iteration of Lentz's algorithm.*/

    d = 1.0 + numerator * d;
    if (fabs(d) < TINY) {
      d = TINY;
    }
    d = 1.0 / d;

    c = 1.0 + numerator / c;
    if (fabs(c) < TINY) {
      c = TINY;
    }

    const double cd = c * d;
    f *= cd;

    /*Check for stop.*/
    if (fabs(1.0 - cd) < STOP) {
      return front * (f - 1.0);
    }
  }

  return 1.0 / 0.0; /*Needed more loops, did not converge.*/
}

inline double get_p_value(const double r, const double n) {
  constexpr double TINY = std::numeric_limits<double>::min();
  double df, t, p;

  df = n - 2;
  t = r * std::sqrt(df / ((1.0 - r + TINY) * (1.0 + r + TINY)));
  p = incbeta(0.5 * df, 0.5, df / (df + t * t));

  return p;
}

#endif
