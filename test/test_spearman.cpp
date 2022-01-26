#define CATCH_CONFIG_MAIN

#include <vector>
#include "../src/util/base.hpp"
#include "../src/util/pearson.hpp"
#include "third_party/catch.hpp"

TEST_CASE("Testing the calculation of Spearman's correlation coefficient",
          "[default]") {
  std::vector<double> vec1 = {1, 2, 3, 4, 5};
  std::vector<double> vec2 = {5, 6, 7, 8, 7};
  std::vector<double> ranked_vec1 = get_rank(vec1);
  std::vector<double> ranked_vec2 = get_rank(vec2);
  double corr = pearson(ranked_vec1, ranked_vec2);
  double p = get_p_value(corr, vec1.size());

  REQUIRE(corr == Approx(0.8207826816681233));
  REQUIRE(p == Approx(0.08858700531354381));
}
