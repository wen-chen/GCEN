#define CATCH_CONFIG_MAIN


#include <vector>
#include "../src/util/pearson.hpp"
#include "third_party/catch.hpp"

TEST_CASE("Testing the calculation of Pearson's correlation coefficient", "[default]") {
  std::vector <double> vec1 = {1.2, 3.4, 5.6, 7.8, 6.9};
  std::vector <double> vec2 = {4.6, 5.8, 7.9, 10.4, 11.8};
  double corr = pearson(vec1, vec2);
  double p = get_p_value(corr, vec1.size());

  REQUIRE(corr == Approx(0.9345719941493214));
  REQUIRE(p == Approx(0.019891640636871065));
}
