#define CATCH_CONFIG_MAIN


#include <vector>
#include "third_party/catch.hpp"
#include "../src/util/fdr.hpp"


TEST_CASE("Testing the calculation of FDR", "[default]") {
  std::vector <double> p_values = {0.28757752, 0.40897693, 0.05953272, 0.47189451, 0.44856499, 0.04316666,
                                   0.32242936, 0.10292469, 0.24608773, 0.32792072};
  std::vector <double> fdrs_R = {0.7026873, 0.7078418, 0.4464954, 0.7078418, 0.7078418, 0.4464954, 0.7026873,
                                 0.5146234, 0.7026873, 0.7026873};
  std::vector <double> fdrs(p_values.size(), 0);
  calc_fdr(p_values, fdrs, 15);
  for (unsigned int i = 0; i < fdrs.size(); ++i) {
    REQUIRE(fdrs[i]  == Approx(fdrs_R[i]));
  }
}
