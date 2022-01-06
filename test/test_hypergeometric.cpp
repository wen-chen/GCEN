#define CATCH_CONFIG_MAIN


#include "third_party/catch.hpp"
#include "../src/util/hypergeometric_p_value.hpp"


TEST_CASE("Testing the calculation of p-value of hypergeometric", "[default]") {
  REQUIRE(calc_p_hypergeometric(100, 10, 20, 4)  == Approx(0.10957185095927705));
  REQUIRE(calc_p_hypergeometric(100, 10, 20, 6)  == Approx(0.003933076466791367));
  REQUIRE(calc_p_hypergeometric(100, 10, 20, 8)  == Approx(2.3782749640379144e-05));
  REQUIRE(calc_p_hypergeometric(100, 10, 20, 10)  == Approx(1.0673177187555421e-08));
}
