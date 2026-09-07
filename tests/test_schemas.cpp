#include "schema_checks.hpp"

#include "filaments/filaments_schema.hpp"
#include "interpol/interpol_schema.hpp"
#include "partrac/partrac_schema.hpp"
#include "spatial/static_space_stepper_schema.hpp"

#include "filaments_felbmRK4/filaments_felbmRK4_schema.hpp"
#include "filaments_triangleRK4/filaments_triangleRK4_schema.hpp"
#include "omp_test/omp_test_schema.hpp"
#include "tracers_triangleRK4/tracers_triangleRK4_schema.hpp"
#include "tracertensors_triangleRK4/tracertensors_triangleRK4_schema.hpp"
#include "tracervectors_analyticRK4/tracervectors_analyticRK4_schema.hpp"
#include "tracervectors_triangleRK4/tracervectors_triangleRK4_schema.hpp"
#include "tracervectors_triangle_spatial/tracervectors_triangle_spatial_schema.hpp"
#include "tracervectors_trianglefreqRK4/tracervectors_trianglefreqRK4_schema.hpp"
#include "weighted_walkers/weighted_walkers_schema.hpp"

// One translation unit for every app; the namespace keeps the two families apart
TEST_CASE("every app declares a usable schema", "[schema]") {
  for (auto s : {filaments_schema(), interpol_schema(),
                 partrac_schema(), spatial_schema(),
                 filaments_felbmRK4_schema(), filaments_triangleRK4_schema(),
                 omp_test_schema(), tracers_triangleRK4_schema(),
                 tracertensors_triangleRK4_schema(),
                 tracervectors_analyticRK4_schema(),
                 tracervectors_triangleRK4_schema(),
                 tracervectors_triangle_spatial_schema(),
                 tracervectors_trianglefreqRK4_schema(),
                 weighted_walkers_schema()})
    check_schema(s);
}
