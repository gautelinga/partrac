#include "schema_checks.hpp"

#include "filaments/filaments_schema.hpp"
#include "interpol/interpol_schema.hpp"
#include "partrac/partrac_schema.hpp"
#include "spatial/static_space_stepper_schema.hpp"

#include "tracers/tracers_schema.hpp"
#include "tracertensors/tracertensors_schema.hpp"
#include "tracervectors/tracervectors_schema.hpp"
#include "tracervectors_spatial/tracervectors_spatial_schema.hpp"
#include "weighted_walkers/weighted_walkers_schema.hpp"

// One translation unit for every app; the namespace keeps the two families apart
TEST_CASE("every app declares a usable schema", "[schema]") {
  for (auto s : {filaments_schema(), interpol_schema(),
                 partrac_schema(), spatial_schema(),
                 tracers_schema(),
                 tracertensors_schema(),
                 tracervectors_schema(),
                 tracervectors_spatial_schema(),
                 weighted_walkers_schema()})
    check_schema(s);
}
