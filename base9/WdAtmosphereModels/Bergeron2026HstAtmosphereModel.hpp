#ifndef BERG2026HSTATMOS_HPP
#define BERG2026HSTATMOS_HPP

#include <array>
#include <map>
#include <string>
#include <utility>

#include "constants.hpp"
#include "WdAtmosphereModel.hpp"

class Bergeron2026HstAtmosphereModel : public BergeronAtmosphereModel
{
  public:
    Bergeron2026HstAtmosphereModel()
    {
        dirName = "bergeron_2026_hst/";

        availableFilters = { "F435W", "F555W", "F606W", "F814W" }; // HST
        hasBc = false;

        files = {
            {"Table_Mass_0.2"},
            {"Table_Mass_0.3"},
            {"Table_Mass_0.4"},
            {"Table_Mass_0.5"},
            {"Table_Mass_0.6"},
            {"Table_Mass_0.7"},
            {"Table_Mass_0.8"},
            {"Table_Mass_0.9"},
            {"Table_Mass_1.0"},
            {"Table_Mass_1.1"},
            {"Table_Mass_1.2"},
            {"Table_Mass_1.3"}
        };
    }
};

#endif
