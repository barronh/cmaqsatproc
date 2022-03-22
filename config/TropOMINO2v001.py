{
    "datagrp": "PRODUCT",
    "datakeys": [
        "qa_value", "nitrogendioxide_tropospheric_column", "air_mass_factor_total",
        "air_mass_factor_troposphere", "time", "tm5_constant_a", "latitude", "longitude"
    ],
    "Latitude": "latitude",
    "Longitude": "longitude",
    "datadims": {
        "scanline": "nTimes",
        "ground_pixel": "nXtrack",
        "layer": "nPresLevels"
    },
    "geogrp": "PRODUCT",
    "geokeys": ["time", "latitude", "longitude"],
    "geodims": {
        "scanline": "nTimes",
        "ground_pixel": "nXtrack"
    },
    "grndfilterexpr": "(qa_value > 0.75)",
    "datafilterexpr": (
        "(qa_value > 0.75)"
    ),
    "renamevars": {
       "nitrogendioxide_tropospheric_column": "TropVCD",
       "air_mass_factor_total": "AmfTot",
       "air_mass_factor_troposphere": "AmfTrop",
    }
}
