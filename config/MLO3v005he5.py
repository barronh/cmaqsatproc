{
    "datagrp": "HDFEOS/SWATHS/O3 column/Data Fields",
    "datakeys": [
        "L2gpValue", "Quality", "Status",
    ],
    "geogrp": "HDFEOS/SWATHS/O3 column/Geolocation Fields",
    "geokeys": ["Time", "SolarZenithAngle", "Latitude", "Longitude"],
    "datadims": {
    },
    "geodims": {
    },
    "Time": "Time",
    "Latitude": "Latitude",
    "Longitude": "Longitude",
    "pressurekey": None,
    "grndfilterexpr": "(SolarZenithAngle >= 70)",
    "datafilterexpr": (
        "(Status[:] != 0)"
    ),
    # Reduce name lengths to 14 or lower
    # IOAPI can be 16, but 15 is better
    # N<name> is one more
    "renamevars": {
        "L2gpValue": "VCD_O3_Strat",
        "Quality": "VCD_O3_Quality",
        "Status": "VCD_O3_Status",
    }
}
