{
    "datagrp": "",
    "datakeys": [
        "O3_column_L2gpValue", "O3_column_Quality", "O3_column_Status",
    ],
    "geogrp": "",
    "geokeys": [
        "O3_column_Time", "O3_column_SolarZenithAngle",
        "O3_column_Latitude", "O3_column_Longitude"
    ],
    "opendapdims": {
        "O3_column_Latitude": "nTimes",
    },
    "datadims": {
    },
    "geodims": {
    },
    "Time": "O3_column_Time",
    "Latitude": "O3_column_Latitude",
    "Longitude": "O3_column_Longitude",
    "pressurekey": None,
    "grndfilterexpr": "(O3_column_SolarZenithAngle >= 70)",
    "datafilterexpr": (
        "(O3_column_Status[:] != 0)"
    ),
    # Reduce name lengths to 14 or lower
    # IOAPI can be 16, but 15 is better
    # N<name> is one more
    "renamevars": {
        "O3_column_L2gpValue": "VCD_O3_Strat",
        "O3_column_Quality": "VCD_O3_Quality",
        "O3_column_Status": "VCD_O3_Status",
    }
}
