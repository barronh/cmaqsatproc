{
    "datagrp": "HDFEOS/SWATHS/ColumnAmountNO2/Data Fields",
    "datakeys": [
        "ScatteringWeight", "ScatteringWtPressure", "AmfTrop", "AmfStrat",
        "ColumnAmountNO2", "ColumnAmountNO2Trop", "ColumnAmountNO2Strat",
        "VcdApTrop", "VcdApStrat", "VcdApBelowCloud",
        "TropopausePressure", "TerrainPressure", "CloudFraction",
        "XTrackQualityFlags", "VcdQualityFlags"
    ],
    "datadims": {
        "phony_dim_0": "nTimes",
        "phony_dim_1": "nXtrack",
        "phony_dim_2": "nLevels"
    },
    "geogrp": "HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields",
    "geokeys": ["Time", "SolarZenithAngle", "Latitude", "Longitude"],
    "geodims": {
        "phony_dim_7": "nTimes",
        "phony_dim_5": "nXtrack"
    },
    "pressurekey": "ScatteringWtPressure",
    "level_center_dim": "nLevels",
    "level_edge_dim": "nPresLevels",
    "grndfilterexpr": "(SolarZenithAngle >= 70)",
    "datafilterexpr": (
        "(CloudFraction[:] >= (300 if CloudFraction.dtype.char == 'h' else 0.30)) | " +
        "(np.bitwise_and(XTrackQualityFlags[:].filled(0), 1) == 1) | " +
        "(np.bitwise_and(VcdQualityFlags, 1) == 1)"
    ),
    "renamevars": {
        "ScatteringWeight": "ScatWgt",
        "ColumnAmountNO2": "VCD_NO2_Total",
        "ColumnAmountNO2Trop": "VCD_NO2_Trop",
        "ColumnAmountNO2Strat": "VCD_NO2_Strat",
        "TropopausePressure": "TropopausePres",
        "XTrackQualityFlags": "XTrackQualFlag",
        "VcdQualityFlags": "VcdQualFlags",
        "VcdApBelowCloud": "VcdApBelowCld"
    }
}
