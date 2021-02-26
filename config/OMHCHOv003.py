{
    "datagrp": "HDFEOS/SWATHS/OMI Total Column Amount HCHO/Data Fields",
    "datakeys": [
        "ScatteringWeights", "AirMassFactor",
        "ColumnAmount", "ColumnAmountDestriped", 
        "AverageFittingRMS", "FittingRMS", "GasProfile",
        "ReferenceSectorCorrectedVerticalColumn",
        "MainDataQualityFlag", "ClimatologyLevels",
        "AMFCloudFraction"
    ],
    "datadims": {
        "phony_dim_0": "nTimes",
        "phony_dim_1": "nXtrack",
        "phony_dim_3": "nLevels"
    },
    "geogrp": "HDFEOS/SWATHS/OMI Total Column Amount HCHO/Geolocation Fields",
    "geokeys": [
        "Time", "SolarZenithAngle", "Latitude", "Longitude",
        "XtrackQualityFlags", 
    ],
    "geodims": {
        "phony_dim_7": "nTimes",
        "phony_dim_8": "nXtrack"
    },
    "pressurekey": "ClimatologyLevels",
    "grndfilterexpr": "(SolarZenithAngle >= 70)",
    "datafilterexpr": (
        "(AMFCloudFraction >= 0.3) | " +
        "(np.bitwise_and(XtrackQualityFlags[:].filled(0), 1) == 1) | " +
        "(MainDataQualityFlag[:] != 0)"
    ),
    "renamevars": {
        "ScatteringWeights": "ScatWgt",
        "ColumnAmount": "VCD_HCHO_Total",
        "ColumnAmountDestriped": "VCD_HCHO_Dstrp",
        "AverageFittingRMS": "AvgFitRMS",
        "ReferenceSectorCorrectedVerticalColumn": "VCD_HCHO_RefCr",
        "XtrackQualityFlags": "XTrackQualFlag",
        "MainDataQualityFlag": "MainQualFlag",
        "ClimatologyLevels": "ClimLevels",
        "AMFCloudFraction": "AMFCloudFrac"
    }
}
