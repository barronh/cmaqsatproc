{
    "datagrp": "HDFEOS/SWATHS/OMI Vertical Ozone Profile/Data Fields",
    "datakeys": [
        "O3TroposphericColumn", "O3APrioriProfile", "O3RetrievedProfile",
        "O3AveragingKernel", "TropopauseIndex", "ProfileLevelPressure", "RMS",
        "ExitStatus", "AverageResiduals", "EffectiveCloudFraction"
    ],
    "geogrp": "HDFEOS/SWATHS/OMI Vertical Ozone Profile/Geolocation Fields",
    "geokeys": ["Time", "SolarZenithAngle", "Latitude", "Longitude"],
    "datadims": {
        "phony_dim_0": "nTimes",
        "phony_dim_1": "nXtrack",
        "phony_dim_4": "nLevels",
        "phony_dim_7": "nLevelEdges"
    },
    "geodims": {
        "phony_dim_9": "nTimes",
        "phony_dim_10": "nXtrack"
    },
    "pressurekey": "ProfileLevelPressure",
    # Data filtering is unclear.
    #
    # According to the README file[1]: "The Exit Status (ES) is one of the most
    # important variables needed to assess the data quality. It is recommended
    # to use 0 < ES < 10. Other useful parameters that can be used to assess
    # the retrieval quality include the fitting RMS, fitting residuals, cloud
    # fraction, aerosol index, and glint probability. It is recommended to use
    # retrievals with RMS <2.0 and average fitting residuals <2% (filter ~4% of
    # the retrievals for SZA < 75degrees and 8% of the retrievals for larger
    # SZA)."
    #
    # However, the postprocessing tools online[2] use a maxrms value of 3 or 4
    # and maxavgres of 3 or 4.
    #
    # Due to data sparsity after 2010 when using 2, we opt for 3. If your data
    # is still too sparse, consider modfifying the maxrms to 4.
    #
    # [1] https://avdc.gsfc.nasa.gov/pub/data/satellite/Aura/OMI/V03/L2/
    #     OMPROFOZ/OMPROFOZ_readme-v3.pdf
    # [2] https://avdc.gsfc.nasa.gov/pub/data/satellite/Aura/OMI/V03/L2/
    #      OMPROFOZ/programs/extract_omi_PROFOZhe5.pro
    "grndfilterexpr": "(SolarZenithAngle >= 70)",
    "datafilterexpr": (
        "(ExitStatus[:] <= 0) | (ExitStatus[:] >= 10) | " +
        "(RMS[:].max(-1) > 3) | (AverageResiduals[:].max(-1) >= 3) | "
        "(EffectiveCloudFraction >= 0.3)"
    ),
    # Reduce name lengths to 14 or lower
    # IOAPI can be 16, but 15 is better
    # N<name> is one more
    "renamevars": {
        "O3TroposphericColumn": "VCD_O3_Trop",
        "O3APrioriProfile": "O3APriori",
        "O3RetrievedProfile": "O3Retrieved",
        "O3AveragingKernel": "AvgKernel",
        "TropopauseIndex": "TropopauseIdx",
        "ProfileLevelPressure": "ProfilePress",
        "AverageResiduals": "AvgResiduals",
        "EffectiveCloudFraction": "EffectCldFrac"
    }
}
