varkeys = ['O3TroposphericColumn', 'O3APrioriProfile', 'O3RetrievedProfile', 'O3AveragingKernel', 'TropopauseIndex']

datagrp = 'HDFEOS/SWATHS/OMI Vertical Ozone Profile/Data Fields'
geogrp = 'HDFEOS/SWATHS/OMI Vertical Ozone Profile/Geolocation Fields'
geovars = ['SolarZenithAngle', 'Latitude', 'Longitude']

datatgtdim = None # 'phony_dim_0'
geotgtdim = None # 'phony_dim_9'
pressurekey = 'ProfileLevelPressure'
grndfilterexpr = (
    '(SolarZenithAngle >= 70)'
)
datafilterexpr = (
    '(RMS[..., 0] >= 2) | ' +
    '(RMS[..., 1] >= 2) | ' +
    '(EffectiveCloudFraction >= 0.3)'
)

