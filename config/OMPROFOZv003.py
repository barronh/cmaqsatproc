varkeys = ['O3APrioriProfile', 'O3RetrievedProfile', 'O3AveragingKernel', 'TropopauseIndex']

datagrp = 'HDFEOS/SWATHS/OMI Vertical Ozone Profile/Data Fields'
geogrp = 'HDFEOS/SWATHS/OMI Vertical Ozone Profile/Geolocation Fields'
geovars = ['SolarZenithAngle', 'Latitude', 'Longitude']

datatgtdim = 'phony_dim_0'
geotgtdim = 'phony_dim_8'
pressurekey = 'ProfileLevelPressure'
grndfilterexpr = (
    '(SolarZenithAngle >= 70)'
)
datafilterexpr = (
    '(RMS[..., 0] >= 2) | ' +
    '(RMS[..., 1] >= 2) | ' +
    '(EffectiveCloudFraction >= 0.3)'
)

