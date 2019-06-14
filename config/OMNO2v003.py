varkeys = ['ScatteringWeight', 'AmfTrop', 'AmfStrat', 'ColumnAmountNO2', 'ColumnAmountNO2Trop', 'ColumnAmountNO2Strat', 'VcdApTrop', 'VcdApStrat', 'VcdApBelowCloud', 'TropopausePressure']
datagrp = 'HDFEOS/SWATHS/ColumnAmountNO2/Data Fields'
datatgtdim = 'phony_dim_0'
geotgtdim = 'phony_dim_7'
geogrp = 'HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields'
pressurekey = 'ScatteringWtPressure'
grndfilterexpr = (
    '(SolarZenithAngle >= 70)'
)
datafilterexpr = (
    '(CloudFraction >= 300) | ' +
    '(np.bitwise_and(XTrackQualityFlags[:].filled(0), 1) == 1) | ' +
    '(np.bitwise_and(VcdQualityFlags, 1) == 1)'
)

