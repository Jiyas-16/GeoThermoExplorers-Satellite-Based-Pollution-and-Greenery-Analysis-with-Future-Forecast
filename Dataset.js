var aoi = ee.Geometry.Polygon(
    [[[77.348709, 28.842898],
      [77.348709, 28.412356],
      [76.935635, 28.412356],
      [76.935635, 28.842898]]]);

var startDate = '2022-05-01';
var endDate = '2022-12-31';


// Function to apply scaling factors
function applyScaleFactors(image) {
    var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
    var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
    return image.addBands(opticalBands, null, true)
                .addBands(thermalBands, null, true);
}

// Function to perform cloud mask
function maskL8sr(col) {
    var cloudShadowBitMask = (1 << 3);
    var cloudsBitMask = (1 << 5);
    var qa = col.select('QA_PIXEL');
    var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
             .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
    return col.updateMask(mask);
}

// Load LANDSAT image collection
var image = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
            .filterDate(startDate, endDate)
            .filterBounds(aoi)
            .map(applyScaleFactors)
            .map(maskL8sr)
            .median();

// Visualization parameters for true color image
var visualization = {
    bands: ['SR_B4', 'SR_B3', 'SR_B2'],
    min: 0.0,
    max: 0.3,
};

// Add true color image layer to the map
Map.addLayer(image, visualization, 'True Color (432)', false);

// Calculate and visualize NDVI
var ndvi  = image.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');
Map.addLayer(ndvi, {min:-1, max:1, palette: ['blue', 'white', 'green']}, 'NDVI', false);

// Calculate minimum and maximum NDVI values
var ndvi_min = ee.Number(ndvi.reduceRegion({
    reducer: ee.Reducer.min(),
    geometry: aoi,
    scale: 30,
    maxPixels: 1e9
}).values().get(0));

var ndvi_max = ee.Number(ndvi.reduceRegion({
    reducer: ee.Reducer.max(),
    geometry: aoi,
    scale: 30,
    maxPixels: 1e9
}).values().get(0));

// Calculate vegetation fraction and emissivity
var fv = (ndvi.subtract(ndvi_min).divide(ndvi_max.subtract(ndvi_min))).pow(ee.Number(2)).rename('FV');
var em = fv.multiply(ee.Number(0.004)).add(ee.Number(0.986)).rename('EM');

// Select and rename thermal band for land surface temperature calculation
var thermal = image.select('ST_B10').rename('thermal');

// Calculate land surface temperature
var lst = thermal.expression(
    '(tb / (1 + (0.00115 * (tb/0.48359547432)) * log(em))) - 273.15',
    {'tb':thermal.select('thermal'),
    'em': em}).rename('LST');

// Visualization parameters for land surface temperature
var lst_vis = {
    min: 25,
    max: 50,
    palette: [
        '040274', '040281', '0502a3', '0502b8', '0502ce', '0502e6',
        '0602ff', '235cb1', '307ef3', '269db1', '30c8e2', '32d3ef',
        '3be285', '3ff38f', '86e26f', '3ae237', 'b5e22e', 'd6e21f',
        'fff705', 'ffd611', 'ffb613', 'ff8b13', 'ff6e08', 'ff500d',
        'ff0000', 'de0101', 'c21301', 'a71001', '911003'
    ]
};

// Add land surface temperature layer to the map
Map.addLayer(lst, lst_vis, 'LST AOI');

// Calculate mean and standard deviation of land surface temperature in AOI
var lst_mean = ee.Number(lst.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: aoi,
    scale: 30,
    maxPixels: 1e9
}).values().get(0));

var lst_std = ee.Number(lst.reduceRegion({
    reducer: ee.Reducer.stdDev(),
    geometry: aoi,
    scale: 30,
    maxPixels: 1e9
}).values().get(0));

print('Mean LST in AOI:', lst_mean);
print('STD LST in AOI:', lst_std);

// Calculate and visualize Urban Heat Island (UHI)
var uhi = lst.subtract(lst_mean).divide(lst_std).rename('UHI');
var uhi_vis = {
    min: -4,
    max: 4,
    palette:['313695', '74add1', 'fed976', 'feb24c', 'fd8d3c', 'fc4e2a', 'e31a1c', 'b10026']
};
Map.addLayer(uhi, uhi_vis, 'UHI AOI');

// Calculate and visualize Urban Thermal Field Variance Index (UTFVI)
var utfvi = lst.subtract(lst_mean).divide(lst).rename('UTFVI');
var utfvi_vis = {
    min: -1,
    max: 0.3,
    palette:['313695', '74add1', 'fed976', 'feb24c', 'fd8d3c', 'fc4e2a', 'e31a1c', 'b10026']
};
Map.addLayer(utfvi, utfvi_vis, 'UTFVI AOI');



// Load Sentinel-5P NO2 image collection
var collection = ee.ImageCollection("COPERNICUS/S5P/NRTI/L3_NO2")
                .filterDate(startDate, endDate)
                .filterBounds(aoi);

// Select NO2 band
var no2 = collection.select('NO2_column_number_density');

// Visualization parameters for NO2
var no2_vis = {
    min: 0,
    max: 0.0002,
    palette: ['black', 'blue', 'purple', 'cyan', 'green', 'yellow', 'red']
};

// Add NO2 layer to the map
Map.addLayer(no2.mean(), no2_vis, 'NO2 Pollution');

// Calculate mean and standard deviation of NO2 pollution in AOI
var no2_mean = ee.Number(no2.mean().reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: aoi,
    scale: 1000, // You may need to adjust the scale
    maxPixels: 1e9
}).values().get(0));

var no2_std = ee.Number(no2.mean().reduceRegion({
    reducer: ee.Reducer.stdDev(),
    geometry: aoi,
    scale: 1000, // You may need to adjust the scale
    maxPixels: 1e9
}).values().get(0));

print('Mean NO2 in AOI:', no2_mean);
print('STD NO2 in AOI:', no2_std);



// Load the GIMMS NDVI ImageCollection
var dataset = ee.ImageCollection('NASA/GIMMS/3GV0')
              .filter(ee.Filter.date('2013-06-01', '2013-12-31'));

// Calculate the NDVI (Normalized Difference Vegetation Index)
var ndvi = dataset.select('ndvi');

// Visualization parameters for NDVI
var ndvi_vis = {
  min: -1,
  max: 1,
  palette: ['blue', 'white', 'green']
};

// Add NDVI layer to the map
Map.addLayer(ndvi.mean(), ndvi_vis, 'NDVI');

// Calculate mean and standard deviation of NDVI in AOI
var ndvi_mean = ee.Number(ndvi.mean().reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: aoi,
    scale: 1000, // You may need to adjust the scale
    maxPixels: 1e9
}).values().get(0));

var ndvi_std = ee.Number(ndvi.mean().reduceRegion({
    reducer: ee.Reducer.stdDev(),
    geometry: aoi,
    scale: 1000, // You may need to adjust the scale
    maxPixels: 1e9
}).values().get(0));

print('Mean NDVI in AOI:', ndvi_mean);
print('STD NDVI in AOI:', ndvi_std);





// Load the GPW Population Density ImageCollection
var populationDensity = ee.ImageCollection("CIESIN/GPWv411/GPW_Population_Density");

// Visualization parameters for population density
var populationDensity_vis = {
  min: 0,
  max: 1000, // You may need to adjust the maximum value based on your data
  palette: ['black', 'blue', 'purple', 'cyan', 'green', 'yellow', 'red']
};

// Add population density layer to the map
Map.addLayer(populationDensity.mean(), populationDensity_vis, 'Population Density');

// Calculate mean and standard deviation of population density in AOI
var populationDensity_mean = ee.Number(populationDensity.mean().reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: aoi,
    scale: 1000, // You may need to adjust the scale
    maxPixels: 1e9
}).values().get(0));

var populationDensity_std = ee.Number(populationDensity.mean().reduceRegion({
    reducer: ee.Reducer.stdDev(),
    geometry: aoi,
    scale: 1000, // You may need to adjust the scale
    maxPixels: 1e9
}).values().get(0));

print('Mean Population Density in AOI:', populationDensity_mean);
print('STD Population Density in AOI:', populationDensity_std);






// Load the LCZ ImageCollection
var urbanization = ee.ImageCollection('RUB/RUBCLIM/LCZ/global_lcz_map/latest');

// Visualization parameters for urbanization
var urbanization_vis = {
  min: 1,
  max: 17,
  palette: [
    'aec3d4', '152106', '225129', '369b47', '30eb5b', '387242',
    '6a2325', 'c3aa69', 'b76031', 'd9903d', '91af40', '111149',
    'cdb33b', 'cc0013', '33280d', 'd7cdcc', 'f7e084'
  ]
};

// Reduce the image collection to a single image by taking the median
var urbanization_median = urbanization.median();

// Select a single band for visualization
var urbanization_single_band = urbanization_median.select('LCZ_Filter');

// Add urbanization layer to the map
Map.addLayer(urbanization_single_band, urbanization_vis, 'Urbanization');

// Calculate mean and standard deviation of urbanization in AOI
var urbanization_mean = ee.Number(urbanization_single_band.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: aoi,
    scale: 1000, // You may need to adjust the scale
    maxPixels: 1e9
}).values().get(0));

var urbanization_std = ee.Number(urbanization_single_band.reduceRegion({
    reducer: ee.Reducer.stdDev(),
    geometry: aoi,
    scale: 1000, // You may need to adjust the scale
    maxPixels: 1e9
}).values().get(0));

print('Mean Urbanization in AOI:', urbanization_mean);
print('STD Urbanization in AOI:', urbanization_std);







// Load the CO image collection
var collection = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_CO')
                  .filterDate(startDate, endDate)
                  .filterBounds(aoi);

// Select the CO_column_number_density band
var co = collection.select('CO_column_number_density');

// Calculate mean and standard deviation of CO concentration in the area of interest
var meanCO = co.mean().reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: aoi,
    scale: 1113.2,  // Resolution in meters
    maxPixels: 1e9
});

var stdCO = co.reduce(ee.Reducer.stdDev()).reduceRegion({
    reducer: ee.Reducer.stdDev(),
    geometry: aoi,
    scale: 1113.2,  // Resolution in meters
    maxPixels: 1e9
});

// Print mean and standard deviation of CO concentration
print('Mean CO concentration in AOI:', meanCO);
print('STD CO concentration in AOI:', stdCO);

// Visualization parameters for CO
var coVis = {
    min: 0,  // Adjust min value according to your data
    max: 0.1,  // Adjust max value according to your data
    palette: ['black', 'blue', 'green', 'yellow', 'red'] // Adjust the colors as needed
};

// Add CO layer to the map
Map.addLayer(co.mean(), coVis, 'CO Concentration');


// Calculate mean and standard deviation of land surface temperature in AOI
var lst_mean = ee.Number(lst.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: aoi,
    scale: 30,
    maxPixels: 1e9
}).values().get(0));

var lst_std = ee.Number(lst.reduceRegion({
    reducer: ee.Reducer.stdDev(),
    geometry: aoi,
    scale: 30,
    maxPixels: 1e9
}).values().get(0));

// Create a chart for land surface temperature

// Create a chart for NO2 pollution
var no2_chart = ui.Chart.image.seriesByRegion({
    imageCollection: no2,
    regions: aoi,
    reducer: ee.Reducer.mean(),
    scale: 1000,
    xProperty: 'system:time_start'
  })
  .setChartType('LineChart')
  .setOptions({
    title: 'NO2 Pollution Over Time',
    hAxis: {title: 'Date'},
    vAxis: {title: 'NO2 Pollution'},
    lineWidth: 1,
    pointSize: 3,
    series: {0: {color: '0000FF'}},
  });

// Display the chart
print(no2_chart);

// Create a chart for NDVI
var ndvi_chart = ui.Chart.image.seriesByRegion({
    imageCollection: ndvi,
    regions: aoi,
    reducer: ee.Reducer.mean(),
    scale: 1000,
    xProperty: 'system:time_start'
  })
  .setChartType('LineChart')
  .setOptions({
    title: 'NDVI Over Time',
    hAxis: {title: 'Date'},
    vAxis: {title: 'NDVI'},
    lineWidth: 1,
    pointSize: 3,
    series: {0: {color: '00FF00'}},
  });

// Display the chart
print(ndvi_chart);

// Create a chart for Population Density
var populationDensity_chart = ui.Chart.image.seriesByRegion({
    imageCollection: populationDensity,
    regions: aoi,
    reducer: ee.Reducer.mean(),
    scale: 1000,
    xProperty: 'system:time_start'
  })
  .setChartType('LineChart')
  .setOptions({
    title: 'Population Density Over Time',
    hAxis: {title: 'Date'},
    vAxis: {title: 'Population Density'},
    lineWidth: 1,
    pointSize: 3,
    series: {0: {color: 'FFFF00'}},
  });

// Display the chart
print(populationDensity_chart);










var urbanization = ee.ImageCollection('RUB/RUBCLIM/LCZ/global_lcz_map/latest');

// Reduce the image collection to a single image by taking the median
var urbanization_median = urbanization.median();

// Select a single band for visualization
var urbanization_single_band = urbanization_median.select('LCZ');

// Get the histogram of urbanization classes
var histogram = ui.Chart.image.histogram({
  image: urbanization_single_band,
  region: aoi,
  scale: 40, // Adjust the scale as needed
  maxBuckets: 17
}).setOptions({
  title: 'Urbanization Histogram',
  hAxis: {title: 'Urbanization Class'},
  vAxis: {title: 'Frequency'}
});

// Display the histogram
print(histogram);




// Create a chart for CO concentration
var co_chart = ui.Chart.image.seriesByRegion({
    imageCollection: collection,
    regions: aoi,
    reducer: ee.Reducer.mean(),
    scale: 1113.2,
    xProperty: 'system:time_start'
  })
  .setChartType('LineChart')
  .setOptions({
    title: 'CO Concentration Over Time',
    hAxis: {title: 'Date'},
    vAxis: {title: 'CO Concentration'},
    lineWidth: 1,
    pointSize: 3,
    series: {0: {color: '000000'}},
  });

// Display the chart
print(co_chart);


// Load LANDSAT image collection
var imageCollection = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
            .filterDate(startDate, endDate)
            .filterBounds(aoi)
            .map(applyScaleFactors)
            .map(maskL8sr);

// Create a chart for mean LST over time
var lstTimeSeries = imageCollection.map(function(image) {
    var lst = image.select('ST_B10').rename('thermal');
    var lstStats = lst.reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: aoi,
        scale: 30,
        maxPixels: 1e9
    });
    
    // Return feature with time start and temperature
    return ee.Feature(null, lstStats)
            .set('system:time_start', image.get('system:time_start'));
}).filter(ee.Filter.notNull(['thermal'])); // Filter out null features

// Convert temperature from Kelvin to Celsius
lstTimeSeries = lstTimeSeries.map(function(feature) {
    var kelvin = ee.Number(feature.get('thermal'));
    var celsius = kelvin.subtract(273.15);
    return feature.set('thermal_celsius', celsius);
});

// Create a chart for mean LST over time
var meanLSTChart = ui.Chart.feature.byFeature(lstTimeSeries, 'system:time_start', 'thermal_celsius')
    .setChartType('LineChart')
    .setOptions({
        title: 'Mean Land Surface Temperature (LST) Over Time',
        hAxis: {title: 'Time'},
        vAxis: {title: 'Temperature (°C)', viewWindow: {min: 0, max: 80}},
        lineWidth: 1,
        pointSize: 3
});

// Display the chart
print(meanLSTChart);
