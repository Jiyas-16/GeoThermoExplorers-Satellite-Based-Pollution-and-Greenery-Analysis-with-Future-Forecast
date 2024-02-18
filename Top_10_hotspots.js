var aoi = ee.Geometry.Polygon(
    [[[77.2090, 28.6139],
      [77.2090, 28.2000],
      [77.7500, 28.2000],
      [77.7500, 28.6139]]]);

var startDate = '2022-05-01';
var endDate = '2022-12-31';

var landsatCollection = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterDate(startDate, endDate)
    .filterBounds(aoi);

function applyScaleFactors(image) {
    var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
    var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
    return image.addBands(opticalBands, null, true)
                .addBands(thermalBands, null, true);
}

function maskL8sr(col) {
    var cloudShadowBitMask = (1 << 3);
    var cloudsBitMask = (1 << 5);
    var qa = col.select('QA_PIXEL');
    var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
             .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
    return col.updateMask(mask);
}

var image = landsatCollection.map(applyScaleFactors)
                             .map(maskL8sr)
                             .median();

var ndvi  = image.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');

var no2_collection = ee.ImageCollection("COPERNICUS/S5P/NRTI/L3_NO2")
                .filterDate(startDate, endDate)
                .filterBounds(aoi);

var no2 = no2_collection.select('NO2_column_number_density');

var ndvi_collection = ee.ImageCollection('NASA/GIMMS/3GV0')
              .filter(ee.Filter.date('2013-06-01', '2013-12-31'));

var ndvi_gimms = ndvi_collection.select('ndvi');

var populationDensity_collection = ee.ImageCollection("CIESIN/GPWv411/GPW_Population_Density");

var urbanization_collection = ee.ImageCollection('RUB/RUBCLIM/LCZ/global_lcz_map/latest');

var co_collection = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_CO')
                  .filterDate(startDate, endDate)
                  .filterBounds(aoi);

var randomPoints = ee.FeatureCollection.randomPoints({
  region: aoi,
  points: 70
});

var calculateStats = function(feature) {
  var point = ee.Geometry.Point(feature.geometry().coordinates());
  var aoi = point.buffer(1000);
  
  var lst_mean = ee.Number(image.select('ST_B10').reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: aoi,
    scale: 30,
    maxPixels: 1e9
  }).values().get(0));

  var no2_mean = ee.Number(no2.mean().reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: aoi,
    scale: 1000,
    maxPixels: 1e9
  }).values().get(0));

  var ndvi_mean = ee.Number(ndvi_gimms.mean().reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: aoi,
    scale: 1000,
    maxPixels: 1e9
  }).values().get(0));

  var populationDensity_mean = ee.Number(populationDensity_collection.mean().reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: aoi,
    scale: 1000,
    maxPixels: 1e9
  }).values().get(0));

  var urbanization_mean = ee.Number(urbanization_collection.median().reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: aoi,
    scale: 1000,
    maxPixels: 1e9
  }).values().get(0));

  var meanCO = ee.Number(co_collection.mean().reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: aoi,
    scale: 1113.2,
    maxPixels: 1e9
  }).values().get(0));
  
  return feature.set({
    'LST_mean': lst_mean,
    'NO2_mean': no2_mean,
    'NDVI_mean': ndvi_mean,
    'PopulationDensity_mean': populationDensity_mean,
    'Urbanization_mean': urbanization_mean,
    'CO_mean': meanCO
  });
};

var locationStats = randomPoints.map(calculateStats);

var sortedLocations = locationStats.sort('LST_mean', false)
                                  .sort('NO2_mean', false)
                                  .sort('NDVI_mean', true)
                                  .sort('PopulationDensity_mean', false)
                                  .sort('Urbanization_mean', false)
                                  .sort('CO_mean', false);

var top10Locations = sortedLocations.limit(10);

var features = top10Locations.getInfo().features;

var locationNames = [];

var top10FeatureCollection = ee.FeatureCollection(features.map(function(feature, index) {
  var coordinates = feature.geometry.coordinates;
  var locationName = "Location " + (index + 1) + ": " + coordinates[1] + ", " + coordinates[0];
  locationNames.push(locationName);
  return ee.Feature(ee.Geometry.Point(coordinates), {});
}));

print("List of top 10 locations:");
print(locationNames);

var map = ui.Map();

map.addLayer(top10FeatureCollection, {color: 'red'}, 'Top 10 Locations');

var firstLocation = features[0].geometry.coordinates;
map.setCenter(firstLocation[0], firstLocation[1], 12);

print(map);
