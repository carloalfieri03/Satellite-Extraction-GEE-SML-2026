// radius
var radius = 30; 
var bufferedTable = table.map(function(feature) {
  return feature.buffer(radius);
});

// red dots are cohordinates points
Map.addLayer(table, {color: 'red'}, 'Transect Points', true);

// in blue the 30m buffers to check extraction zone
Map.addLayer(bufferedTable, {color: 'blue'}, '30m Extraction Buffers', true);


var aoi = table.geometry().bounds(); 
Map.addLayer(aoi, {}, 'Study Area', false);
Map.centerObject(aoi, 10);


// masking for sentinel 2 using scene classification layer

function cloudMask(image){
  var scl = image.select('SCL');
  var mask = scl.eq(3).or(scl.gte(7).and(scl.lte(10))); // cloud shadow(3) or between 7-10 (uncl,medium,high proba, cirrus)
  var opticalScaled = image.select('B.*').divide(10000); // digital numbers conversion
  return image.updateMask(mask.eq(0)) // applies mask and keeps only 0 valued pixels
              .addBands(opticalScaled, null, true) // scales image and applies
              .copyProperties(image, ["system:time_start"]); // keeps date and time metadata of the image
}



// Import Sentinel-2 Surface Reflectance Collection

var sentinelimage = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
  .filterBounds(aoi)
  .filterDate('2026-03-01', '2026-05-20')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 10)) // max cloud percentage
  .map(cloudMask) // apply cleaning function
  .mean()
  .clip(aoi);


//  Indexes

// Ndvi, savi, ndwi , ndbi
var ndvi = sentinelimage.normalizedDifference(['B8', 'B4']).rename('NDVI');

var savi = sentinelimage.expression('1.5*((NIR-RED)/(NIR+RED+0.5))',
  {'NIR':sentinelimage.select('B8'),
  'RED':sentinelimage.select('B4')}).rename('SAVI') // soil brightness correction of 0.5  ( 1.0 for low vegetation and 0 for area with dense vegetation)

var ndwi = sentinelimage.normalizedDifference(['B3', 'B8']).rename('NDWI');

var ndbi = sentinelimage.normalizedDifference(['B11', 'B8']).rename('NDBI');
// we know that sentinel swir band is at 20 m resolution but then this will be upsampled

Map.addLayer(ndvi.clip(aoi), {min:-1, max:1, palette:['blue','white','green']}, 'NDVI', false);
Map.addLayer(ndbi.clip(aoi), {min:-1, max:1, palette:['green','white','red']}, 'NDBI', false);

var vegfrac = sentinelimage.expression(
    '(NDVI - NDVI_min) / (NDVI_max - NDVI_min)',
    {
      'NDVI': ndvi,    
      'NDVI_min': 0.15, 
      'NDVI_max': 0.80  
    }
  )
  .clamp(0, 1)         // normalize percentage
  .rename('Veg_Frac')
  .toFloat();
  
var buildfrac = sentinelimage.expression(
    '(NDBI - NDBI_min) / (NDBI_max - NDBI_min)',
    {
      'NDBI': ndbi,       
      'NDBI_min': -0.20,  
      'NDBI_max': 0.30    
    }
  )
  .clamp(0, 1)            
  .rename('Build_Frac')
  .toFloat();

//  google dynamic since most recent data available

var dwld = ee.ImageCollection("GOOGLE/DYNAMICWORLD/V1")
  .filterBounds(aoi)
  .filterDate('2026-03-01', '2026-05-20');
  

// mode for the pixel value
var dynimg = dwld.select('label').mode();


// mask the image according to water that we want to extract
var waterMask = dynimg.eq(0);

// distance from water 
// 20 km distance set
var distWater = waterMask.not().fastDistanceTransform(2000).sqrt()
                         .multiply(10) // * 10 to convert pixels to meters
                         .rename('Dist_Water');


// topographic values are at 30m so will be resampled

//  Copernicus Global instead of Nasa since Copernicus is more recent, 2011-2015.
var dem = ee.ImageCollection('COPERNICUS/DEM/GLO30')
  .select('DEM')
  .mosaic()
  .resample('bilinear') // upsampling
  .rename('Elevation');
var slope = ee.Terrain.slope(dem).rename('Slope');
var aspect = ee.Terrain.aspect(dem).rename('Aspect');

// tpi with 100 m radius to calculate topographic index
var focalMean = dem.focal_mean(100, 'circle', 'meters');
var tpi = dem.subtract(focalMean).rename('TPI');

Map.addLayer(dem.clip(aoi), {min:0, max:300, palette:['blue','green','yellow','brown','white']}, 'Elevation', false);

// for images

var Hillshade= ee.Terrain.hillshade(dem);

// albedo

var albedo = sentinelimage.expression(
'((0.356 * B2) + (0.130 * B4) + (0.373 * B5) + (0.085 * B6) + (0.072 * B7) - 0.0018) / 1.016',
    {
      'B2':  sentinelimage.select('B2'),  // Blue (10m)
      'B4':  sentinelimage.select('B4'),  // Red (10m)
      'B8':  sentinelimage.select('B8'),  // NIR (10m)
      'B11': sentinelimage.select('B11'), // SWIR 1 (20m -> auto-scaled to 10m on export)
      'B12': sentinelimage.select('B12')  // SWIR 2 (20m -> auto-scaled to 10m on export)
    }
  )
  .clamp(0, 1)        
  .rename('Albedo')
  .toFloat();

// LST 

//function to do additional filtering and resampling at 10m with bilinear sampling.
function prepLandsatLST(image) {
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 3).eq(0).and(qa.bitwiseAnd(1 << 4).eq(0));
  var lstCelsius = image.select('ST_B10')
                        .multiply(0.00341802).add(149.0)
                        .subtract(273.15)
                        .rename('LST_10m');
                        
  return lstCelsius.updateMask(mask).resample('bilinear');
}


var lst10m = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
  .filterBounds(aoi)
  .filterDate('2026-03-01', '2026-05-20')
  .filter(ee.Filter.lt('CLOUD_COVER', 10))
  .map(prepLandsatLST) 
  .mean();
  



Map.addLayer(lst10m.clip(aoi), {min: 15, max: 40, palette: ['blue', 'yellow', 'red']}, 'LST 10m', false);

//  csv export

var finalStack = ee.Image([
  ndvi, savi, ndwi, ndbi,       
  vegfrac, buildfrac, 
  distWater,                   
  dem, slope, aspect, tpi,albedo,lst10m     
]).toFloat();

var extractedData = finalStack.reduceRegions({
  collection: bufferedTable, // buffered areas defined at the top
  reducer: ee.Reducer.mean(), 
  scale: 10                   // Forces extraction at 10m resolution
});


Export.table.toDrive({
  collection: extractedData,
  description: 'Rome_trainingSentinel_30radius_CSV',
  folder: 'GEE_Exports_10m', 
  fileFormat: 'CSV'
});


// Images exports

//  Spectral Indices 
Export.image.toDrive({
  image: ee.Image([ndvi, savi, ndwi, ndbi]).toFloat(),
  description: 'S2_Spectral_Indices_30rad',
  folder: 'GEE_Exports_10m',
  region: aoi,
  scale: 10,
  maxPixels: 1e13
});

//  Topography 
Export.image.toDrive({
  image: ee.Image([dem, slope, aspect, tpi]).toFloat(),
  description: 'Topography_30rad',
  folder: 'GEE_Exports_10m',
  region: aoi,
  scale: 10,
  maxPixels: 1e13
});

// Land coverage metrics 

Export.image.toDrive({
  image: ee.Image([vegfrac, buildfrac, distWater]).toFloat(),
  description: 'LandCover_Metrics_30rad',
  folder: 'GEE_Exports_10m',
  region: aoi,
  scale: 10,
  maxPixels: 1e13
});


// LST

Export.image.toDrive({
  image: lst10m.toFloat(),
  description: 'LST_10m_30rad',
  folder: 'GEE_Exports_10m',
  region: aoi,
  scale: 10,
  maxPixels: 1e13
});
