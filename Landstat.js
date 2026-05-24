// Calculate & Export NDVI, SAVI & NDWI Indices from Landsat 8 Imagery in Google Earth Engine

// radius
var radius = 90; 

var bufferedTable = table.map(function(feature) {
  return feature.buffer(radius);
});

// red dots are cohordinates points
Map.addLayer(table, {color: 'red'}, 'Transect Points', true);

// in blue the 30m buffers to check extraction zone
Map.addLayer(bufferedTable, {color: 'blue'}, 'Extraction Buffers', true);


var aoi = table.geometry().bounds(); 
Map.addLayer(aoi, {}, 'Study Area', false);
Map.centerObject(aoi, 10);

// 2. Apply scaling factors to convert digital numbers to physical values
function applyImageScaling(image) {
  // Scale optical bands (SR_B) to reflectance values
  var opticalScaled = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  // Scale thermal bands (ST_B*) to temperature values in Kelvin
  var thermalScaled = image.select('ST_B.*').multiply(0.00341802).add(149.0);
  // Replace original bands with scaled bands
  return image.addBands(opticalScaled, null, true)
              .addBands(thermalScaled, null, true);
}

// 3. Function to mask clouds and cloud shadows using the QA_PIXEL band
function applyCloudMask(imageCollection) {
  // Define bit masks for cloud shadow (bit 3) and cloud (bit 5)
  var cloudShadowMaskBit = (1 << 3);
  var cloudMaskBit = (1 << 4);  // !! to check
  // Select the quality assessment band
  var qualityAssessment = imageCollection.select('QA_PIXEL');
  // Create mask where both cloud and cloud shadow bits are 0 (clear conditions)
  var clearSkyMask = qualityAssessment.bitwiseAnd(cloudShadowMaskBit).eq(0)
                 .and(qualityAssessment.bitwiseAnd(cloudMaskBit).eq(0));
  // Apply the mask to the image
  return imageCollection.updateMask(clearSkyMask);
}

// 4. Import Landsat 8 Surface Reflectance Collection
var landsat = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
  .filterBounds(aoi)
  .filterDate('2026-03-01', '2026-05-20')
  .filter(ee.Filter.lt('CLOUD_COVER', 10)) 
  .map(applyImageScaling)                  
  .map(applyCloudMask);                   

// 5. Select the least cloudy image
var image = landsat.sort('CLOUD_COVER').mean();

Map.addLayer(image.clip(aoi), {bands: ['SR_B4','SR_B3','SR_B2'], min:0, max:0.3}, "True Color");


var ndvi = image.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');
//L = 0.5 (soil brightness correction factor)
var savi = image.expression(
  '((NIR - RED) / (NIR + RED + L)) * (1 + L)', {
    'NIR': image.select('SR_B5'),
    'RED': image.select('SR_B4'),
    'L': 0.5
}).rename('SAVI');

var ndwi = image.normalizedDifference(['SR_B3', 'SR_B5']).rename('NDWI');

// 9. Add layers to the Map
Map.addLayer(ndvi.clip(aoi), {min:-1, max:1, palette:['blue','white','green']}, 'NDVI');
Map.addLayer(savi.clip(aoi), {min:-1, max:1, palette:['brown','yellow','green']}, 'SAVI');
Map.addLayer(ndwi.clip(aoi), {min:-1, max:1, palette:['white','blue']}, 'NDWI');

var ndbi = image.normalizedDifference(['SR_B6', 'SR_B5']).rename('NDBI');
var vegfrac = image.expression(
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
  
var buildfrac = image.expression(
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




//  Copernicus Global instead of Nasa since Copernicus is more recent, 2011-2015.
var dem = ee.ImageCollection('COPERNICUS/DEM/GLO30')
  .select('DEM')
  .mosaic()
  .rename('Elevation');
var slope = ee.Terrain.slope(dem).rename('Slope');
var aspect = ee.Terrain.aspect(dem).rename('Aspect');

// tpi with 100 m radius to calculate topographic index
var focalMean = dem.focal_mean(100, 'circle', 'meters');
var tpi = dem.subtract(focalMean).rename('TPI');

Map.addLayer(dem.clip(aoi), {min:0, max:300, palette:['blue','green','yellow','brown','white']}, 'Elevation', false);

// for images

var Hillshade= ee.Terrain.hillshade(dem);

var albedo = image.expression(
'((0.356 * B2) + (0.130 * B4) + (0.373 * B5) + (0.085 * B6) + (0.072 * B7) - 0.0018) / 1.016',
{
      'B2': image.select('SR_B2'), // Blue
      'B4': image.select('SR_B4'), // Red
      'B5': image.select('SR_B5'), // NIR
      'B6': image.select('SR_B6'), // SWIR 1
      'B7': image.select('SR_B7')  // SWIR 2
    }
  )
  .clamp(0, 1)        // normalizing 0-1 percentage value
  .rename('Albedo')
  .toFloat();

// LST 

//function to do additional filtering and resampling at 30m with bilinear sampling.
// LST 

// Function to filter, scale, mask, and smoothly resample 100m thermal to 30m
function prepLST(img) {
  var qa = img.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 3).eq(0).and(qa.bitwiseAnd(1 << 4).eq(0));
  var lstCelsius = img.select('ST_B10')
                        .multiply(0.00341802).add(149.0)
                        .subtract(273.15)
                        .rename('LST_30m');
                        
  // mask and bilinear resampling to go from 100m to 30m
  return lstCelsius.updateMask(mask).resample('bilinear');
}

var lst = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
  .filterBounds(aoi)
  .filterDate('2026-03-01', '2026-05-20')
  .filter(ee.Filter.lt('CLOUD_COVER', 10))
  .map(prepLST) 
  .mean();
  
Map.addLayer(lst.clip(aoi), {min: 15, max: 40, palette: ['blue', 'yellow', 'red']}, 'LST 30m', false);

//  csv export

var finalStack = ee.Image([
  ndvi, savi, ndwi, ndbi,       
  vegfrac, buildfrac, 
  distWater,                   
  dem, slope, aspect, tpi,albedo, lst     
]).toFloat();

var extractedData = finalStack.reduceRegions({
  collection: bufferedTable, // buffered areas defined at the top
  reducer: ee.Reducer.mean(), 
  scale: 30                   // Forces extraction at 30m resolution
});


Export.table.toDrive({
  collection: extractedData,
  description: 'Rome_trainingLandsat_90radius_CSV',
  folder: 'GEE_Exports_30m', 
  fileFormat: 'CSV'
});


// Images exports

//  Spectral Indices 
Export.image.toDrive({
  image: ee.Image([ndvi, savi, ndwi, ndbi]).toFloat(),
  description: 'Landsat_Spectral_Indices_90rad',
  folder: 'GEE_Exports_30m',
  region: aoi,
  scale: 30,
  maxPixels: 1e13
});

//  Topography 
Export.image.toDrive({
  image: ee.Image([dem, slope, aspect, tpi]).toFloat(),
  description: 'Topography_90rad',
  folder: 'GEE_Exports_30m',
  region: aoi,
  scale: 30,
  maxPixels: 1e13
});

// Land coverage metrics 

Export.image.toDrive({
  image: ee.Image([vegfrac, buildfrac, distWater]).toFloat(),
  description: 'LandCover_Metrics_90rad',
  folder: 'GEE_Exports_30m',
  region: aoi,
  scale: 30,
  maxPixels: 1e13
});


// LST

Export.image.toDrive({
  image: lst.toFloat(),
  description: 'LST_30m_90rad',
  folder: 'GEE_Exports_30m',
  region: aoi,
  scale: 30,
  maxPixels: 1e13
});
