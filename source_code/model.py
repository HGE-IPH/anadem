import ee

# Import collections 
col_palsar ="JAXA/ALOS/PALSAR/YEARLY/SAR"
col_sentinel = 'COPERNICUS/S1_GRD'
col_sentinel2 = 'COPERNICUS/S2_HARMONIZED'
col_gedi = "LARSE/GEDI/GEDI02_A_002_MONTHLY"
copdem = ee.ImageCollection("COPERNICUS/DEM/GLO30").select('DEM').mean()
l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')

def input_model(start_date,end_date,geometry):

    '''
    Digital Elevation Model vegetation removal

    Function to generate the input image 

    Version 0.1

    '''

    #palsar_hh_mosaic = get_palsar_mosaic(col_palsar,'HH',start_date,end_date,'mean')
    #palsar_hv_mosaic = get_palsar_mosaic(col_palsar,'HV',start_date,end_date,'mean')
  
    #ratio_palsar = radar_ratio(palsar_hv_mosaic,palsar_hh_mosaic,'ratio_palsar')
    
    slope = ee.Terrain.slope(copdem).rename('slope')

    #sen_vh_mosaic = get_sentinel_mosaic(col_sentinel,'VH',start_date,end_date,'median')
    #sen_vv_mosaic = get_sentinel_mosaic(col_sentinel,'VV',start_date,end_date,'median')
  
    s2_col = ee.ImageCollection(col_sentinel2).filterDate(start_date, end_date)\
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',50))\
    .select(['B.*','QA60'])\
    .map(maskS2clouds)\
    .map(sentinel_addNDVI)\
    .map(sentinel_addNDWI)\
    .map(sentinel_addNDMI)\
    .map(sentinel_addEVI)\
    .map(sentinel_addSR)\
    .map(sentinel_addMSR)\
    .map(sentinel_addSAVI)\
    .map(sentinel_addMSAVI)

    l8_col = l8.filterBounds(geometry)\
      .filterDate(start_date,end_date)\
      .filterMetadata('CLOUD_COVER','less_than',50)\
      .sort('CLOUD_COVER')\
      .map(prepSrL8)\
      .map(l8_addNDVI)\
      .map(l8_addNDWI)\
      .map(l8_addNDMI)\
      .map(l8_addEVI)\
      .map(l8_addSR)\
      .map(l8_addMSR)\
      .map(l8_addSAVI)\
      .map(l8_addMSAVI)
  
# Selects bands from each dataset
    l8_med = l8_col.median().select(['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7',
        'ST_B.*','EVI_L8','NDVI_L8','NDMI_L8','NDWI_L8','SAVI_L8','MSAVI_L8','SR_L8','MSR_L8' ])
 
    s2_med = s2_col.median().select(['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12',
         'EVI_S2','NDVI_S2','NDMI_S2','NDWI_S2','SAVI_S2','MSAVI_S2','SR_S2','MSR_S2'])
  
    #l8_iqr = get_iqr(l8_col.select(['EVI_L8','NDVI_L8','NDMI_L8','NDWI_L8','SAVI_L8','MSAVI_L8','SR_L8','MSR_L8']),
    #             ee.List(['EVI_L8','NDVI_L8','NDMI_L8','NDWI_L8','SAVI_L8','MSAVI_L8','SR_L8','MSR_L8']))
    #s2_iqr = get_iqr(s2_col.select(['EVI_S2','NDVI_S2','NDMI_S2','NDWI_S2','SAVI_S2','MSAVI_S2','SR_S2','MSR_S2']),
    #            ee.List(['EVI_S2','NDVI_S2','NDMI_S2','NDWI_S2','SAVI_S2','MSAVI_S2','SR_S2','MSR_S2']))
    
# Concatenates images as bands to a single image
  
    water = ee.ImageCollection('COPERNICUS/DEM/GLO30').select('WBM').mean()
    water_mask = water.gt(0)


    def convert(image):
        output = image.expression(
        '10*logDN - 83.0',{
        'logDN':image.pow(2).log10()})
        return output
    
    def radar_black_scatter(hv_band,hh_band,rename):
        return hv_band.expression(
            '(vv-hh)**2',{'vv':hv_band,'hh':hh_band}
            ).rename(rename)

    alos = ee.ImageCollection("JAXA/ALOS/PALSAR-2/Level2_2/ScanSAR")
    alos_hh = alos.select("HH").map(convert).mean().unmask(0)
    alos_hv = alos.select("HV").map(convert).mean().unmask(0)

    ratio = alos_hv.divide(alos_hh)
    radar = radar_black_scatter(alos_hv,alos_hh,'ratio').unmask(0)

    water_mask = water_mask.unmask(0).rename("water_mask")

  
    image = ee.Image.cat([
    l8_med,
    s2_med,
    #l8_iqr,
    #s2_iqr,
    #palsar_hv_mosaic,
    #palsar_hh_mosaic,
    #ratio_palsar,
    #sen_vv_mosaic,
    #sen_vh_mosaic,
    copdem.rename('copernicus'),
    slope,
    water_mask
    ])

    return image

# Ancillary functions

def get_iqr(col, bands_list):
    
    iqr= col.reduce(ee.Reducer.percentile([25, 75]))
    
    #var bands_list = col.first().bandNames()

    def band_fun(band):
        band = ee.String(band);\
        band_name = ee.String(band.split('_').get(1))
        return iqr.select(band.cat('_p75')).subtract(iqr.select(band.cat('_p25'))).rename(band.cat('_IQR'))

    
    iqr_image_list = bands_list.map(band_fun).flatten()

    def merge(image,previous):
      return ee.Image(previous).addBands(image)
    

    return ee.ImageCollection(iqr_image_list).iterate(merge,ee.Image([]))


# Landsat

def l8_addEVI(image):
  return image.addBands(ee.Image(image)
  .expression('2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)', {
        'nir': image.select('SR_B5'),
        'red': image.select('SR_B4'),
        'blue': image.select('SR_B2')}).rename('EVI_L8'))


def l8_addNDWI(image):
  return image.addBands(image.normalizedDifference(['SR_B5', 'SR_B2']).rename('NDWI_L8'))


def l8_addNDVI(image):
  ndvi = image.addBands(image.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI_L8'))
  return ndvi


def l8_addNDMI(image):
  ndvi = image.addBands(image.normalizedDifference(['SR_B5', 'SR_B6']).rename('NDMI_L8'))
  return ndvi


def l8_addSR(image):
  return image.addBands(ee.Image(image)
  .expression('nir/red', {
        'nir': image.select('SR_B5'),
        'red': image.select('SR_B4')
  })
  .rename('SR_L8'))

def l8_addMSR(image):
  return image.addBands(ee.Image(image)
  .expression('((nir/red) - 1) / ((nir/red)**0.5+1)', {
        'nir': image.select('SR_B5'),
        'red': image.select('SR_B4')
  })
  .rename('MSR_L8'))

def l8_addMSAVI(image):
  return image.addBands(ee.Image(image)
  .expression('(2 *nir + 1 - ((2 *nir + 1)**2 - 8*(nir-red))**0.5)/2', {
        'nir': image.select('SR_B5'),
        'red': image.select('SR_B4')
  })
  .rename('MSAVI_L8'))

def l8_addSAVI(image):
  return image.addBands(ee.Image(image)
  .expression(
    '((NIR - RED) / (NIR + RED + 0.5)) * 1.5', {
      'NIR': image.select('SR_B5'),
      'RED': image.select('SR_B4'),
    }).rename('SAVI_L8'))


def prepSrL8(image):
  #Develop masks for unwanted pixels (fill, cloud, cloud shadow).
  #qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0)
  saturationMask = image.select('QA_RADSAT').eq(0)
  c01 = image.select('QA_PIXEL').eq(21824); #Clear, low confidence cloud
  c02 = image.select('QA_PIXEL').eq(21888); #water, low confidence cloud
  mask = c01.Or(c02)

  #pply the scaling factors to the appropriate bands.
  def getFactorImg (factorNames):
    factorList = image.toDictionary().select(factorNames).values()
    return ee.Image.constant(factorList)
  
  scaleImg = getFactorImg([
    'REFLECTANCE_MULT_BAND_.|TEMPERATURE_MULT_BAND_ST_B10'])
  offsetImg = getFactorImg([
    'REFLECTANCE_ADD_BAND_.|TEMPERATURE_ADD_BAND_ST_B10'])
  scaled = image.select('SR_B.|ST_B10').multiply(scaleImg).add(offsetImg)

  # Replace original bands with scaled bands and apply masks.
  return image.addBands(scaled, None, True).updateMask(saturationMask).updateMask(mask)


# Sentinel

def maskS2clouds(image):
  qa = image.select('QA60')

  # Bits 10 and 11 are clouds and cirrus, respectively.
  cloudBitMask = 1 << 10
  cirrusBitMask = 1 << 11

  #Both flags should be set to zero, indicating clear conditions.
  mask = qa.bitwiseAnd(cloudBitMask).eq(0)\
      .And(qa.bitwiseAnd(cirrusBitMask).eq(0))

  return image.updateMask(mask).divide(10000)



def sentinel_addEVI(image):
  return image.addBands(ee.Image(image).expression('2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)', {
        'nir': image.select('B8'),
        'red': image.select('B4'),
        'blue': image.select('B2')}).rename('EVI_S2'))


def sentinel_addNDWI(image):
  return image.addBands(image.normalizedDifference(['B8', 'B2']).rename('NDWI_S2'))


def sentinel_addNDVI(image):
  return image.addBands(image.normalizedDifference(['B8', 'B4']).rename('NDVI_S2'))

def sentinel_addNDMI(image):
  return image.addBands(image.normalizedDifference(['B8A', 'B11']).rename('NDMI_S2'))


def sentinel_addSR(image):
  return image.addBands(ee.Image(image)
  .expression('nir/red', {
        'nir': image.select('B8'),
        'red': image.select('B4')
  }).rename('SR_S2'))


def sentinel_addMSR(image):
  return image.addBands(ee.Image(image)
  .expression('((nir/red) - 1) / ((nir/red)**0.5+1)', {
        'nir': image.select('B8'),
        'red': image.select('B4')
  }).rename('MSR_S2'))


def sentinel_addSAVI(image):
  return image.addBands(ee.Image(image)
  .expression(
    '((NIR - RED) / (NIR + RED + 0.5)) * 1.5', {
      'NIR': image.select('B8'),
      'RED': image.select('B4'),
    }).rename('SAVI_S2'))


def sentinel_addMSAVI(image):
  return image.addBands(ee.Image(image)
  .expression('(2 *nir + 1 - ((2 *nir + 1)**2 - 8*(nir-red))**0.5)/2', {
        'nir': image.select('B8'),
        'red': image.select('B4')
  })
  .rename('MSAVI_S2'))



def sentinel_indexes(start_date,end_date):
  s2_col = ee.ImageCollection(col_sentinel2).filterDate(start_date, end_date)\
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',20))\
    .select(['B.*','QA60'])\
    .map(maskS2clouds)

  
  s2_mean = s2_col.mean()
  return sentinel_addNDWI(sentinel_addEVI(sentinel_addNDVI(s2_mean)))


def get_gedi_data(collection,geometry):
  
  # ex:'LARSE/GEDI/GEDI02_A_002/GEDI02_A_2019114183224_O02064_04_T05177_02_003_01_V002'
    gedi = ee.FeatureCollection(collection)\
            .filterBounds(geometry)\
            .filter(ee.Filter.eq('quality_flag',1))\
            .filter(ee.Filter.eq('degrade_flag',0))\
            .filter(ee.Filter.eq('elevation_bias_flag',0))
    
    return gedi
 

def calc_forest_offset(image, bands, classifier_name, numberoftrees, samples, max_samples, predictor):
  
  # Samples 
  
    random_samples = samples.randomColumn('random') #.filterBounds(image.geometry());

    training = random_samples.filter(ee.Filter.lt('random', 0.7)).limit(max_samples,'random')

    def conditional(feat):
        return ee.Algorithms.If(ee.Number(feat.get(predictor)).lte(5),
        feat.set({predictor: 0}),
        feat)

    #var training = training.map(conditional);

    # Classifier
    if classifier_name == "random_forest":

        classifier = ee.Classifier.smileGradientTreeBoost(
            numberOfTrees=numberoftrees,
            variablesPerSplit=2,
            minLeafPopulation=10,
            bagFraction=0.5,
            #maxNodes=,
            #  seed=2
        )


    if classifier_name == "tree_boost":

        classifier = ee.Classifier.smileGradientTreeBoost(
        numberOfTrees=numberoftrees,
        shrinkage= 0.1, #0.1
        samplingRate= 0.6, #0.6,
        #maxNodes= 8,
        loss = 'LeastAbsoluteDeviation'
        )


    trained = classifier.train(training, predictor, bands)
    classification = image.select(bands).classify(trained.setOutputMode('REGRESSION'))

    return ee.List([classification.rename('classification'), trained])

  
def fvc_func(geometry):

    col_sentinel2 = 'COPERNICUS/S2_HARMONIZED'

    def calc_fvc(image):

        ndvi =image.normalizedDifference(['B8', 'B4']).rename("ndvi")

        #var red = ndvi.reduceRegion({
        #  reducer: ee.Reducer.minMax(),
        #  scale: 30,
        #  maxPixels: 10e12,
        #  geometry: ndvi.geometry()}
        #  )

        fvc = ndvi.clamp(0,1).expression(
        '(ndvi - ndvi_min)/(ndvi_max - ndvi)',{
        'ndvi': ndvi,
        'ndvi_min': 0.45,#red.getNumber('ndvi_min'),
        'ndvi_max':  0.9,#red.getNumber('ndvi_max'),

        })

        return fvc.clamp(0,1).rename("fvc")


    s2_col = ee.ImageCollection(col_sentinel2).filterDate('2020-01-01','2022-12-31')\
    .filterBounds(geometry)\
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',5))\
    .select(['B.*','QA60'])\
    .map(maskS2clouds)\
    .map(calc_fvc)

    s2_fvc = s2_col.select('fvc').mean()

    return s2_fvc