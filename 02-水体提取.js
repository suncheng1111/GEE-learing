Map.centerObject(geometry);
var startDate = "2020-4-1";
var endDate = "2020-5-1";
var otsu = function(histogram) {
    var counts = ee.Array(ee.Dictionary(histogram).get('histogram'));
    var means = ee.Array(ee.Dictionary(histogram).get('bucketMeans'));
    var size = means.length().get([0]);
    var total = counts.reduce(ee.Reducer.sum(), [0]).get([0]);
    var sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]);
    var mean = sum.divide(total);
    var indices = ee.List.sequence(1, size);
    // Compute between sum of squares, where each mean partitions the data.
    var bss = indices.map(function(i) {
      var aCounts = counts.slice(0, 0, i);
      var aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0]);
      var aMeans = means.slice(0, 0, i);
      var aMean = aMeans.multiply(aCounts)
          .reduce(ee.Reducer.sum(), [0]).get([0])
          .divide(aCount);
      var bCount = total.subtract(aCount);
      var bMean = sum.subtract(aCount.multiply(aMean)).divide(bCount);
      return aCount.multiply(aMean.subtract(mean).pow(2)).add(
            bCount.multiply(bMean.subtract(mean).pow(2)));
    });
    return means.sort(bss).get([-1]);
  };
  
// // /************************indexa and function********************* */
//NDWI: (G - NIR)/(G + NIR)
  function NDWIl8(image) {
    return image.addBands(
      image.normalizedDifference(["B3", "B5"])
          .rename("NDWI")
    );
  }
  //MNDWI: (G - SWIR)/(G + SWIR)
  function MNDWIl8(image) {
    return image.addBands(
      image.normalizedDifference(["B3", "B7"])
          .rename("MNDWI")
    );
  }
  //AWEI_nsh: 4 * (G - SWIR1) - (0.25 * NIR + 2.75 * SWIR2)
  function AWEI_nshl8(image) {
    var awei_nsh = image.expression(
      "4 * (G - SWIR1) - (0.25 * NIR + 2.75 * SWIR2)",
      {
        "G": image.select("B3"),
        "NIR": image.select("B5"),
        "SWIR1": image.select("B6"),
        "SWIR2": image.select("B7")
      }
    );
    return image.addBands(awei_nsh.rename("AWEI_nsh"));
  }
  //AWEI_sh: B + 2.5 * G - 1.5 * (NIR + SWIR1) - 0.25 * SWIR2
  function AWEI_shl8(image) {
    var awei_sh = image.expression(
      "B + 2.5 * G - 1.5 * (NIR + SWIR1) - 0.25 * SWIR2",
      {
        "B": image.select("B2"),
        "G": image.select("B3"),
        "NIR": image.select("B5"),
        "SWIR1": image.select("B6"),
        "SWIR2": image.select("B7")
      }
    );
    return image.addBands(awei_sh.rename("AWEI_sh"));
  }    
  //WI_2015: 1.7204 + 171*G + 3*R - 70*NIR - 45*SWIR1 - 71*SWIR2
  function WI_2015l8(image) {
    var wi_2015 = image.expression(
      "1.7204 + 171*G + 3*R - 70*NIR - 45*SWIR1 - 71*SWIR2",
      {
        "G": image.select("B3"),
        "R": image.select("B4"),
        "NIR": image.select("B5"),
        "SWIR1": image.select("B6"),
        "SWIR2": image.select("B7")
      }
    );
    return image.addBands(wi_2015.rename("WI_2015"));
  }

  
// // /************************datacollection********************* */
  
function maskL8sr(image) {
  var cloudShadowBitMask = 1 << 3;
  var cloudsBitMask = 1 << 5;
  var qa = image.select('pixel_qa');
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
      .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return image.updateMask(mask)
      .select("B[0-9]*")
      .copyProperties(image, ["system:time_start"]);
}
var l8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
    .filterBounds(geometry)   
    .filterDate(startDate,endDate)
    .map(NDWIl8)
    .map(MNDWIl8)
    .map(AWEI_nshl8)
    .map(AWEI_shl8)
    .map(WI_2015l8);

print('l8',l8);
// var l8= maskL8sr(l8);
var l8 = l8.mosaic().clip(table)

var visParamsl8 = {bands:['B4','B3','B2'],min:0,max:2000};

Map.addLayer(l8,visParamsl8,'l8');
  //===========================阈值提取 ==================

  print(ui.Chart.image.histogram({
    image: l8.select("NDWI"),
    region: table,
    scale:500
  }));
  print(ui.Chart.image.histogram({
    image: l8.select("MNDWI"),
    region: table,
    scale:500
  }));
  print(ui.Chart.image.histogram({
    image: l8.select("AWEI_nsh"),
    region: table,
    scale:500
  }));
  print(ui.Chart.image.histogram({
    image: l8.select("AWEI_sh"),
    region: table,
    scale:500
  }));
  print(ui.Chart.image.histogram({
    image: l8.select("WI_2015"),
    region: table,
    scale:500
  }));

  //add2*
  var histogramNDWIl8 = l8.select("NDWI")
                        .reduceRegion({
                          reducer: ee.Reducer.histogram(), 
                          geometry: table, 
                          scale: 30,
                          maxPixels: 1e13
                        });
  var histogramMNDWIl8 = l8.select("MNDWI")
                        .reduceRegion({
                          reducer: ee.Reducer.histogram(), 
                          geometry: table, 
                          scale: 30,
                          maxPixels: 1e13
                        });
  var histogramAWEI_nshl8 = l8.select("AWEI_nsh")
                        .reduceRegion({
                          reducer: ee.Reducer.histogram(), 
                          geometry:table , 
                          scale: 30,
                          maxPixels: 1e13
                        }); 
  var histogramAWEI_shl8 = l8.select("AWEI_sh")
                        .reduceRegion({
                          reducer: ee.Reducer.histogram(), 
                          geometry: table, 
                          scale: 30,
                          maxPixels: 1e13
                        });           
  var histogramWI_2015l8 = l8.select("WI_2015")
                        .reduceRegion({
                          reducer: ee.Reducer.histogram(), 
                          geometry: table, 
                          scale: 30,
                          maxPixels: 1e13
                        });                        
  //add3*                      
  var thresholdNDWIl8 = otsu(histogramNDWIl8.get("NDWI"));
  var thresholdMNDWIl8 = otsu(histogramMNDWIl8.get("MNDWI"));
  var thresholdAWEI_nshl8= otsu(histogramAWEI_nshl8.get("AWEI_nsh"));
  var thresholdAWEI_shl8 = otsu(histogramAWEI_shl8.get("AWEI_sh"));
  var thresholdWI_2015l8 = otsu(histogramWI_2015l8.get("WI_2015"));
  //add4*
  print("thresholdNDWIl8", thresholdNDWIl8);
  print("thresholdMNDWIl8", thresholdMNDWIl8);
  print("thresholdAWEI_nshl8", thresholdAWEI_nshl8);
  print("thresholdAWEI_shl8", thresholdAWEI_shl8);
  print("thresholdWI_2015l8", thresholdWI_2015l8);
  //add5*
  Map.addLayer(
    l8.select("NDWI")
          .updateMask(l8.select("NDWI").gte(thresholdNDWIl8)), 
    {palette: "ff0000"}, 
    "NDWIl8",
    false
  );
  Map.addLayer(
    l8.select("MNDWI")
          .updateMask(l8.select("MNDWI").gte(thresholdNDWIl8)), 
    {palette: "ff0000"}, 
    "MNDWIl8",
    false
  );
  Map.addLayer(
    l8.select("AWEI_nsh")
          .updateMask(l8.select("AWEI_nsh").gte(thresholdAWEI_nshl8)), 
    {palette: "ff0000"}, 
    "AWEI_nsh",
    false
  );
  Map.addLayer(
    l8.select("AWEI_sh")
          .updateMask(l8.select("AWEI_sh").gte(thresholdAWEI_shl8)), 
    {palette: "ff0000"}, 
    "AWEI_shl8",
    false
  );
  Map.addLayer(
    l8.select("WI_2015")
          .updateMask(l8.select("WI_2015").gte(thresholdWI_2015l8)), 
    {palette: "ff0000"}, 
    "WI_2015l8",
    false
  );
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

  function NDWIs2(image) {
    return image.addBands(
      image.normalizedDifference(["B3", "B8"])
          .rename("NDWI")
    );
  }
  //MNDWI: (G - SWIR)/(G + SWIR)
  function MNDWIs2(image) {
    return image.addBands(
      image.normalizedDifference(["B3", "B12"])
          .rename("MNDWI")
    );
  }
  //AWEI_nsh: 4 * (G - SWIR1) - (0.25 * NIR + 2.75 * SWIR2)
  function AWEI_nshs2(image) {
    var awei_nsh = image.expression(
      "4 * (G - SWIR1) - (0.25 * NIR + 2.75 * SWIR2)",
      {
        "G": image.select("B3"),
        "NIR": image.select("B8"),
        "SWIR1": image.select("B11"),
        "SWIR2": image.select("B12")
      }
    );
    return image.addBands(awei_nsh.rename("AWEI_nsh"));
  }
  //AWEI_sh: B + 2.5 * G - 1.5 * (NIR + SWIR1) - 0.25 * SWIR2
  function AWEI_shs2(image) {
    var awei_sh = image.expression(
      "B + 2.5 * G - 1.5 * (NIR + SWIR1) - 0.25 * SWIR2",
      {
        "B": image.select("B2"),
        "G": image.select("B3"),
        "NIR": image.select("B8"),
        "SWIR1": image.select("B11"),
        "SWIR2": image.select("B12")
      }
    );
    return image.addBands(awei_sh.rename("AWEI_sh"));
  }    
  //WI_2015: 1.7204 + 171*G + 3*R - 70*NIR - 45*SWIR1 - 71*SWIR2
  function WI_2015s2(image) {
    var wi_2015 = image.expression(
      "1.7204 + 171*G + 3*R - 70*NIR - 45*SWIR1 - 71*SWIR2",
      {
        "G": image.select("B3"),
        "R": image.select("B4"),
        "NIR": image.select("B8"),
        "SWIR1": image.select("B11"),
        "SWIR2": image.select("B12")
      }
    );
    return image.addBands(wi_2015.rename("WI_2015"));
  }
  
function maskS2clouds(image) {
  var qa = image.select('QA60');

  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

  return image.updateMask(mask);
}
var s2 = ee.ImageCollection('COPERNICUS/S2_SR')
                .filterBounds(geometry)   
                .filterDate(startDate,endDate)
                  //Pre-filter to get less cloudy granules.
                .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',20))
                .map(maskS2clouds)
                .map(NDWIs2)
                .map(MNDWIs2)
                .map(AWEI_nshs2)
                .map(AWEI_shs2)
                .map(WI_2015s2);;
var s2 = s2.median().clip(table);
print('s2',s2);
var visParams2 = {bands:['B4','B3','B2'],min:0,max:2000};
Map.addLayer(s2, visParams2, 's2');


    print(ui.Chart.image.histogram({
      image: s2.select("NDWI"),
      region: table,
      scale:500
    }));
    print(ui.Chart.image.histogram({
      image: s2.select("MNDWI"),
      region: table,
      scale:500
    }));
    print(ui.Chart.image.histogram({
      image: s2.select("AWEI_nsh"),
      region: table,
      scale:500
    }));
    print(ui.Chart.image.histogram({
      image: s2.select("AWEI_sh"),
      region: table,
      scale:500
    }));
    print(ui.Chart.image.histogram({
      image: s2.select("WI_2015"),
      region: table,
      scale:500
    }));
    var histogramNDWIs2 = s2.select("NDWI")
                          .reduceRegion({
                            reducer: ee.Reducer.histogram(), 
                            geometry: table, 
                            scale: 30,
                            maxPixels: 1e13
                          });
    var histogramMNDWIs2 = s2.select("MNDWI")
                          .reduceRegion({
                            reducer: ee.Reducer.histogram(), 
                            geometry: table, 
                            scale: 30,
                            maxPixels: 1e13
                          });
    var histogramAWEI_nshs2 = s2.select("AWEI_nsh")
                          .reduceRegion({
                            reducer: ee.Reducer.histogram(), 
                            geometry:table , 
                            scale: 30,
                            maxPixels: 1e13
                          }); 
    var histogramAWEI_shs2 = s2.select("AWEI_sh")
                          .reduceRegion({
                            reducer: ee.Reducer.histogram(), 
                            geometry: table, 
                            scale: 30,
                            maxPixels: 1e13
                          });           
    var histogramWI_2015s2 = s2.select("WI_2015")
                          .reduceRegion({
                            reducer: ee.Reducer.histogram(), 
                            geometry: table, 
                            scale: 30,
                            maxPixels: 1e13
                          });                        
    //add3*                      
    var thresholdNDWIs2 = otsu(histogramNDWIs2.get("NDWI"));
    var thresholdMNDWIs2 = otsu(histogramMNDWIs2.get("MNDWI"));
    var thresholdAWEI_nshs2= otsu(histogramAWEI_nshs2.get("AWEI_nsh"));
    var thresholdAWEI_shs2 = otsu(histogramAWEI_shs2.get("AWEI_sh"));
    var thresholdWI_2015s2 = otsu(histogramWI_2015s2.get("WI_2015"));
    //add4*
    print("thresholdNDWIs2", thresholdNDWIs2);
    print("thresholdMNDWIs2", thresholdMNDWIs2);
    print("thresholdAWEI_nshs2", thresholdAWEI_nshs2);
    print("thresholdAWEI_shs2", thresholdAWEI_shs2);
    print("thresholdWI_2015s2", thresholdWI_2015s2);
    //add5*
    Map.addLayer(
      s2.select("NDWI")
            .updateMask(s2.select("NDWI").gte(thresholdNDWIs2)), 
      {palette: "ff0000"}, 
      "NDWIs2",
      false
    );
    Map.addLayer(
      s2.select("MNDWI")
            .updateMask(s2.select("MNDWI").gte(thresholdNDWIs2)), 
      {palette: "ff0000"}, 
      "MNDWIs2",
      false
    );
    Map.addLayer(
      s2.select("AWEI_nsh")
            .updateMask(s2.select("AWEI_nsh").gte(thresholdAWEI_nshs2)), 
      {palette: "ff0000"}, 
      "AWEI_nshs2",
      false
    );
    Map.addLayer(
      s2.select("AWEI_sh")
            .updateMask(s2.select("AWEI_sh").gte(thresholdAWEI_shs2)), 
      {palette: "ff0000"}, 
      "AWEI_shs2",
      false
    );
    Map.addLayer(
      s2.select("WI_2015")
            .updateMask(s2.select("WI_2015").gte(thresholdWI_2015s2)), 
      {palette: "ff0000"}, 
      "WI_2015s2",
      false
    );





















//modis
  //NDWI: (G - NIR)/(G + NIR)
    function NDWImod(image) {
      return image.addBands(
        image.normalizedDifference(["sur_refl_b04", "sur_refl_b02"])
            .rename("NDWI")
      );
    }
    //MNDWI: (G - SWIR)/(G + SWIR)
    function MNDWImod(image) {
      return image.addBands(
        image.normalizedDifference(["sur_refl_b04", "sur_refl_b05"])
        //sur_refl_b05
// 1230-1250nm
// sur_refl_b06
// 1628-1652nm
            .rename("MNDWI")
      );
    }
    //AWEI_nsh: 4 * (G - SWIR1) - (0.25 * NIR + 2.75 * SWIR2)
    function AWEI_nshmod(image) {
      var awei_nsh = image.expression(
        "4 * (G - SWIR1) - (0.25 * NIR + 2.75 * SWIR2)",
        {
          "G": image.select("sur_refl_b04"),
          "NIR": image.select("sur_refl_b02"),
          "SWIR1": image.select("sur_refl_b05"),
          "SWIR2": image.select("sur_refl_b07")
        }
      );
      return image.addBands(awei_nsh.rename("AWEI_nsh"));
    }
    //AWEI_sh: B + 2.5 * G - 1.5 * (NIR + SWIR1) - 0.25 * SWIR2
    function AWEI_shmod(image) {
      var awei_sh = image.expression(
        "B + 2.5 * G - 1.5 * (NIR + SWIR1) - 0.25 * SWIR2",
        {
          "B": image.select("sur_refl_b03"),
          "G": image.select("sur_refl_b04"),
          "NIR": image.select("sur_refl_b02"),
          "SWIR1": image.select("sur_refl_b05"),
          "SWIR2": image.select("sur_refl_b07")
        }
      );
      return image.addBands(awei_sh.rename("AWEI_sh"));
    }    
    //WI_2015: 1.7204 + 171*G + 3*R - 70*NIR - 45*SWIR1 - 71*SWIR2
    function WI_2015mod(image) {
      var wi_2015 = image.expression(
        "1.7204 + 171*G + 3*R - 70*NIR - 45*SWIR1 - 71*SWIR2",
        {
          "G": image.select("sur_refl_b04"),
          "R": image.select("sur_refl_b01"),
          "NIR": image.select("sur_refl_b02"),
          "SWIR1": image.select("sur_refl_b05"),
          "SWIR2": image.select("sur_refl_b07")
        }
      );
      return image.addBands(wi_2015.rename("WI_2015"));
    }
var mod = ee.ImageCollection('MODIS/006/MOD09A1')
                .filterBounds(geometry)   
                .filterDate(startDate,endDate)
                      .map(NDWImod)
      .map(MNDWImod)
      .map(AWEI_nshmod)
      .map(AWEI_shmod)
      .map(WI_2015mod);
                
var mod=mod.median().clip(table);

print('mod',mod);
var visParamsmod = {bands:['sur_refl_b01','sur_refl_b04','sur_refl_b03'],min:0,max:2000};
Map.addLayer(mod, visParamsmod, 'mod');





    print(ui.Chart.image.histogram({
      image: mod.select("NDWI"),
      region: table,
      scale:500
    }));
    print(ui.Chart.image.histogram({
      image: mod.select("MNDWI"),
      region: table,
      scale:500
    }));
    print(ui.Chart.image.histogram({
      image: mod.select("AWEI_nsh"),
      region: table,
      scale:500
    }));
    print(ui.Chart.image.histogram({
      image: mod.select("AWEI_sh"),
      region: table,
      scale:500
    }));
    print(ui.Chart.image.histogram({
      image: mod.select("WI_2015"),
      region: table,
      scale:500
    }));
  
    //add2*
    var histogramNDWImod = mod.select("NDWI")
                          .reduceRegion({
                            reducer: ee.Reducer.histogram(), 
                            geometry: table, 
                            scale: 30,
                            maxPixels: 1e13
                          });
    var histogramMNDWImod= mod.select("MNDWI")
                          .reduceRegion({
                            reducer: ee.Reducer.histogram(), 
                            geometry: table, 
                            scale: 30,
                            maxPixels: 1e13
                          });
    var histogramAWEI_nshmod = mod.select("AWEI_nsh")
                          .reduceRegion({
                            reducer: ee.Reducer.histogram(), 
                            geometry:table , 
                            scale: 30,
                            maxPixels: 1e13
                          }); 
    var histogramAWEI_shmod = mod.select("AWEI_sh")
                          .reduceRegion({
                            reducer: ee.Reducer.histogram(), 
                            geometry: table, 
                            scale: 30,
                            maxPixels: 1e13
                          });           
    var histogramWI_2015mod= mod.select("WI_2015")
                          .reduceRegion({
                            reducer: ee.Reducer.histogram(), 
                            geometry: table, 
                            scale: 30,
                            maxPixels: 1e13
                          });                        
    //add3*                      
    var thresholdNDWImod = otsu(histogramNDWImod.get("NDWI"));
    var thresholdMNDWImod = otsu(histogramMNDWImod.get("MNDWI"));
    var thresholdAWEI_nshmod= otsu(histogramAWEI_nshmod.get("AWEI_nsh"));
    var thresholdAWEI_shmod = otsu(histogramAWEI_shmod.get("AWEI_sh"));
    var thresholdWI_2015mod = otsu(histogramWI_2015mod.get("WI_2015"));
    //add4*
    print("thresholdNDWImod", thresholdNDWImod);
    print("thresholdMNDWImod", thresholdMNDWImod);
    print("thresholdAWEI_nshmod", thresholdAWEI_nshmod);
    print("thresholdAWEI_shmod", thresholdAWEI_shmod);
    print("thresholdWI_2015mod", thresholdWI_2015mod);
    //add5*
    Map.addLayer(
      mod.select("NDWI")
            .updateMask(mod.select("NDWI").gte(thresholdNDWImod)), 
      {palette: "ff0000"}, 
      "NDWImod",
      false
    );
    Map.addLayer(
      mod.select("MNDWI")
            .updateMask(mod.select("MNDWI").gte(thresholdNDWImod)), 
      {palette: "ff0000"}, 
      "MNDWImod",
      false
    );
    Map.addLayer(
      mod.select("AWEI_nsh")
            .updateMask(mod.select("AWEI_nsh").gte(thresholdAWEI_nshmod)), 
      {palette: "ff0000"}, 
      "AWEI_nshmod",
      false
    );
    Map.addLayer(
      mod.select("AWEI_sh")
            .updateMask(mod.select("AWEI_sh").gte(thresholdAWEI_shmod)), 
      {palette: "ff0000"}, 
      "AWEI_shmod",
      false
    );
    Map.addLayer(
      mod.select("WI_2015")
            .updateMask(mod.select("WI_2015").gte(thresholdWI_2015mod)), 
      {palette: "ff0000"}, 
      "WI_2015mod",
      false
    );





















//landsat457===
//NDWI: (G - NIR)/(G + NIR)
    function NDWIl7(image) {
      return image.addBands(
        image.normalizedDifference(["B2", "B4"])
            .rename("NDWI")
      );
    }
    //MNDWI: (G - SWIR)/(G + SWIR)
    function MNDWIl7(image) {
      return image.addBands(
        image.normalizedDifference(["B2", "B5"])
            .rename("MNDWI")
      );
    }
    //AWEI_nsh: 4 * (G - SWIR1) - (0.25 * NIR + 2.75 * SWIR2)
    function AWEI_nshl7(image) {
      var awei_nsh = image.expression(
        "4 * (G - SWIR1) - (0.25 * NIR + 2.75 * SWIR2)",
        {
          "G": image.select("B2"),
          "NIR": image.select("B4"),
          "SWIR1": image.select("B5"),
          "SWIR2": image.select("B7")
        }
      );
      return image.addBands(awei_nsh.rename("AWEI_nsh"));
    }
    //AWEI_sh: B + 2.5 * G - 1.5 * (NIR + SWIR1) - 0.25 * SWIR2
    function AWEI_shl7(image) {
      var awei_sh = image.expression(
        "B + 2.5 * G - 1.5 * (NIR + SWIR1) - 0.25 * SWIR2",
        {
          "B": image.select("B1"),
          "G": image.select("B2"),
          "NIR": image.select("B4"),
          "SWIR1": image.select("B5"),
          "SWIR2": image.select("B7")
        }
      );
      return image.addBands(awei_sh.rename("AWEI_sh"));
    }    
    //WI_2015: 1.7204 + 171*G + 3*R - 70*NIR - 45*SWIR1 - 71*SWIR2
    function WI_2015l7(image) {
      var wi_2015 = image.expression(
        "1.7204 + 171*G + 3*R - 70*NIR - 45*SWIR1 - 71*SWIR2",
        {
          "G": image.select("B2"),
          "R": image.select("B3"),
          "NIR": image.select("B4"),
          "SWIR1": image.select("B5"),
          "SWIR2": image.select("B7")
        }
      );
      return image.addBands(wi_2015.rename("WI_2015"));
    }
    
var cloudMaskL457 = function(image) {
  var qa = image.select('pixel_qa');
  var cloud = qa.bitwiseAnd(1 << 5)
                  .and(qa.bitwiseAnd(1 << 7))
                  .or(qa.bitwiseAnd(1 << 3));
  var mask2 = image.mask().reduce(ee.Reducer.min());
  return image.updateMask(cloud.not()).updateMask(mask2);
};

var l7 = ee.ImageCollection("LANDSAT/LE07/C01/T1_SR")
    .filterBounds(geometry)
    .map(cloudMaskL457)
      .map(NDWIl7)
      .map(MNDWIl7)
      .map(AWEI_nshl7)
      .map(AWEI_shl7)
      .map(WI_2015l7)
    .filterDate(startDate,endDate);
var l7 = l7.median().clip(table);
print('l7',l7);
var visParamsl457 = {
  bands: ['B3', 'B2', 'B1'],
  min: 0,
  max: 3000,
  gamma: 1.4,
};
Map.addLayer(l7,visParamsl457,'l7');
    print(ui.Chart.image.histogram({
      image: l7.select("NDWI"),
      region: table,
      scale:500
    }));

    print(ui.Chart.image.histogram({
      image: l7.select("MNDWI"),
      region: table,
      scale:500
    }));
    print(ui.Chart.image.histogram({
      image: l7.select("AWEI_nsh"),
      region: table,
      scale:500
    }));
    print(ui.Chart.image.histogram({
      image: l7.select("AWEI_sh"),
      region: table,
      scale:500
    }));
    print(ui.Chart.image.histogram({
      image: l7.select("WI_2015"),
      region: table,
      scale:500
    }));
  
    // add2*
    var histogramNDWIl7 = l7.select("NDWI")
                          .reduceRegion({
                            reducer: ee.Reducer.histogram(), 
                            geometry: table, 
                            scale: 30,
                            maxPixels: 1e13
                          });
    var histogramMNDWIl7 = l7.select("MNDWI")
                          .reduceRegion({
                            reducer: ee.Reducer.histogram(), 
                            geometry: table, 
                            scale: 30,
                            maxPixels: 1e13
                          });
    var histogramAWEI_nshl7 = l7.select("AWEI_nsh")
                          .reduceRegion({
                            reducer: ee.Reducer.histogram(), 
                            geometry:table , 
                            scale: 30,
                            maxPixels: 1e13
                          }); 
    var histogramAWEI_shl7 = l7.select("AWEI_sh")
                          .reduceRegion({
                            reducer: ee.Reducer.histogram(), 
                            geometry: table, 
                            scale: 30,
                            maxPixels: 1e13
                          });           
    var histogramWI_2015l7 = l7.select("WI_2015")
                          .reduceRegion({
                            reducer: ee.Reducer.histogram(), 
                            geometry: table, 
                            scale: 30,
                            maxPixels: 1e13
                          });                        
    //add3*                      
    var thresholdNDWIl7= otsu(histogramNDWIl7.get("NDWI"));
    var thresholdMNDWIl7 = otsu(histogramMNDWIl7.get("MNDWI"));
    var thresholdAWEI_nshl7= otsu(histogramAWEI_nshl7.get("AWEI_nsh"));
    var thresholdAWEI_shl7 = otsu(histogramAWEI_shl7.get("AWEI_sh"));
    var thresholdWI_2015l7 = otsu(histogramWI_2015l7.get("WI_2015"));
    //add4*
    print("thresholdNDWIl7", thresholdNDWIl7);
    print("thresholdMNDWIl7", thresholdMNDWIl7);
    print("thresholdAWEI_nshl7", thresholdAWEI_nshl7);
    print("thresholdAWEI_shl7", thresholdAWEI_shl7);
    print("thresholdWI_2015l7", thresholdWI_2015l7);

    Map.addLayer(
      l7.select("NDWI")
            .updateMask(l7.select("NDWI").gte(thresholdNDWIl7)), 
      {palette: "ff0000"}, 
      "NDWIl7",
      false
    );
    Map.addLayer(
      l7.select("MNDWI")
            .updateMask(l7.select("MNDWI").gte(thresholdNDWIl7)), 
      {palette: "ff0000"}, 
      "MNDWIl7",
      false
    );
    Map.addLayer(
      l7.select("AWEI_nsh")
            .updateMask(l7.select("AWEI_nsh").gte(thresholdAWEI_nshl7)), 
      {palette: "ff0000"}, 
      "AWEI_nshl7",
      false
    );
    Map.addLayer(
      l7.select("AWEI_sh")
            .updateMask(l7.select("AWEI_sh").gte(thresholdAWEI_shl7)), 
      {palette: "ff0000"}, 
      "AWEI_shl7",
      false
    );
    Map.addLayer(
      l7.select("WI_2015")
            .updateMask(l7.select("WI_2015").gte(thresholdWI_2015l7)), 
      {palette: "ff0000"}, 
      "WI_2015l7",
      false
    );

    //first load geometry.
//using geometry center cordination and display.
Map.centerObject(geometry, 8);
//select VH band.


//load layers.
// Map.addLayer(s1, {min: -30, max: 5}, 's1', true);


  function SDWI1(image) {
    var SDWI1 = image.expression(
      "abs(VV*VH)",
      {
        "VV": image.select("VV"),
        "VH": image.select("VH"),
      }
      );
      return image.addBands(SDWI1.rename("SDWI1"));
    }
    
  function SDWI(image) {
  var SDWI = image.expression(
    "log(10*SDWI1)-8",
    {
      "SDWI1": image.select("SDWI1"),
    }
    );
    return image.addBands(SDWI.rename("SDWI"));
  }
  var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
        .filter(ee.Filter.eq('instrumentMode', 'IW'))
        .filterDate(startDate,endDate)
        .filterBounds(geometry)
        Map.addLayer(s1,{color:'blue'},'s1or');
var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
        .filter(ee.Filter.eq('instrumentMode', 'IW'))
        .filterDate(startDate,endDate)
        .filterBounds(geometry)
        .map(SDWI1)
       .map(SDWI);

print(s1,'s1')
//   print(ui.Chart.image.histogram({
//       image: s1.select("SDWI"),
//       region: table,
//       scale:500
//     }));
  
//     // add2*

//     //add4*
//     print("thresholdSDWIs1", thresholdSDWIs1);         
//mosaic.
var s1 = s1.mosaic().clip(table).select("SDWI");
Map.addLayer(s1,{color:'blue'},'s1');
Map.addLayer(table2,{color:'blue'},'rivernet');
