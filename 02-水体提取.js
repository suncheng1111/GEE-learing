  
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
  var l8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR")
      .filterBounds(geometry)   
      .filterDate(startDate,endDate)
      // .map(maskL8sr)
      .map(NDWIl8)
      .map(MNDWIl8)
      .map(AWEI_nshl8)
      .map(AWEI_shl8)
      .map(WI_2015l8);
      
     l8 = l8.mosaic().clip(table);
    // l8=l8.multiply(10000); 
  print('l8',l8);
  // var l8 = l8.mosaic().clip(table);
  // var l8 = l8.median().clip(table);
  
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
  Map.addLayer(table2,{color:'blue'},'rivernet');
  
  
