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

var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
        .filter(ee.Filter.eq('instrumentMode', 'IW'))
        .filterDate(startDate,endDate)
        .filterBounds(geometry)
        .map(SDWI1)
       .map(SDWI);


// print(s1,'s1')
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
// Map.addLayer(s1,{color:'blue'},'s1');


  // print(ui.Chart.image.histogram({
  //   image: s1.select("SDWI"),
  //   region: table,
  //   scale:500
  // }));
  
  
var histogramSDWIs1 = s1.select("SDWI")
                        .reduceRegion({
                          reducer: ee.Reducer.histogram(), 
                          geometry: table, 
                          scale: 30,
                          maxPixels: 1e13
                        });
 // var thresholdSDWIs1 = otsu(histogramSDWIs1.get("SDWI"));
    // print("thresholdSDWIs1", thresholdSDWIs1);

  Map.addLayer(
    s1.select("SDWI")
          .updateMask(s1.select("SDWI").gte(0)), 
    {palette: "ff0000"}, 
    "s1SDWI",
    false
  );
  
  //尝试将SDWI小于0.2的部分mask，首先区分为1和0然后把等于0的部分使用
  //updataMask函数掩模去掉 solution:updateMask会将影像上为0的区域掩膜掉
  //或者直接采用函数
  //water=water.updateMask(s2.mask().not)
  //首先重分类，在用重分的结果去掩
  
  
var water=s1.where(s1.lte(0.2),0).where(s1.gt(0.2),1)
var myMask = water.eq(0).not();
water = water.updateMask(myMask);
 var viz = {min:0, max:50, palette: 'red'};
  Map.addLayer(water,viz,false)


//面积的计算
var areas1=water.reduceRegion({
  reducer:ee.Reducer.sum(),
  geometry:table3,
  scale:30,
  maxPixels:1E13
});
print(areas1)
