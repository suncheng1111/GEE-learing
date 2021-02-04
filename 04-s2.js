
//****************************************************************************************************************
/                                                     VERSION01
//****************************************************************************************************************
//02-根据点数据得到对应的图像
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
//按照时间和位置筛选影像，加载的T1_SR影像，已经过大气校正。
//按照shape筛选
var collection = ee.ImageCollection('COPERNICUS/S2').filterBounds(table3)
.filterDate('2020-04-01','2020-05-01')
// .map(maskS2clouds)
.sort('CLOUD_COVER')
.min()
.clip(table3);

print('collection',collection);

var visParams = {bands:['B4','B3','B2'],min:0,max:2000};
Map.addLayer(collection,visParams,'432');

//****************************************************************************************************************
/                                                     VERSION02
//****************************************************************************************************************
//02-根据点数据得到对应的图像
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
//按照时间和位置筛选影像，加载的T1_SR影像，已经过大气校正。
//按照shape筛选
var s2 = ee.ImageCollection('COPERNICUS/S2').filterBounds(table3)
.filterDate('2020-04-01','2020-05-01')
// .map(maskS2clouds)
.sort('CLOUD_COVER')
.map(WI_2015s2)
.min()
.clip(table3)
;

print('s2',s2);
var visParams = {bands:['B4','B3','B2'],min:0,max:2000};
Map.addLayer(s2,visParams,'432');

    Map.addLayer(
      s2.select("WI_2015")
            .updateMask(s2.select("WI_2015").gte(0.1)), 
      {palette: "ff0000"}, 
      "WI_2015s2",
      false
    );
    
    
    
var s2 = s2.select("WI_2015");

//******************
//面积计算
//********************

//01-先完成掩膜
var water=s2.where(s2.lte(0),0).where(s2.gt(0),1)//01-区分01
var myMask=water.eq(1);
water=water.updateMask(myMask);
var viz = {min:0, max:50, palette: 'red'};
Map.addLayer(water,viz,'water',false)


//02-面积计算
//面积的计算
var areas2=water.reduceRegion({
  reducer:ee.Reducer.sum(),
  geometry:table3,
  scale:30,
  maxPixels:1E13
});
print(areas2)
