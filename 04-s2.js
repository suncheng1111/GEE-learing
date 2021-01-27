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
