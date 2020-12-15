//按照时间和位置筛选影像，加载的T1_SR影像，已经过大气校正。
//按照shape筛选
var collection = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR').filterBounds(table)
.filterDate('2020-06-01','2020-09-17');
// 按照云量排序
var sorted = collection.sort('CLOUD_COVER');
var scene = sorted.first();
print('collection',collection);
print('sorted',sorted);
print('scene',scene);
var visParams = {bands:['B4','B3','B2'],min:0,max:2000};
Map.addLayer(scene,visParams,'432');
var image1=ee.Image('LANDSAT/LC08/C01/T1_SR/LC08_119039_20200908')
Map.addLayer(image1)
Map.centerObject(image1)
