
//作者:快多宝，微信kitmyfaceplease， email：kitmyfaceplease2@gmail.com

function get_yearly_water(year) {
  
    //设置需要提取的区域，由于是上传的shp文件，需要转为geometry的格式
    var yantze_down_region = yantze_down.geometry();

    //设置需要提取的年份
    var startDate = ee.Date.fromYMD(year, 4, 1);
    var endDate = ee.Date.fromYMD(year, 5, 31);
    
    //筛选JRC水体数据
    var myjrc = jrc.filterBounds(yantze_down_region).filterDate(startDate, endDate);
    
    //在每个月份的影像中添加一个obs属性的波段，一个像素如果有数据，则为1，没有数据则为0
    myjrc = myjrc.map(function(img){
      var obs = img.gt(0);
      return img.addBands(obs.rename('obs').set('system:time_start', img.get('system:time_start')));
    });

    //在每个月份的影像中添加一个onlywater属性的波段，一个像素如果有水则为1，没有水则为0
    myjrc = myjrc.map(function(img){
      var water = img.select('water').eq(2);
      return img.addBands(water.rename('onlywater').set('system:time_start', img.get('system:time_start')));
    });
    
    //计算每个像素点在一年12景影像中， 有数据的次数
    var totalObs = ee.ImageCollection(myjrc.select('obs')).sum().toFloat();
    
    //计算每个像素点在一年12景影像中， 有水的次数
    var totalWater = ee.ImageCollection(myjrc.select('onlywater')).sum().toFloat();
    
    //统计每个像素点在一年中有水的比例
    var floodfreq = totalWater.divide(totalObs).multiply(100);
    
    //删除没有值的像素
    var myMask = floodfreq.eq(0).not();
    floodfreq = floodfreq.updateMask(myMask);
    
    var viz = {min:0, max:50, palette: 'red'};
    var floodfreq1=floodfreq.clip(yantze_down_region);
    var year_folder=year+"folder_gte";
    
    //如果某个像素一年有7个月有水，则为水体
    var gte60=floodfreq1.gte(20)
    //加载范围
  //  Map.addLayer(yantze_down_region,false)
    //加载影像
    Map.addLayer(floodfreq.clip(yantze_down_region),viz,year_folder,false)
    
    //导出影像
    Export.image.toDrive({
      image: gte60,
      region: yantze_down_region,
      // fileDimensions:2560,
      scale: 30,
      maxPixels : 1e13,
      folder:year_folder,
      description:year_folder});

     
    //计算计算长江下游水体面积
    var stats2 = gte60.reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: yantze_down_region,
      scale: 30,
      maxPixels: 1E13
    });

    print(year_folder);
    print(stats2);
    }

//获取哪一年的，如果你想获取2000年到2019年，将条件改为i<20
// for(var i=1;i<5;i++){
//   if (i<10){  var year='201'+i;}
//   if (i>10||i==10){  var year='20'+i;}
//   //parseInt(12,16)=16+2=18
//   //parse() 方法可解析一个日期时间字符串，并返回 1970/1/1 午夜距离该日期时间的毫秒数。
//   var yearn = parseInt(JSON.parse(year));
//   get_yearly_water(yearn);
// }
for(var i=2014;i<2020;i++){
 var year=i;
  //parseInt(12,16)=16+2=18
  //parse() 方法可解析一个日期时间字符串，并返回 1970/1/1 午夜距离该日期时间的毫秒数。
  var yearn = parseInt(JSON.parse(year));
  get_yearly_water(yearn);
}



    Map.addLayer(yantze_down,false)






