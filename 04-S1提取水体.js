Map.centerObject(geometry, 8);
//select VH band.


//load layers.
// Map.addLayer(s1, {min: -30, max: 5}, 's1', true);

var startDate = "2020-4-1";
var endDate = "2020-5-1";
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
var s1 = s1.mosaic().clip(table).select("SDWI");
Map.addLayer(s1,{color:'blue'},'s1');
// Map.addLayer(table2,{color:'blue'},'rivernet');
