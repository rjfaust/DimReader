var drData = null;
var width = 800;
var padding = 20;
var size = width-padding;
var contOn = true;
var vectOn = false;

btnDiv=d3.select("#buttons");
btnDiv.append("label")
.html("Select DimReader File");

btnDiv.append("input")
.attr("type","file")
.attr("id","drFile")
.attr("accept",".dimreader")
.style("margin","5px")
.on("change",function(){
    var file = d3.event.target.files[0];
    if (file) {
    var reader = new FileReader();
    reader.onloadend = function(evt) {
    var text = evt.target.result;
    // The following call results in an "Access denied" error in IE.
    loadInputFile(text);
    };
    reader.readAsText(file);
    }
    });

btnDiv.append("label")
.html("Show Contours");

btnDiv.append("input")
.attr("type","checkbox")
.attr("id","contours")
.style("margin","5px")
.attr("checked",true)
.on("change",function(){
   if(contOn){
      contOn=false;
      plot(drData,svg);
   }
   else{
      contOn = true;
      plot(drData,svg);
   }
   
});



btnDiv.append("label")
.html("Show Vectors");

btnDiv.append("input")
.attr("type","checkbox")
.attr("id","vectors")
.style("margin","5px")
.on("change",function(){
    if(vectOn){
      vectOn=false;
      plot(drData,svg);
   }
   else{
      vectOn = true;
      plot(drData,svg);
   }
});

var svg = d3.select("#plots")
.append("svg")
.attr("width",size+padding)
.attr("height",size+padding);

function loadInputFile(text){
   drData = JSON.parse(text);
   plot(drData,svg);
   }


