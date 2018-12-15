var points = null;
var scalarField = null;

function plot(drData,svg){
    points=drData.points;
    scalarField = drData.scalarField;
    var n = Math.sqrt(scalarField.length);

    values = scalarField;
    
    var maxVal = Math.max.apply(Math,values);
    var minVal = Math.min.apply(Math,values);
    var numLines = 11;
    var contours = d3.contours()
    .size([n, n])
    .thresholds(d3.range(1, numLines).map(p => (maxVal-minVal)*p/(numLines*1.0) + minVal))
    (values);
    
    
    var color = d3.scaleSequential(d3.interpolateGreys)
    .domain([minVal, maxVal*2]);
    
    svg.selectAll("*").remove();
    if (contOn){
        svg.selectAll("path")
        .attr("viewBox", [0, 0, size, size])
        .data(contours)
        .enter().append("path")
        .attr("d", d3.geoPath(d3.geoIdentity().scale(size / n)))
        .attr("fill", function(d) { return color(d.value); });
    }
    
    var x = d3.scaleLinear();
    var y = d3.scaleLinear();
    x.range([padding / 2, size - padding / 2]);
    y.range([size - padding / 2, padding / 2]);
    x.domain(d3.extent(points,function(d){return d.range[0];}));
    y.domain(d3.extent(points,function(d){return d.range[1];}));
    
    svg.selectAll("circle")
    .data(points)
    .enter().append("circle")
    .attr("cx",function(d){return x(d.range[0]);})
    .attr("cy",function(d){return y(d.range[1]);})
    .attr("r",5)
    .on("mouseover",function(){
        //show bar chart of perturbation
    });
    
    if(vectOn){
        svg.selectAll("path")
        .data(points)
        .enter().append("path")
        .attr("stroke","black")
        .attr("d",function(d){
            var factor = 10;
            var coords = [x(d.range[0]),y(d.range[1])];
            var normConst =Math.sqrt(Math.pow(d.outputPert[0],2)+Math.pow(d.outputPert[1],2));
            console.log(d.outputPert);
            var normV = [d.outputPert[0]/normConst,d.outputPert[1]/(normConst)];
            
            return "M"+coords[0]+" "+coords[1]+" L"+(coords[0]+factor*normV[0])+" "+(coords[1]+factor*normV[1])+" Z";
            });
    }
    
    
}

