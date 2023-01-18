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
    .attr('id',function(d,i){return "circ_"+i;})
    .attr("cx",function(d){return x(d.range[0]);})
    .attr("cy",function(d){return y(d.range[1]);})
    .attr("r",5)
    .on("mouseover",function(){
        //show bar chart of perturbation
    });

    if(vectOn){
        svg.selectAll("path.vects")
        .data(points)
        .enter().append("path")
        .attr("class","vects")
        .attr("stroke","black")
        .attr('id',function(d,i){return "v_"+i;})
        .attr("d",function(d){
            var factor = (Math.max(x.domain()[1],y.domain()[1]) - Math.min(x.domain()[0],y.domain()[0]))/20
            var coords = [x(d.range[0]),y(d.range[1])];
            dx = x(d.outputPert[0])
            dy = y(d.outputPert[1])

            var normConst =Math.sqrt(Math.pow(d.outputPert[0],2)+Math.pow(d.outputPert[1],2));
            // var normConst = 1
            console.log(d.outputPert);

            var normV = [factor*d.outputPert[0]/normConst,factor*d.outputPert[1]/(normConst)];
            coords_end = [x(d.range[0] + normV[0]), y(d.range[1]+normV[1])]

            return "M"+coords[0]+" "+coords[1]+" L"+coords_end[0]+" "+coords_end[1]+" Z";
            });
    }


}
