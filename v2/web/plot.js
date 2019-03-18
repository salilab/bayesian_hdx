function addViolin(svg, results, height, width, domain, imposeMax, violinColor){


        var data = d3.layout.histogram()
                        .bins(resolution)
                        .frequency(0)
                        (results);

        var y = d3.scale.linear()
                    .range([width/2, 0])
                    .domain([0, Math.max(imposeMax, d3.max(data, function(d) { return d.y; }))]);

        var x = d3.scale.linear()
                    .range([height, 0])
                    .domain(domain)
                    .nice();


        var area = d3.svg.area()
            .interpolate(interpolation)
            .x(function(d) {
                   if(interpolation=="step-before")
                        return x(d.x+d.dx/2)
                   return x(d.x);
                })
            .y0(width/2)
            .y1(function(d) { return y(d.y); });

        var line=d3.svg.line()
            .interpolate(interpolation)
            .x(function(d) {
                   if(interpolation=="step-before")
                        return x(d.x+d.dx/2)
                   return x(d.x);
                })
            .y(function(d) { return y(d.y); });

        var gPlus=svg.append("g")
        var gMinus=svg.append("g")

        gPlus.append("path")
          .datum(data)
          .attr("class", "area")
          .attr("d", area)
          .style("fill", violinColor);

        gPlus.append("path")
          .datum(data)
          .attr("class", "violin")
          .attr("d", line)
          .style("stroke", violinColor);


        gMinus.append("path")
          .datum(data)
          .attr("class", "area")
          .attr("d", area)
          .style("fill", violinColor);

        gMinus.append("path")
          .datum(data)
          .attr("class", "violin")
          .attr("d", line)
          .style("stroke", violinColor);

        var x=width;

        gPlus.attr("transform", "rotate(90,0,0)  translate(0,-"+width+")");//translate(0,-200)");
        gMinus.attr("transform", "rotate(90,0,0) scale(1,-1)");


}

var margin={top:10, bottom:30, left:30, right:10};

var width=600;
var height=200;
var boxWidth=100;
var boxSpacing=10;
var domain=[-10,10];

var resolution=20;
var d3ObjId="svgElement1";
var interpolation='step-before';

var y = d3.scale.linear()
            .range([height-margin.bottom, margin.top])
            .domain(domain);

var yAxis = d3.svg.axis()
                .scale(y)
                .ticks(5)
                .orient("left")
                .tickSize(5,0,5);


var svg = d3.select("#"+d3ObjId)
            .append("svg")
                .attr("width", width)
                .attr("height", height);

svg.append("line")
    .attr("class", "boxplot")
    .attr("x1", margin.left)
    .attr("x2", width-margin.right)
    .attr("y1", y(0))
    .attr("y2", y(0));

for(var i=0; i<results.length; i++){
    results[i]=results[i].sort(d3.ascending)
    var g=svg.append("g").attr("transform", "translate("+(i*(boxWidth+boxSpacing)+margin.left)+",0)");
    addViolin(g, results[i], height, boxWidth, domain, 0.25, "#cccccc");
    addBoxPlot(g, results[i], height, boxWidth, domain, .15, "black", "white");

}

svg.append("g")
    .attr('class', 'axis')
    .attr("transform", "translate("+margin.left+",0)")
    .call(yAxis);

var d3ObjId="svgElement2";
var interpolation='basis';

var svg = d3.select("#"+d3ObjId)
            .append("svg")
                .attr("width", width)
                .attr("height", height);

svg.append("line")
    .attr("class", "boxplot")
    .attr("x1", margin.left)
    .attr("x2", width-margin.right)
    .attr("y1", y(0))
    .attr("y2", y(0));

for(var i=0; i<results.length; i++){
    results[i]=results[i].sort(d3.ascending)
    var g=svg.append("g").attr("transform", "translate("+(i*(boxWidth+boxSpacing)+margin.left)+",0)");
    addViolin(g, results[i], height, boxWidth, domain, 0.25, "#cccccc");

}

svg.append("g")
    .attr('class', 'axis')
    .attr("transform", "translate("+margin.left+",0)")
    .call(yAxis);
