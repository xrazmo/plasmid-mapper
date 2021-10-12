$(document).ready(function() {
    plotPlasmids()

    function plotPlasmids() {
        var size = 500;
        var radius = 150,
            innerRadius = radius - 10,
            outterRadius = radius + 10;
        var plasmid_len = 21500;

        var svg = d3.select('#main-svg');

        svg.attr("viewBox", "0 0 " + size + ' ' + size)
            .attr('xmlns', "http://www.w3.org/2000/svg")
            .attr('version', "1.1");
        var focus = svg.append('g').attr("transform", "translate(" + size / 2 + "," + size / 2 + ")");

        tick_values = d3.range(0, plasmid_len, 1e3)

        var x = d3.scaleBand()
            .range([0, 2 * Math.PI])
            .domain(tick_values);

        var y = d3.scaleRadial()
            .range([radius, radius + 2]) // Domain will be define later.
            .domain([0, 2]);

        // Add the labels
        // svg.append("g")
        //     .selectAll("g")
        //     .data(data)
        //     .enter()
        //     .append("g")
        //     .attr("text-anchor", function(d) { return (x(d.Country) + x.bandwidth() / 2 + Math.PI) % (2 * Math.PI) < Math.PI ? "end" : "start"; })
        //     .attr("transform", function(d) { return "rotate(" + ((x(d.Country) + x.bandwidth() / 2) * 180 / Math.PI - 90) + ")" + "translate(" + (y(d['Value']) + 10) + ",0)"; })
        //     .append("text")
        //     .text(function(d) { return (d.Country) })
        //     .attr("transform", function(d) { return (x(d.Country) + x.bandwidth() / 2 + Math.PI) % (2 * Math.PI) < Math.PI ? "rotate(180)" : "rotate(0)"; })
        //     .style("font-size", "11px")
        //     .attr("alignment-baseline", "middle")
        console.log(x.bandwidth());
        var xAxis = focus.append("g")
            .selectAll(".axis")
            .data(tick_values)
            .enter()
            .append("g")
            .attr("class", "axis")
            .attr("text-anchor", function(d) { return (x(d) + Math.PI) % (2 * Math.PI) < Math.PI ? "end" : "start"; })
            .attr("transform", function(d) { return "rotate(" + (x(d) * 180 / Math.PI - 90) + ")" + "translate(" + y(1.5) + ",0)"; });
        xAxis.append('line')
            .attr("x2", -4);
        xAxis.attr("stroke", "#bdbdbd")
            .append("text")
            .text(function(d) { return d + " bp" })
            .attr("transform", function(d) {
                var sign = (x(d) + x.bandwidth() / 2 + Math.PI) % (2 * Math.PI) < Math.PI ? 1 : -1;
                return "translate(0," + sign * 5 * (d.toString().length + 2) / 5 + ")rotate(90)"
            });


        focus.append("path")
            .attr("d", d3.arc()
                .innerRadius(radius)
                .outerRadius(radius + 0.5)
                .startAngle(0) // It's in radian, so Pi = 3.14 = bottom.
                .endAngle(2 * Math.PI) // 2*Pi = 6.28 = top
            ).attr('strok', '#bdbdbd');

        focus.append("path")
            .attr("d", d3.arc()
                .innerRadius(innerRadius)
                .outerRadius(innerRadius + 1.5)
                .startAngle(0 * Math.PI) // It's in radian, so Pi = 3.14 = bottom.
                .endAngle(1.5 * Math.PI) // 2*Pi = 6.28 = top
            ).attr('fill', '#69b3a2');

        focus.append("path")
            .attr("d", d3.arc()
                .innerRadius(outterRadius)
                .outerRadius(outterRadius - 1.5)
                .startAngle(0 * Math.PI) // It's in radian, so Pi = 3.14 = bottom.
                .endAngle(0.5 * Math.PI) // 2*Pi = 6.28 = top
            ).attr('fill', '#f903a2');

    }
});