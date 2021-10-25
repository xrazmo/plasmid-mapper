$(document).ready(function() {
    plotPlasmids()

    function plotPlasmids() {
        var size = 500;
        var radius = 150,
            innerRadius = radius - 10,
            outterRadius = radius + 10;
        var plasmid_len = 25000;

        var svg = d3.select('#main-svg');

        svg.attr("viewBox", "0 0 " + size + ' ' + size)
            .attr('xmlns', "http://www.w3.org/2000/svg")
            .attr('version', "1.1");
        var focus = svg.append('g').attr("transform", "translate(" + size / 2 + "," + size / 2 + ")");

        var tick_values = d3.range(0, plasmid_len, 1e3)

        var coord2Angle = d3.scaleLinear().range([0, 2 * Math.PI]).domain([0, plasmid_len]);

        var x = d3.scaleBand()
            .range([0, 2 * Math.PI])
            .domain(d3.range(0, plasmid_len));

        var y = d3.scaleRadial()
            .range([radius, radius + 2]) // Domain will be define later.
            .domain([0, 2]);

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
                .outerRadius(innerRadius + 5)
                .startAngle(coord2Angle(500)) // It's in radian, so Pi = 3.14 = bottom.
                .endAngle(coord2Angle(20000)) // 2*Pi = 6.28 = top
            ).attr('fill', '#69b3a2');

        focus.append("path")
            .attr("d", getORF(outterRadius - 5, outterRadius, coord2Angle(15232), coord2Angle(20545), -1)) // 2*Pi = 6.28 = top  d3.arc()
            .attr('fill', '#f903a2');

        console.log(x(0) + ", " + x(2000));

    }

    function getORF(innerRadius, outerRadius, startAngle, endAngle, strand) {

        var half_pi = Math.PI / 2.0;

        startAngle = startAngle - half_pi;
        endAngle = endAngle - half_pi;

        var deltaAngle = Math.abs(endAngle - startAngle),
            ar = Math.abs(Math.min(0.02, 0.1 * deltaAngle)),
            arrowAngle = endAngle;

        var lgflag = deltaAngle > Math.PI ? 1 : 0;
        endAngle = endAngle - ar;

        if (strand == -1) {
            arrowAngle = startAngle
            startAngle = startAngle + ar;
        }
        var midR = innerRadius + Math.abs(innerRadius - outerRadius) / 2;

        var x0 = innerRadius * Math.cos(startAngle),
            y0 = innerRadius * Math.sin(startAngle),
            x1 = innerRadius * Math.cos(endAngle),
            y1 = innerRadius * Math.sin(endAngle),
            x2 = outerRadius * Math.cos(endAngle),
            y2 = outerRadius * Math.sin(endAngle),
            x3 = outerRadius * Math.cos(startAngle),
            y3 = outerRadius * Math.sin(startAngle),
            xA = midR * Math.cos(arrowAngle),
            yA = midR * Math.sin(arrowAngle);

        var d = ["M", x0, y0,
            "A", innerRadius, innerRadius, 1, lgflag, 1, x1, y1,
            "L", xA, yA,
            "L", x2, y2,
            "A", outerRadius, outerRadius, 1, lgflag, 0, x3, y3,
            "Z"
        ];
        if (strand == -1) {
            console.log('reverse');
            d = ["M", x2, y2,
                "A", outerRadius, outerRadius, 1, lgflag, 0, x3, y3,
                "L", xA, yA,
                "L", x0, y0,
                "A", innerRadius, innerRadius, 1, lgflag, 1, x1, y1,
                "Z"
            ];
        }
        return d.join(' ');

    }

});