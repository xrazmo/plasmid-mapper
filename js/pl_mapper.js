$(document).ready(function() {
    drawSequences()


    function drawSequences() {

        // $(".seqviewer").map(function() {
        //     var alignID = $(this).attr('id');
        //     var data = MAP_DATA[alignID]['ranges']
        //         // console.log(alignID, data);
        //     var isRef = (alignID == 'ref' ? true : false) || IS_MAINPAGE;
        //     plotPlasmids(alignID, data,
        //         MAP_DATA[alignID]['accession'], isRef);
        // });
        var alignID = 'aln-1'
        plotPlasmids(MAP_DATA, alignID);
    }

    function plotPlasmids(data, alignID) {
        var size = 500;
        var radius = 120,
            innerRadius = radius - 10,
            outterRadius = radius + 10;
        var plasmid_len = MAP_DATA.slen;

        var svg = d3.select('#main-svg');

        svg.attr("viewBox", "0 0 " + size + ' ' + size)
            .attr('xmlns', "http://www.w3.org/2000/svg")
            .attr('version', "1.1");
        var focus = svg.append('g').attr("transform", "translate(" + size / 2 + "," + size / 2 + ")");

        //--------AXIS-----------
        var tick_values = d3.range(0, plasmid_len, 15e3)
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
            .text(function(d) { return (d / 1000.0).toFixed(0) + " kb" })
            .attr("transform", function(d) {
                var sign = (x(d) + x.bandwidth() / 2 + Math.PI) % (2 * Math.PI) < Math.PI ? 1 : -1;
                return "translate(0," + sign * 5 * (d.toString().length) / 5 + ")rotate(90)"
            });


        focus.append("path")
            .attr("d", d3.arc()
                .innerRadius(radius)
                .outerRadius(radius + 0.5)
                .startAngle(0) // It's in radian, so Pi = 3.14 = bottom.
                .endAngle(2 * Math.PI) // 2*Pi = 6.28 = top
            ).attr('strok', '#bdbdbd');

        // -------------------------------

        // background
        focus.append("path")
            .attr("d", d3.arc()
                .innerRadius(outterRadius)
                .outerRadius(outterRadius + 50)
                .startAngle(0) // It's in radian, so Pi = 3.14 = bottom.
                .endAngle(2 * Math.PI) // 2*Pi = 6.28 = top
            ).attr('fill', '#f7fbff');

        var alignRad1 = outterRadius + 2,
            alignRad2 = outterRadius + 5;

        focus.selectAll('.alignbox')
            .data(MAP_DATA[alignID]["ranges"])
            .enter()
            .append("path")
            .attr("d", function(d) {
                return getArrowedArc(alignRad1, alignRad2, coord2Angle(d.sbj_index[0]), coord2Angle(d.sbj_index[1]), d.sstrand == 'Plus' ? 1 : -1)
            })
            .attr('fill', function(d) { return '#' + Math.floor(Math.random() * 16777215).toString(16) });



        $.each(MAP_DATA[alignID]["una_reg"], function(i, reg) {
            var unaligR1 = outterRadius + 13,
                unaligR2 = outterRadius + 14;

            console.log(reg);
            var arcSidx = reg.ssidx,
                arcEidx = reg.seidx,
                baseIndex = reg.seidx,
                isDeletion = false;
            if (Math.abs(reg.qsidx - reg.qeidx) < 1000) {
                return;
            }
            if (Math.abs(reg.ssidx - reg.seidx) < 1000) {
                isDeletion = true;
                qryLen = Math.abs(reg.qsidx - reg.qeidx)
                arcSidx = reg.ssidx - Math.min(2 * qryLen, plasmid_len / 8);
                arcEidx = reg.seidx + Math.min(2 * qryLen, plasmid_len / 8);
                baseIndex = reg.seidx;
                if (qryLen > 1e5) {
                    unaligR2 += 10;
                    unaligR2 += 10;
                }

            }
            var ticks = d3.range(reg.qsidx, reg.qeidx, 1e3)
            var tx = d3.scaleBand()
                .range([coord2Angle(arcSidx), coord2Angle(arcEidx)])
                .domain(d3.range(reg.qsidx, reg.qeidx));

            var ty = d3.scaleRadial()
                .range([unaligR1, unaligR2]) // Domain will be define later.
                .domain([0, 2]);

            var txAxis = focus.append("g");

            var ticks = txAxis.selectAll(".taxis")
                .data(ticks)
                .enter()
                .append("g")
                .attr("class", "taxis")
                // .attr("text-anchor", function(d) { return (tx(d) + Math.PI) % (2 * Math.PI) < Math.PI ? "end" : "start"; })
                .attr("transform", function(d) { return "rotate(" + (tx(d) * 180 / Math.PI - 90) + ")" + "translate(" + ty(2) + ",0)"; });
            ticks.append('line')
                .attr("x2", -2);

            txAxis.append("path")
                .attr("d", isDeletion ? getUnaligned_deletion(alignRad1, unaligR2, coord2Angle(baseIndex), coord2Angle(arcSidx), coord2Angle(arcEidx)) :
                    getUnaligned_insertion(alignRad1, unaligR2, coord2Angle(arcSidx), coord2Angle(arcEidx)))
                .style('stroke', '#252525')
                .style("stroke-dasharray", ("1,1"))
                .style('stroke-width', '0.3')
                .attr('fill', 'None');

            txAxis.append("path")
                .attr("id", "unaligned_" + i) //Unique id of the path
                .attr("d", getUnaligned_insertion(alignRad1, unaligR2 - 7, coord2Angle(arcSidx), coord2Angle(arcEidx)))
                .style("fill", "none");

            txAxis.append("text")
                .append("textPath")
                .attr("xlink:href", "#unaligned_" + i)
                .attr("startOffset", "60%")
                .text((Math.abs(reg.qsidx - reg.qeidx) / 1000.0).toFixed(1) + " kb")
                .style("text-anchor", "middle")
                .style("font-size", "0.3rem")
                .style("font-weight", "600")



        });



        // focus.append("path")
        //     .attr("d", getORF(outterRadius - 5, outterRadius, coord2Angle(15232), coord2Angle(20545), -1)) // 2*Pi = 6.28 = top  d3.arc()
        //     .attr('fill', '#f903a2');


    }

    function getArrowedArc(innerRadius, outerRadius, startAngle, endAngle, strand) {

        var half_pi = Math.PI / 2.0;

        startAngle = startAngle - half_pi;
        endAngle = endAngle - half_pi;

        var deltaAngle = Math.abs(endAngle - startAngle) % 2 * Math.PI,
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

    function getUnaligned_insertion(innerRadius, outerRadius, startAngle, endAngle) {

        var half_pi = Math.PI / 2.0;

        startAngle = startAngle - half_pi;
        endAngle = endAngle - half_pi;
        var deltaAngle = Math.abs(endAngle - startAngle);
        var lgflag = deltaAngle > Math.PI ? 1 : 0;

        var x0 = innerRadius * Math.cos(startAngle),
            y0 = innerRadius * Math.sin(startAngle),
            x1 = outerRadius * Math.cos(startAngle),
            y1 = outerRadius * Math.sin(startAngle),
            x2 = outerRadius * Math.cos(endAngle),
            y2 = outerRadius * Math.sin(endAngle),
            x3 = innerRadius * Math.cos(endAngle),
            y3 = innerRadius * Math.sin(endAngle);

        var d = ["M", x0, y0,
            "L", x1, y1,
            "A", outerRadius, outerRadius, 1, lgflag, 1, x2, y2,
            "L", x3, y3,
        ];
        if (startAngle < half_pi && startAngle > -half_pi) {
            d = ["M", x3, y3,
                "L", x2, y2,
                "A", outerRadius, outerRadius, 1, lgflag, 0, x1, y1,
                "L", x0, y0,
            ];
        }

        return d.join(' ');

    }

    function getUnaligned_deletion(innerRadius, outerRadius, baseAngle, startAngle, endAngle) {

        var half_pi = Math.PI / 2.0;

        startAngle = startAngle - half_pi;
        endAngle = endAngle - half_pi;
        baseAngle = baseAngle - half_pi;
        var deltaAngle = Math.abs(endAngle - startAngle);

        var lgflag = deltaAngle > Math.PI ? 1 : 0;

        var x0 = innerRadius * Math.cos(baseAngle),
            y0 = innerRadius * Math.sin(baseAngle),
            x1 = outerRadius * Math.cos(baseAngle),
            y1 = outerRadius * Math.sin(baseAngle),
            x2 = outerRadius * Math.cos(startAngle),
            y2 = outerRadius * Math.sin(startAngle),
            x3 = outerRadius * Math.cos(endAngle),
            y3 = outerRadius * Math.sin(endAngle);

        var d = ["M", x2, y2,
            "A", outerRadius, outerRadius, 1, lgflag, 1, x1, y1,
            "L", x0, y0,
            "L", x1, y1,
            "A", outerRadius, outerRadius, 1, lgflag, 1, x3, y3,
        ];


        // if (baseAngle < half_pi) {
        //     var d = ["M", x3, y3,
        //         "A", outerRadius, outerRadius, 1, lgflag, 0, x1, y1,
        //         "L", x0, y0,
        //         "L", x1, y1,
        //         "A", outerRadius, outerRadius, 1, lgflag, 0, x2, y2,
        //     ];

        // }

        return d.join(' ');

    }

});