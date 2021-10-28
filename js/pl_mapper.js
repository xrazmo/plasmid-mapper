$(document).ready(function() {


    var controls = initSVG();


    function initSVG() {
        var size = 500;
        var radius = 120,
            sRadius = radius - 10,
            qRadius = radius + 50;
        var slen = MAP_DATA.slen,
            qlen = MAP_DATA.slen;


        var svg = d3.select('#main-svg');

        svg.attr("viewBox", "0 0 " + size + ' ' + size)
            .attr('xmlns', "http://www.w3.org/2000/svg")
            .attr('version', "1.1");
        svg.append('g')
            .attr('id', 'focus')
            .attr("transform", "translate(" + size / 2 + "," + size / 2 + ")");


        plotPlasmid(qRadius, qlen, 'qry');
        plotPlasmid(sRadius, slen, 'sbj');

        return { 'sRadius': sRadius, 'qRadius': qRadius, "size": size }
        // $.each(MAP_DATA["una_reg"], function(i, reg) {
        //     var unaligR1 = outterRadius + 13,
        //         unaligR2 = outterRadius + 14;

        //     console.log(reg);
        //     var arcSidx = reg.ssidx,
        //         arcEidx = reg.seidx,
        //         baseIndex = reg.seidx,
        //         isDeletion = false;
        //     if (Math.abs(reg.qsidx - reg.qeidx) < 1000) {
        //         return;
        //     }
        //     if (Math.abs(reg.ssidx - reg.seidx) < 1000) {
        //         isDeletion = true;
        //         qryLen = Math.abs(reg.qsidx - reg.qeidx)
        //         arcSidx = reg.ssidx - Math.min(2 * qryLen, subject_len / 8);
        //         arcEidx = reg.seidx + Math.min(2 * qryLen, subject_len / 8);
        //         baseIndex = reg.seidx;
        //         if (qryLen > 1e5) {
        //             unaligR2 += 10;
        //             unaligR2 += 10;
        //         }

        //     }
        //     var ticks = d3.range(reg.qsidx, reg.qeidx, 1e3)
        //     var tx = d3.scaleBand()
        //         .range([coord2Angle(arcSidx), coord2Angle(arcEidx)])
        //         .domain(d3.range(reg.qsidx, reg.qeidx));

        //     var ty = d3.scaleRadial()
        //         .range([unaligR1, unaligR2]) // Domain will be define later.
        //         .domain([0, 2]);

        //     var txAxis = focus.append("g");

        //     var ticks = txAxis.selectAll(".taxis")
        //         .data(ticks)
        //         .enter()
        //         .append("g")
        //         .attr("class", "taxis")
        //         // .attr("text-anchor", function(d) { return (tx(d) + Math.PI) % (2 * Math.PI) < Math.PI ? "end" : "start"; })
        //         .attr("transform", function(d) { return "rotate(" + (tx(d) * 180 / Math.PI - 90) + ")" + "translate(" + ty(2) + ",0)"; });
        //     ticks.append('line')
        //         .attr("x2", -2);

        //     txAxis.append("path")
        //         .attr("d", isDeletion ? getUnaligned_deletion(alignRad1, unaligR2, coord2Angle(baseIndex), coord2Angle(arcSidx), coord2Angle(arcEidx)) :
        //             getUnaligned_insertion(alignRad1, unaligR2, coord2Angle(arcSidx), coord2Angle(arcEidx)))
        //         .style('stroke', '#252525')
        //         .style("stroke-dasharray", ("1,1"))
        //         .style('stroke-width', '0.3')
        //         .attr('fill', 'None');

        //     txAxis.append("path")
        //         .attr("id", "unaligned_" + i) //Unique id of the path
        //         .attr("d", getUnaligned_insertion(alignRad1, unaligR2 - 7, coord2Angle(arcSidx), coord2Angle(arcEidx)))
        //         .style("fill", "none");

        //     txAxis.append("text")
        //         .append("textPath")
        //         .attr("xlink:href", "#unaligned_" + i)
        //         .attr("startOffset", "60%")
        //         .text((Math.abs(reg.qsidx - reg.qeidx) / 1000.0).toFixed(1) + " kb")
        //         .style("text-anchor", "middle")
        //         .style("font-size", "0.3rem")
        //         .style("font-weight", "600")



        // });



        // focus.append("path")
        //     .attr("d", getORF(outterRadius - 5, outterRadius, coord2Angle(15232), coord2Angle(20545), -1)) // 2*Pi = 6.28 = top  d3.arc()
        //     .attr('fill', '#f903a2');


    }

    function plotPlasmid(radius, pLen, type) {

        var focus = d3.select('#focus');
        var subFocus = focus.append('g').attr('class', type == 'sbj' ? 's-focus' : 'q-focus');
        var stick_values = d3.range(0, pLen, 15e3)
        var coord2Angle = d3.scaleLinear().range([0, 2 * Math.PI]).domain([0, pLen]);
        var x = d3.scaleBand()
            .range([0, 2 * Math.PI])
            .domain(d3.range(0, pLen));

        var y = d3.scaleRadial()
            .range([radius, radius + 2]) // Domain will be define later.
            .domain([0, 2]);

        var xAxis = subFocus.append("g")
            .selectAll(".axis")
            .data(stick_values)
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


        subFocus.append("path")
            .attr("d", d3.arc()
                .innerRadius(radius)
                .outerRadius(radius + 0.5)
                .startAngle(0) // It's in radian, so Pi = 3.14 = bottom.
                .endAngle(2 * Math.PI) // 2*Pi = 6.28 = top
            ).attr('strok', '#bdbdbd');

        subFocus.selectAll(type == 'sbj' ? '.s-alignbox' : '.q-alignbox')
            .data(MAP_DATA["ranges"])
            .enter()
            .append("path")
            .attr('id', d => (type == 'sbj' ? 's-' : 'q-') + d.id)
            .attr('class', type == 'sbj' ? 'alignbox s-alignbox' : 'alignbox q-alignbox')
            .attr("d", function(d) {

                return type == 'sbj' ? getArrowedArc(radius + 12, radius + 20, coord2Angle(d.ssidx), coord2Angle(d.seidx), true) :
                    getArrowedArc(radius - 10, radius - 2, coord2Angle(d.qsidx), coord2Angle(d.qeidx), d.strand);
            })
            .attr('fill', function() { return '#' + Math.floor(Math.random() * 16777215).toString(16) })
            .on('mouseover', function() {
                d3.selectAll('.alignbox').attr('opacity', '10%');
                d3.select(this).attr('opacity', '100%');
                var prefix = type == 'sbj' ? "#q-" : "#s-"
                var sid = prefix + this.id.split('-')[1];
                d3.select(sid).attr('opacity', '100%')

            })
            .on('mouseout', function() {

                d3.selectAll('.alignbox').attr('opacity', '100%');

            });
    }


    function getArrowedArc(innerRadius, outerRadius, startAngle, endAngle, strand) {

        var half_pi = Math.PI / 2.0;

        startAngle = startAngle - half_pi;
        endAngle = endAngle - half_pi;

        var deltaAngle = Math.abs(endAngle - startAngle) % (2 * Math.PI),
            ar = Math.abs(Math.min(0.02, 0.1 * deltaAngle)),
            arrowAngle = endAngle;

        if (startAngle > endAngle) {
            strand = !strand;
            t = endAngle
            endAngle = startAngle;
            startAngle = t

        }

        var lgflag = deltaAngle > Math.PI ? 1 : 0;

        endAngle = endAngle - ar;

        if (!strand) {
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
        if (!strand) {

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


    $('.form-check-input').change(function() {
        type = $(this).prop('value');
        if (type == 'qry') {
            var qlen = MAP_DATA['qlen']
            $.each(MAP_DATA['ranges'], function(idx, val) {

                val.qsidx = qlen - val.qsidx;
                val.qeidx = qlen - val.qeidx;

            });
            d3.selectAll('.q-focus').selectAll('*').remove();
            plotPlasmid(controls.qRadius, qlen, 'qry')


        } else if (type = 'sbj') {

        }
    });

});