$(document).ready(function() {


    var controls, touched, ringNr;
    const half_pi = Math.PI / 2.0;
    initForm()

    function initForm() {

        var qrySelect = document.getElementById('qryselect');
        $.each(Contig_ref, function(k, d) {
            qrySelect.options[qrySelect.options.length] = new Option(k, k);
        });
    }

    function initSVG(qryId) {
        var size = 700;
        var radius = 260,
            radiusStep = -20;

        controls = { 'radius': radius, 'radiusStep': radiusStep, "size": size };
        var svg = d3.select('#main-svg');

        svg.attr("viewBox", "0 0 " + size + ' ' + size)
            .attr('xmlns', "http://www.w3.org/2000/svg")
            .attr('version', "1.1");
        svg.append('g')
            .attr('id', 'focus')
            .attr("transform", "translate(" + size / 2 + "," + size / 2 + ")");

        plotPlasmid(Contig_ref[qryId], radius);
        var exR = radius;
        ringNr = 0
        $.each(MAP_DATA, function(key, data) {

            if (data.qseqid == qryId) {
                ringNr += 1
                exR += radiusStep;
                plotBlastRings(data, exR);

            }
        });


    }

    function plotBlastRings(data, radius) {
        var qLen = data.qlen;
        var focus = d3.select('#focus');
        var bl_focus = focus.append('g').attr('class', 'blast-focus');
        var coord2Angle = d3.scaleLinear().range([0, 2 * Math.PI]).domain([0, qLen])
        var arcW = 5,
            panelW = 20,
            p_inR = radius,
            p_outR = radius + panelW,
            innerR = p_inR + 2,
            outterR = innerR + arcW;

        bl_focus
            .append('path')
            .attr('class', ringNr % 2 == 0 ? 'e-ring' : 'o-ring')
            .attr('d', d3.arc()
                .innerRadius(p_inR)
                .outerRadius(p_outR)
                .startAngle(0)
                .endAngle(2 * Math.PI));


        $.each(data.ranges, function(i, rng) {

            bl_focus
                .append('path')
                .attr('class', 'sbj')
                .attr('d', d3.arc()
                    .innerRadius(innerR)
                    .outerRadius(outterR)
                    .startAngle(coord2Angle(rng.qstart))
                    .endAngle(coord2Angle(rng.qend)))
                .attr('fill', '#c6dbef')
                //d => "#" + Math.floor(Math.random() * 16777215).toString(16));
            bl_focus.selectAll('.miss_line-' + i)
                .data(rng.line_annot)
                .enter()
                .append('path')
                .attr('class', d => 'miss_line-' + i + ' ' + d.t)
                .attr('d', function(d) {

                    var angle = coord2Angle(d.v) - half_pi;
                    var x0 = innerR * Math.cos(angle),
                        y0 = innerR * Math.sin(angle),
                        x1 = outterR * Math.cos(angle),
                        y1 = outterR * Math.sin(angle);
                    return ["M", x0, y0, "L", x1, y1].join(' ')
                })




        });

        var header = bl_focus.append('text');
        bl_focus.append('path')
            .attr('id', 'hp' + data.sseqid)
            .attr('d', d3.arc()
                .innerRadius(p_outR - 10)
                .outerRadius(p_outR - 9.9)
                .startAngle(-Math.PI / 8)
                .endAngle(2 * Math.PI)).attr('fill', 'none');


        header.append("textPath")
            .attr("xlink:href", "#hp" + data.sseqid)
            .text(data.stitle)
            .attr("startOffset", "0%")
            .style('font', '6px sans-serif');


    }

    function plotPlasmid(data, radius) {

        var qLen = data.qlen;
        var focus = d3.select('#focus');
        var qryfocus = focus.append('g').attr('class', 'qry-focus');
        var stick_values = d3.range(0, qLen, 15e3)
        var coord2Angle = d3.scaleLinear().range([0, 2 * Math.PI]).domain([0, qLen]);
        var x = d3.scaleBand()
            .range([0, 2 * Math.PI])
            .domain(d3.range(0, qLen));

        var y = d3.scaleRadial()
            .range([radius - 5, radius - 1]) // Domain will be define later.
            .domain([0, 2]);


        var xAxis = qryfocus.append("g")
            .selectAll(".axis")
            .data(stick_values)
            .enter()
            .append("g")
            .attr("class", "axis")
            .attr("text-anchor", function(d) { return (x(d) + Math.PI) % (2 * Math.PI) < Math.PI ? "end" : "start"; })
            .attr("transform", function(d) { return "rotate(" + (x(d) * 180 / Math.PI - 90) + ")" + "translate(" + y(0) + ",0)"; });
        xAxis.append('line')
            .attr("x2", 6);

        xAxis.attr("stroke", "#bdbdbd")
            .append("text")
            .text(function(d) { return (d / 1000.0).toFixed(0) + " kb" })
            .attr("transform", function(d) {
                var sign = (x(d) + x.bandwidth() / 2 + Math.PI) % (2 * Math.PI) < Math.PI ? 1 : -1;
                return "translate(0," + sign * 5 * (d.toString().length) / 5 + ")rotate(90)"
            });


        qryfocus.append("path")
            .attr("d", d3.arc()
                .innerRadius(radius)
                .outerRadius(radius + 0.1)
                .startAngle(0) // It's in radian, so Pi = 3.14 = bottom.
                .endAngle(2 * Math.PI) // 2*Pi = 6.28 = top
            ).style('stroke', '#969696')
            .style('stroke-width', '0.3');



        textg = qryfocus.append('g');

        // textg.append('text')
        //     .attr('x', (radius / 4) * Math.cos(Math.PI))
        //     .attr('y', (radius / 3) * Math.cos(Math.PI))
        //     .text(data.accession)
        //     .style("font-size", "1rem")
        //     .style('font-weight', 600);

        var half_pi = Math.PI / 2,
            pi2 = 2 * Math.PI,
            orfR = [radius + 3, radius + 8],
            orfLblR = radius - 20;
        var secondRadius,
            tcoord2Angle;
        var sColor = '#dd3497';
        var recR = [orfR[0] - 10, radius + 15];
        secondRadius = recR[0] + 30;

        $.each(data.annotations, function(i, d) {

            qryfocus.append("path")
                .attr("d", d3.arc()
                    .innerRadius(recR[0])
                    .outerRadius(recR[1])
                    .startAngle(coord2Angle(d.sidx))
                    .endAngle(coord2Angle(d.eidx))
                ).style('fill', 'none')
                .style('stroke', sColor)
                .style('stroke-width', '0.5')
                .style("stroke-dasharray", ("2,1"));

            // #------------------
            var qryLen = Math.abs(d.sidx - d.eidx);
            var arcSidx = d.sidx - Math.min(2 * qryLen, qLen / 15),
                arcEidx = d.eidx + Math.min(2 * qryLen, qLen / 15);

            var ticks = d3.range(d.sidx, d.eidx, 1e3)
            tcoord2Angle = d3.scaleLinear().range([coord2Angle(arcSidx) % pi2, coord2Angle(arcEidx) % pi2]).domain([d.sidx, d.eidx]);
            var tx = d3.scaleBand()
                .range([coord2Angle(arcSidx), coord2Angle(arcEidx)])
                .domain(d3.range(d.sidx, d.eidx));

            var ty = d3.scaleRadial()
                .range([secondRadius - 2, secondRadius]) // Domain will be define later.
                .domain([0, 2]);

            var txAxis = qryfocus.append("g");

            var ticks = txAxis.selectAll(".taxis")
                .data(ticks)
                .enter()
                .append("g")
                .attr("class", "taxis")
                .attr("transform", function(d) { return "rotate(" + (tx(d) * 180 / Math.PI - 90) + ")" + "translate(" + ty(2) + ",0)"; });
            ticks.append('line')
                .attr("x2", -2).style('stroke', sColor).style('stroke-width', 0.5);

            txAxis.append("path")
                .attr("d", getUnaligned_deletion(recR[1], secondRadius, coord2Angle(d.sidx + (qryLen / 2)), coord2Angle(arcSidx), coord2Angle(arcEidx)))
                .style('stroke', sColor)
                .style("stroke-dasharray", ("1,1"))
                .style('stroke-width', '0.5')
                .attr('fill', 'None');

        });

        $.each(data.orfs, function(i, d) {

            var isAnn = false;
            $.each(data.annotations, function(i, ann) {
                if (d.eidx >= ann.sidx && d.sidx <= ann.eidx) {
                    isAnn = true;
                    return 1;
                }
            });

            qryfocus.append("path")
                .attr('class', "orf " + d.type)
                .attr("d", getArrowedArc(orfR[0], orfR[1], coord2Angle(d.sidx),
                    coord2Angle(d.eidx), d.strand == 1))
                .style('stroke', '#737373')
                .style('stroke-width', 0.3);

            // Also plot the ORF on the second/inner/zoomed axis    
            if (isAnn) {
                qryfocus.append("path")
                    .attr('class', "orf " + d.type)
                    .attr("d", getArrowedArc(secondRadius + 2, secondRadius + 8, tcoord2Angle(d.sidx),
                        tcoord2Angle(d.eidx), d.strand == 1))
                    .style('stroke', '#737373')
                    .style('stroke-width', 0.3);
            }
            if (d.type == 'hypothetical') return;

            if (isAnn) {
                orfLblR = secondRadius + 10
                textg.append('path')
                    .attr('id', 'line-' + d.id)
                    .attr("d", getORFLables(secondRadius + 8, secondRadius + 18,
                        tcoord2Angle(d.sidx), tcoord2Angle(d.eidx)))
                    .style('stroke', '#000')
                    .style("stroke-dasharray", ("1,1"))
                    .style('stroke-width', '0.1')
                    .style('fill', 'none');


                var midPoint = d.sidx + Math.abs(d.sidx - d.eidx),
                    x = (secondRadius + 18) * Math.cos(tcoord2Angle(midPoint) - half_pi),
                    y = (secondRadius + 18) * Math.sin(tcoord2Angle(midPoint) - half_pi);
                textg.append('g').append('text')
                    .attr('x', x)
                    .attr('y', y)
                    .attr('transform', 'rotate(0,' + x + ',' + y + ')')
                    .style("font-size", "0.35rem")
                    .style('font-weight', 600).style('font-style', 'italic')
                    .text(d.dscr.replace('family transposase', ''))
                    .on("mousedown", function(event) {
                        event.preventDefault();

                        this.style.cursor = "grabbing";
                        touched = true;
                        d3.select(this).style('font-size', '0.7rem');
                    })
                    .on('mouseleave mouseup', function(event) {
                        touched = false; // signals mouse up for (D) and (E)
                        this.style.cursor = "grab";
                        d3.select(this).style('font-size', '0.35rem');

                    })
                    .on("mousemove", function(event) {
                        event.preventDefault();
                        if (!touched) return; // mousemove with the mouse up

                        var t = d3.pointer(event),
                            x1 = t[0],
                            y1 = t[1];

                        var line = qryfocus.select('#line-' + d.id);
                        var sp = line.attr("d").split(" ");
                        var x, y, x0 = sp[1],
                            y0 = sp[2];
                        var txt = $(this).text();
                        x = x1, y = y1;
                        var bias = 0; // txt.length < 20 ? 0 : Math.max(5, 2 * txt.length);
                        if (x0 > x1) {
                            x = x + bias;
                        }
                        if (y1 > y0) {
                            y = y + 1;
                        }
                        sp[sp.length - 2] = x
                        sp[sp.length - 1] = y

                        line.attr('d', sp.join(" "));
                        $(this).attr('x', x1 - 5)
                            .attr('y', y1 + 2);

                        var tr = d3.select(this).attr("transform");
                        pp = tr.replace('rotate(', '').replace(');', '').split(',')
                        $(this).attr('transform', 'rotate(' + pp[0] +
                            ',' + (x1 - 5) + ',' + (y1 + 2) + ')');

                    }).on('mousewheel', function(event) {

                        event.preventDefault();
                        var sig = event.wheelDelta > 0 ? 1 : -1;
                        var x1 = parseInt($(this).attr('x')),
                            y1 = parseInt($(this).attr('y'));

                        var line = qryfocus.select('#line-' + d.id);
                        var sp = line.attr("d").split(" ");

                        var tr = d3.select(this).attr("transform");
                        pp = tr.replace('rotate(', '').replace(');', '').split(',')
                        $(this).attr('transform', 'rotate(' + (parseInt(pp[0]) + sig * 5) +
                            ',' + x1 + ',' + y1 + ')');
                    });

            }

        });
    }


    function getORFLables(innerRadius, outerRadius, startAngle, endAngle) {


        startAngle = startAngle - half_pi;
        endAngle = endAngle - half_pi;

        midAngle = startAngle + Math.abs(startAngle - endAngle) / 2;
        midRadius = innerRadius + 3;
        var x0 = innerRadius * Math.cos(midAngle),
            y0 = innerRadius * Math.sin(midAngle),
            x1 = midRadius * Math.cos(midAngle),
            y1 = midRadius * Math.sin(midAngle),
            x2 = outerRadius * Math.cos(midAngle),
            y2 = outerRadius * Math.sin(midAngle);

        var d = ["M", x0, y0, "L", x1, y1, "L", x2, y2]

        return d.join(' ');

    }

    function getArrowedArc(innerRadius, outerRadius, startAngle, endAngle, strand) {



        startAngle = startAngle - half_pi;
        endAngle = endAngle - half_pi;

        var deltaAngle = Math.abs(endAngle - startAngle) % (2 * Math.PI),
            ar = Math.abs(Math.min(0.02, 0.3 * deltaAngle)),
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

    function getUnaligned_deletion(innerRadius, outerRadius, baseAngle, startAngle, endAngle) {



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



        return d.join(' ');

    }

    $("#qryselect").on('change', function() {
        d3.select("#main-svg").selectAll('*').remove();
        controls = initSVG(this.value);
    });
    $('#qryselect').val('p004KP_6').change();


});