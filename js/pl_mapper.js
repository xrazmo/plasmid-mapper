$(document).ready(function() {


    var controls, touched, ringNr, qryId, selected_alignments = [];
    const half_pi = Math.PI / 2.0;
    initForm()

    function initForm() {

        var qrySelect = document.getElementById('qryselect');
        $.each(Contig_ref, function(k, d) {
            qrySelect.options[qrySelect.options.length] = new Option(k, k);
        });
    }

    function update_page(qryId) {
        var size = 800;
        var radius = 260,
            radiusStep = -5;

        controls = { 'radius': radius, 'radiusStep': radiusStep, "size": size };
        var svg = d3.select('#main-svg');

        svg.attr("viewBox", "0 0 " + size + ' ' + size)
            .attr('xmlns', "http://www.w3.org/2000/svg")
            .attr('version', "1.1")
            .style('font', 'avenir next, sans-serif');
        svg.append('g')
            .attr('id', 'focus')
            .attr("transform", "translate(" + size / 2 + "," + size / 2 + ")");

        var data = [];
        var effectiveData = PlasmidMapperEdits.mergeEdits(qryId, Contig_ref[qryId]);
        plotPlasmid(effectiveData, radius);
        var columns = ["#", "select", "qcov", "sseqid", "stitle", 'qseqid']

        $.each(MAP_DATA, function(key, d) {

            if (d.qseqid == qryId) {
                var tmpDic = { id: key }
                columns.map(function(col) {
                    tmpDic[col] = d[col];
                });

                data.push(tmpDic);

            }
        });
        tabulate(data, columns);
        plotLegend(qryId);

        $('.big-checkbox').change(function() {
            var id = $(this).attr('id');
            if (this.checked) {
                selected_alignments.push(id);
            } else {
                index = selected_alignments.indexOf(id);
                if (index > -1) {
                    selected_alignments.splice(index, 1);
                }
            }

        });
        $.each($('.big-checkbox'), function(i, chk) {
            if (i < 40) {
                $(chk).prop("checked", true);
                $(chk).trigger('change');
            }
        });
        $('#uptBtn').trigger('click');
        return controls;

    }

    function plotLegend(qryId) {
        var deg = Math.PI / 180,
            pi2 = 2 * Math.PI;
        var legend = d3.select('#focus').append('g')
        var orf_labels = {
            'ARGs': { 'kl': 'args', 'coef': 1.5 },
            'Insertion sequences': { 'kl': 'isel', 'coef': 4.1 },
            "Transposons": { 'kl': "transposase", 'coef': 2.8 },
            'Virulence factors': { 'kl': 'virulence', 'coef': 3.5 },
            'Biocide and metal resistance': { 'kl': 'biocidemetal', 'coef': 5.8 },
            "Integron": { 'kl': "integrase", 'coef': 2 },
            "Hypothetical proteins": { 'kl': 'hypothetical', 'coef': 4.5 },
            "Other": { 'kl': "other", 'coef': 1.5 },

        }
        var lengendAngle = {
            "p004KP_6": 65 * deg,
            "p165E_3": 115 * deg,
            "s082Km_2": -80 * deg,
            "s164ECL_2": 140 * deg,
            "s257ECL_2": -125 * deg,
            "s304ECL_3": 40 * deg,
            "m481ECL_2": 140 * deg,
            "s202ECL_2": 140 * deg,
        }


        var radius = controls.radius,
            sAngle = lengendAngle[qryId] ? lengendAngle[qryId] : 0,
            k = 0,
            step = 2.7 * deg;
        var ta, tb, bias, sa, sb;

        var lgTxt = legend.append('text');

        var tmp = (sAngle + pi2) % pi2;
        var revert = tmp > half_pi && tmp < 2.5 * half_pi ? true : false;

        $.each(orf_labels, function(txt, d) {


            ta = sAngle + (k * step),
                tb = sAngle + (k + d.coef) * step,
                bias = (tb - step - ta) / 2
            sa = ta + bias;
            sb = sa + step;


            var pR1 = revert ? radius + 20 : radius + 27,
                pR2 = revert ? radius + 21 : radius + 28,
                oR1 = revert ? radius + 23 : radius + 20,
                oR2 = revert ? radius + 28 : radius + 25;
            if (revert) {
                tmp = ta;
                ta = tb;
                tb = tmp;
            }


            legend.append("path")
                .attr('class', "orf " + d.kl)
                .attr("d", getArrowedArc(oR1, oR2, sa, sb, true))
                .style('fill', ORF_COLOR[d.kl])
                .style('stroke', '#737373')
                .style('stroke-width', 0.3);

            legend.append("path")
                .attr('id', 'lgd-' + d.kl)
                .attr("d", d3.arc()
                    .innerRadius(pR1)
                    .outerRadius(pR2)
                    .startAngle(ta)
                    .endAngle(tb))
                .style('stroke', 'none').style('fill', 'none')

            lgTxt.append("textPath")
                .attr("xlink:href", "#lgd-" + d.kl)
                .text(txt)
                .attr("startOffset", "0%")
                .style('font-size', '5px')
                .style('font-weight', 600)
                .style('font-family', 'tahoma');
            k += d.coef
        });


        legend.append("path")
            .attr("d", d3.arc()
                .innerRadius(radius + 15)
                .outerRadius(radius + 37)
                .startAngle(sAngle - deg)
                .endAngle(Math.max(tb, ta) + deg))
            .style('stroke', '#bdbdbd')
            .style('fill', '#cccccc2b')
            .style('stroke-width', '0.5');

    }

    function plotBlastLegend(qryId) {
        var deg = Math.PI / 180,
            pi2 = 2 * Math.PI;

        var legend = d3.select('#focus').append('g')
        var lengendAngle = {
            "p004KP_6": -120 * deg,
            "p165E_3": -70 * deg,
            "s082Km_2": 95 * deg,
            "s164ECL_2": -45 * deg,
            "s257ECL_2": -45 * deg,
            "m481ECL_2": 220 * deg,
            "s202ECL_2": -45 * deg,
            "s304ECL_3": 210 * deg,
        }
        var radius = controls.radius,
            sAngle = lengendAngle[qryId] ? lengendAngle[qryId] : 0,
            k = 0,
            step = 4.5 * Math.PI / 180;
        var ta, tb, bias, sa, sb;
        var lgTxt = legend.append('text');
        var r1 = radius + 18,
            r2 = radius + 28,
            r, coef = 2;

        var tmp = (sAngle + pi2) % pi2;
        var revert = tmp > 0.5 * half_pi && tmp < 2 * half_pi ? true : false;
        $.each(selected_alignments, function(i, key) {

            r = i % 2 == 1 ? r1 : r2;

            var tR1 = revert ? r + -1 : r - 1,
                tR2 = revert ? r + 3 : r + 3,
                pR1 = revert ? r + 4 : r - 1,
                pR2 = revert ? r + 7.5 : r + 2;


            if (i % 2 == 0) {
                ta = sAngle + (k * step),
                    tb = sAngle + (k + coef) * step,
                    bias = (tb - step - ta) / 2
                sa = ta + bias;
                sb = sa + step;
                if (revert) {
                    tmp = ta;
                    ta = tb;
                    tb = tmp;
                }

            }
            legend.append("path")
                .attr("d", d3.arc()
                    .innerRadius(pR1)
                    .outerRadius(pR2)
                    .startAngle(sa)
                    .endAngle(sb))
                .style('stroke', '#ccc')
                .style('stroke-width', 0.5)
                .style('fill', Color_collection[i])

            legend.append("path")
                .attr('id', 'lgdB-' + key)
                .attr("d", d3.arc()
                    .innerRadius(tR1)
                    .outerRadius(tR2)
                    .startAngle(ta)
                    .endAngle(tb))
                .style('stroke', 'none')
                .style('fill', 'none');

            lgTxt.append("textPath")
                .attr("xlink:href", "#lgdB-" + key)
                .text(key.split('$')[1])
                .attr("startOffset", "0%")
                .style('font-size', '5px')
                .style('font-weight', 600)
                .style('font-family', 'tahoma');
            k += 1
        });

        legend.append("path")
            .attr("d", d3.arc()
                .innerRadius(radius + 15)
                .outerRadius(radius + 37)
                .startAngle(sAngle - deg)
                .endAngle(Math.max(ta, tb) + deg))
            .style('stroke', '#bdbdbd')
            .style('fill', '#cccccc2b')
            .style('stroke-width', '0.5');

    }

    function plotBlastRings_withHeader(data, radius) {
        var qLen = data.qlen;
        var bl_focus = d3.select('#bl-focus');
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
                .endAngle(2 * Math.PI))
            .style('fill', ringNr % 2 == 0 ? "#f0f0f0" : "none")
            .style('stroke', ringNr % 2 == 0 ? '#ccc' : 'none')
            .style('stroke-width', 1)
            .style("stroke-dasharray", ("5,4"));


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
            bl_focus.selectAll('.miss_line-' + i + '-' + outterR)
                .data(rng.line_annot)
                .enter()
                .append('path')
                .attr('class', d => 'miss_line-' + i + '-' + outterR + ' ' + d.t)
                .attr('d', function(d) {

                    var angle = coord2Angle(d.v) - half_pi;
                    var x0 = innerR * Math.cos(angle),
                        y0 = innerR * Math.sin(angle),
                        x1 = outterR * Math.cos(angle),
                        y1 = outterR * Math.sin(angle);
                    return ["M", x0, y0, "L", x1, y1].join(' ')
                }).style('stroke-width', 0.2)
                .style('stroke', d => Mismatch_COLOR[d.t]);

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
            .style('font-size', '8px')
            .style('font-family', 'monospace')
            .style('font-weight', 'bold');

    }

    function plotBlastRings(data, radius, color_index) {
        var qLen = data.qlen;
        var bl_focus = d3.select('#bl-focus');
        var coord2Angle = d3.scaleLinear().range([0, 2 * Math.PI]).domain([0, qLen])
        var arcW = 4,
            panelW = 4,
            p_inR = radius,
            p_outR = radius + panelW,
            innerR = p_inR + 2,
            outterR = innerR + arcW;

        var cl = Color_collection[color_index % Color_collection.length]

        $.each(data.ranges, function(i, rng) {

            bl_focus
                .append('path')
                .attr('class', 'sbj')
                .attr('d', d3.arc()
                    .innerRadius(innerR)
                    .outerRadius(outterR)
                    .startAngle(coord2Angle(rng.qstart))
                    .endAngle(coord2Angle(rng.qend)))
                .style('fill', cl + 'cc')
                .style('stroke', '#ccc')
                .style('stroke-width', 0.7)
                //d => "#" + Math.floor(Math.random() * 16777215).toString(16)); '#c6dbef'
            bl_focus.selectAll('.miss_line-' + i + '-' + outterR)
                .data(rng.line_annot)
                .enter()
                .append('path')
                .attr('class', d => 'miss_line-' + i + '-' + outterR + ' ' + d.t)
                .attr('d', function(d) {

                    var angle = coord2Angle(d.v) - half_pi;
                    var x0 = innerR * Math.cos(angle),
                        y0 = innerR * Math.sin(angle),
                        x1 = outterR * Math.cos(angle),
                        y1 = outterR * Math.sin(angle);
                    return ["M", x0, y0, "L", x1, y1].join(' ')
                }).style('stroke-width', 0.2)
                .style('stroke', d => Mismatch_COLOR[d.t]);

        });

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
            .range([radius - 8, radius]) // Domain will be define later.
            .domain([0, 2]);


        var xAxis = qryfocus.append("g")
            .selectAll(".axis")
            .data(stick_values)
            .enter()
            .append("g")
            .attr("class", "axis")
            .attr("text-anchor", function(d) { return (x(d) + Math.PI) % (2 * Math.PI) < Math.PI ? "end" : "start"; })
            .attr("transform", function(d) { return "rotate(" + (x(d) * 180 / Math.PI - 90) + ")" + "translate(" + y(0) + ",0)"; })
            .style('stroke', '#000')
            .style('stroke-width', '0.2px')
            .style('font-size', '7px')
            .style('font-weight', 600)
            .style('font-family', 'sans-serif');

        xAxis.append('line')
            .attr("x2", 8);

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
                .style('stroke-width', 1)
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
                .attr("x2", -2).style('stroke', sColor)
                .style('stroke-width', 0.5);

            txAxis.append("path")
                .attr("d", getUnaligned_deletion(recR[1], secondRadius, coord2Angle(d.sidx + (qryLen / 2)), coord2Angle(arcSidx), coord2Angle(arcEidx)))
                .style('stroke', sColor)
                .style("stroke-dasharray", ("1,1"))
                .style('stroke-width', 1)
                .attr('fill', 'none');

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
                .style('fill', ORF_COLOR[d.type])
                .style('stroke', '#737373')
                .style('stroke-width', 0.3);

            // Also plot the ORF on the second/inner/zoomed axis    
            if (isAnn) {
                qryfocus.append("path")
                    .attr('class', "orf " + d.type)
                    .attr('id', 'orf-' + d.id)
                    .attr("d", getArrowedArc(secondRadius + 2, secondRadius + 8, tcoord2Angle(d.sidx),
                        tcoord2Angle(d.eidx), d.strand == 1))
                    .style('fill', ORF_COLOR[d.type])
                    .style('stroke', '#737373')
                    .style('stroke-width', 0.3).on('click', function(event) {
                        var orfid = d3.select(this).attr('id').split('-')[1];
                        var curStat = d3.select('#line-' + orfid).attr('display')
                        d3.select('#line-' + orfid).attr('display', curStat == 'none' ? 'block' : 'none');
                        d3.select('#txt-' + orfid).attr('display', curStat == 'none' ? 'block' : 'none');
                    });
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
                var labelOverride = d._labelOverride;
                var initRotation = 0;
                if (labelOverride && typeof labelOverride.x === 'number') {
                    x = labelOverride.x;
                    y = labelOverride.y;
                    initRotation = labelOverride.rotation || 0;
                }
                textg.append('g')
                    .append('text')
                    .attr('id', 'txt-' + d.id)
                    .attr('x', x)
                    .attr('y', y)
                    .attr('transform', 'rotate(' + initRotation + ',' + x + ',' + y + ')')
                    .style("font-size", "8px")
                    .style('font-weight', 600)
                    .style('font-family', 'Helvetica')
                    .text((labelOverride && labelOverride.text) || d.dscr.replace('family transposase', ''))
                    .on("dblclick", function(event) {
                        event.preventDefault();
                        openLabelTextEditor(this, qryId, d.id, $(this).text());
                    })
                    .on("mousedown", function(event) {
                        event.preventDefault();

                        this.style.cursor = "grabbing";
                        touched = true;
                        d3.select(this).style('font-size', '10px');
                    })
                    .on('mouseleave mouseup', function(event) {
                        touched = false; // signals mouse up for (D) and (E)
                        this.style.cursor = "grab";
                        d3.select(this).style('font-size', '8px');

                        var finalX = parseFloat($(this).attr('x')),
                            finalY = parseFloat($(this).attr('y'));
                        var tr = d3.select(this).attr("transform");
                        var rotation = parseFloat(tr.replace('rotate(', '').split(',')[0]) || 0;
                        PlasmidMapperEdits.setLabelPosition(qryId, 'orf-' + d.id, finalX, finalY, rotation);
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
                        var newRotation = parseInt(pp[0]) + sig * 5;
                        $(this).attr('transform', 'rotate(' + newRotation +
                            ',' + x1 + ',' + y1 + ')');
                        PlasmidMapperEdits.setLabelPosition(qryId, 'orf-' + d.id, x1, y1, newRotation);
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

    function tabulate(data, columns) {

        var table = d3.select('#tbl-main')
        var thead = table.append('thead')
        var tbody = table.append('tbody');

        // append the header row
        thead.append('tr')
            .selectAll('th')
            .data(columns).enter()
            .append('th')
            .text(function(column) { return column; });

        // create a row for each object in the data

        var rows = tbody.selectAll('tr')
            .data(data)
            .enter()
            .append('tr');

        // create a cell in each row for each column
        var rownr = 0
        rows.selectAll('td')
            .data(function(row) {

                return columns.map(function(column) {
                    if (column == "select") {
                        return { column: "select", id: row['id'] }
                    }
                    return { column: column, value: row[column] };
                });
            })
            .enter()
            .append('td')
            .html(function(d) {

                var colW = 80;
                if (d.column == '#') {
                    rownr = rownr + 1
                    return rownr;
                } else if (d.column == 'sseqid') { return "<a href='https://www.ncbi.nlm.nih.gov/nuccore/" + d.value + "' target='_blank'>" + d.value + "</a>" } else if (d.column == 'select') { return "<input type=\"checkbox\" class=\"big-checkbox\" id=\"" + d.id + "\"></input>" }

                return chunkSubstr(d.value, colW);
            });

    }

    function chunkSubstr(in_str, size) {
        if (typeof in_str == 'undefined' || typeof in_str == 'number') {
            return in_str
        }
        const numChunks = Math.ceil(in_str.length / size)
            //   const chunks = new Array(numChunks)
        var new_str = "";
        for (let i = 0, o = 0; i < numChunks; ++i, o += size) {
            // chunks[i] = str.substr(o, size)
            new_str = new_str + "<br>" + in_str.substr(o, size)

        }

        return new_str;
    }

    function clear_canvas() {
        var canvas = document.getElementById("canvas");
        var ctx = canvas.getContext("2d");
        ctx.clearRect(0, 0, canvas.width, canvas.height);
    }

    $("#qryselect").on('change', function() {
        d3.select("#main-svg").selectAll('*').remove();
        d3.select("#tbl-main").selectAll('*').remove();
        clear_canvas();
        selected_alignments = [];
        qryId = this.value;
        controls = update_page(qryId);

    });



    $('#uptBtn').on('click', function(event) {
        var focus = d3.select('#focus')
        focus.selectAll('.blast-focus').remove();
        focus.append('g')
            .attr('class', 'blast-focus')
            .attr('id', 'bl-focus');

        if (selected_alignments.length > 0) {

            var exR = controls.radius - 10
            ringNr = 0

            $.each(selected_alignments, function(i, key) {
                ringNr += 1
                exR += controls.radiusStep;
                plotBlastRings(MAP_DATA[key], exR, i);

            });

            plotBlastLegend(qryId);

        }

    });


    $('#genBtn').on('click', function(event) {
        // saveSVG("main-svg", "saveLink");
        clear_canvas();
        var svgString = new XMLSerializer().serializeToString(document.querySelector('svg'));
        var canvas = document.getElementById("canvas");
        var ctx = canvas.getContext("2d");
        var DOMURL = self.URL || self.webkitURL || self;
        var img = new Image();
        var svg = new Blob([svgString], { type: "image/svg+xml;charset=utf-8" });
        var url = DOMURL.createObjectURL(svg);
        img.onload = function() {
            ctx.drawImage(img, 0, 0);
            var png = canvas.toDataURL("image/tiff");
            saveAs(png, qryId + ".tiff");
            DOMURL.revokeObjectURL(png);
        };
        img.src = url;
    });


    $('.inpt').on('change', function() {
        var id = $(this).attr('id');
        if (id == 'width-inpt') {
            $('#canvas').attr('width', $(this).val() + 'px')

        } else {
            $('#canvas').attr('height', $(this).val() + 'px')

        }
    });

    var qrySelectEl = document.getElementById('qryselect');
    var defaultQryId = qrySelectEl.options.length > 1 ? qrySelectEl.options[1].value : qrySelectEl.options[0].value;
    $('#qryselect').val(defaultQryId).change();
});