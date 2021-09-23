$(document).ready(function() {

    var columns = ["accession", "location", "organism", "nat_host", "isolation_source", "plasmid_name", "strain", "country", "collection_date", "sequence"]
    var freezedColumns = []

    tabulate(columns)
    drawSequences();

    function tabulate(columns) {
        var data = []
        $.each(alignments_data, function(k, cnt) {
            var row = {}
            $.each(columns, function(idx, col) {
                row[col] = alignments_data[k][col]
            });
            row['id'] = k;
            if (k.startsWith('aln')) {
                data.push(row)
            }

        });


        var table = d3.select('#tbl-main')
        var thead = table.append('thead')
        var tbody = table.append('tbody');

        // append the header row
        thead.append('tr')
            .selectAll('th')
            .data(columns).enter()
            .append('th')
            .append('button')
            .attr('class', 'btn btn-light')
            .attr('data-toggle', 'button')
            .text(function(column) { return column; }).on('click', function() {
                var col_name = $(this).text();
                var idx = freezedColumns.indexOf(col_name)

                if (idx > -1) {
                    columnTh = $("table th:contains('" + col_name + "')");
                    columnIndex = columnTh.index() + 1;
                    $('table tr td:nth-child(' + columnIndex + ')').css("position", "")
                        .css('color', '').css('background', '');;
                    columnTh.css("position", "");
                    freezedColumns.splice(idx, 1);
                } else {
                    freezedColumns.push(col_name);
                    var leftMargin = 2;
                    $.each(columns, function(i, col) {
                        if (freezedColumns.indexOf(col) > -1) {
                            columnTh = $("table th:contains('" + col + "')");
                            columnIndex = columnTh.index() + 1;

                            $('table tr td:nth-child(' + columnIndex + ')')
                                .css("position", "sticky").css("left", leftMargin + "px")
                                .css('color', '#2c7fb8').css('background', '#cccccc47');
                            columnTh.css("position", "sticky").css("left", leftMargin + "px");
                            leftMargin += columnTh.width() + 5;
                        }
                    });


                }


            });

        // create a row for each object in the data
        var rows = tbody.selectAll('tr')
            .data(data)
            .enter()
            .append('tr');

        // create a cell in each row for each column
        rows.selectAll('td')
            .data(function(row) {
                return columns.map(function(column) {
                    if (column == "sequence") {
                        return { column: "sequence", id: row['id'] }
                    }
                    return { column: column, value: row[column] };
                });
            })
            .enter()
            .append('td')
            .attr('class', function(d) { return d.column == 'sequence' ? 'seqviewer' : undefined; })
            .attr('id', function(d) { return d.column == 'sequence' ? d.id : undefined; })
            .text(function(d) { return d.column != 'sequence' ? d.value : undefined; });



    }

    function drawSequences() {

        $(".seqviewer").map(function() {
            var alignID = $(this).attr('id');
            var data = alignments_data[alignID]['ranges']
                // console.log(alignID, data);
            alignmentBox(alignID, data, alignments_data[alignID]['accession']);
        });
    }

    function alignmentBox(alignID, data, accession) {

        var margin = { top: 10, right: 10, bottom: 1, left: 10 };

        var width = 4000,
            height = 50;

        var alignmentLength = 21580;
        var x = d3.scaleLinear().domain([1, alignmentLength]).range([0, width]),
            y = d3.scaleLinear().domain([1, 3]).range([height, 0]);


        var xAxis = d3.axisTop(x).tickSize(height - 2).tickFormat(t => t + ' bp');
        var dummy_xAxis = d3.axisTop(x).tickSize(0).tickValues([]);
        var yAxis = d3.axisLeft(y).tickSize(0).tickValues([]);

        var svg = d3.select('#' + alignID)
            .append('svg')
            .attr('width', width)
            .attr('height', height)
            .attr('xmlns', "http://www.w3.org/2000/svg")
            .attr('version', "1.1");

        var focus = svg.append("g")
            .attr("class", "focus")
            .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

        focus.selectAll('.alignbox')
            .data(data)
            .enter()
            .append("rect")
            .attr('class', 'alignbox')
            .attr('accession', accession)
            .attr('sbjindex', function(d) { return d.sbj_index[0] + '-' + d.sbj_index[1]; })
            .attr('x', function(d) { return x(d.qry_index[0]); })
            .attr('y', 0)
            .attr('width', function(d) { return x(Math.abs(d.qry_index[0] - d.qry_index[1])); })
            .attr('height', 15)
            .on('click', function(e) {
                var acc = $(this).attr('accession');
                var indexes = $(this).attr('sbjindex');
                var pp = indexes.split('-');
                pp = pp.map(x => parseInt(x))

                var ncbiGraphicStr = "https://www.ncbi.nlm.nih.gov/nuccore/" + acc +
                    "?report=graph&amp;from=" + (pp[0] - 1000) +
                    "&amp;to=" + (pp[1] + 1000) + "&amp;mk=" + pp[0] + ":" + pp[1] +
                    "|Aligned region|008000&amp";

                var win = window.open(ncbiGraphicStr, '_blank');
                if (win) {
                    //Browser has allowed it to be opened
                    win.focus();
                } else {
                    //Browser has blocked it
                    alert('Please allow popups for this website');
                }
            });

        var seq_comp = [],
            orflist = [];
        $.each(data, function(k, arr) {
            seq_comp = seq_comp.concat(arr['line_annot']);
            orflist = orflist.concat(arr['orfs']);
        });


        // console.log(seq_comp);
        focus.selectAll('.alig-guid')
            .data(seq_comp)
            .enter()
            .append('line')
            .attr('class', function(d) { return 'alig-guid ' + d.t; })
            .attr('x1', function(d) { return x(d.v) })
            .attr('y1', 0)
            .attr('x2', function(d) { return x(d.v) })
            .attr('y2', 15)
            .style("stroke-width", 1);

        focus.append("g")
            .attr("class", "axis x-axis")
            .attr("transform", "translate(0," + height + ")")
            .call(xAxis);

        focus.append("g")
            .attr("class", "axis x-axis")
            .attr("transform", "translate(0," + 0 + ")").call(dummy_xAxis);


        focus.append("g")
            .attr("class", "axis axis--y")
            .call(yAxis);


        var xMax = 0;

        $.each(orflist, function(i, orf_dic) {
            var orfInfo = alignments_data['orf_collection'][orf_dic['id']]
            if (typeof orfInfo == "undefined" || !orfInfo) {

                orfInfo = {
                    "type": 'unknown',
                    "dbname": 'N/A',
                    "refprotien": 'N/A',
                    "idty": 'N/A',
                    "cov": 'N/A',
                    "gap": 'N/A',
                    "mismatch": 'N/A',
                    "dscr": "unknown "
                }
            }
            orflist[i] = Object.assign({}, orfInfo, orf_dic)

            xMax = Math.max(xMax, orf_dic['sidx'])
            xMax = Math.max(xMax, orf_dic['eidx'])

        });

        focus.selectAll(".orf")
            .data(orflist)
            .enter()
            .append("path")
            .attr('id', function(d) {
                return d.id
            })
            .attr('class', function(d) {
                // console.log(d);
                return "orf " + d.type;
            })
            .attr("d", function(d) {
                return getPath({
                    x: x(d.sidx),
                    y: y(2)
                }, {
                    x: x(d.eidx),
                    y: y(2)
                }, height / 3, height / 3, height / 3);
            })
            .attr("transform", function(d) {
                return getTransform({
                    x: x(d.sidx),
                    y: y(2)
                }, {
                    x: x(d.eidx),
                    y: y(2)
                });
            }).on('mouseover', function(e) {
                // var newpopover = Mustache.render(BLASTX_POPOVER_TEMPLATE, d);
                orfId = $(this).attr('id')
                var d = alignments_data['orf_collection'][orfId]
                if (!d) {
                    d = {
                        "idty": 'N/A',
                        "cov": 'N/A',
                        "dscr": "unknown "
                    }
                }
                $(this).popover({
                    placement: 'auto',
                    trigger: 'hover',
                    "html": true,
                    content: function() {

                        return Mustache.render(BLASTX_POPOVER_TEMPLATE, d);
                    }
                });
                $(this).popover('show');

                d3.select(this)
                    .transition()
                    .attr('style', 'stroke-width:5px;');
            })
            .on('mouseleave', function(d, i) {

                d3.select(this)
                    .transition().delay(100)
                    .attr('style', 'stroke-width:0px;');
            });


        var xVisibleMax = x(xMax);
        var CHAR_SPACE = 4;

        focus.selectAll(".orfLbl")
            .data(orflist).enter()
            .append('text')
            .attr('class', 'orfLbl')
            .attr("transform", function(d) {

                return getTextTransform(d, 1.85, d.dscr.length, xVisibleMax);
            })
            .attr('display', d => textFits(d, xVisibleMax) ? null : 'none').text(d => trimText(d, xVisibleMax));



        function getTextTransform(d, ydt, textlen, maxAxis) {
            var leftPoint = Math.min(x(d.eidx), x(d.sidx));
            var rightPoint = Math.max(x(d.eidx), x(d.sidx));
            leftPoint = Math.max(0, leftPoint);
            rightPoint = Math.min(rightPoint, maxAxis);
            var xtr = leftPoint + Math.abs(rightPoint - leftPoint) / 2
            return "translate(" + (xtr - (0.5 * textlen * CHAR_SPACE)) + "," + y(ydt) + ")";
        }

        function textFits(d, maxAxis) {
            var leftPoint = Math.max(0, Math.min(x(d.eidx), x(d.sidx)));
            var rightPoint = Math.min(Math.max(x(d.eidx), x(d.sidx)), maxAxis);
            var visibleLength = rightPoint - leftPoint;
            return d.dscr.length * CHAR_SPACE < visibleLength;
        }

        function trimText(d, maxAxis) {
            var leftPoint = Math.max(0, Math.min(x(d.eidx), x(d.sidx)));
            var rightPoint = Math.min(Math.max(x(d.eidx), x(d.sidx)), maxAxis);
            var visibleLength = rightPoint - leftPoint;
            var tt = visibleLength - d.dscr.length * CHAR_SPACE
            return d.dscr.substring(0, d.dscr.length + tt);
        }
    }

    function getPath(from, to, lineWidth, arrowheadWidth, arrowheadLength) {
        var dx = to.x - from.x;
        var dy = to.y - from.y;

        // Calculate the length of the line
        var len = Math.sqrt(dx * dx + dy * dy);

        if (len < arrowheadLength) {

            var rx = 0.5 * dx,
                ry = 1.5;
            // return 'M '+dx+' '+dy+' m -'+r+', 0 a '+r+','+r+' 0 1,0 '+(r*2)+',0 a '+r+','+r+' 0 1,0 -'+(r*2)+',0';
            var d = ['M' + (-rx), '0a' + rx, ry + " 0 1", "0 " + (2 * rx), '0a' + rx, ry + " 0 1", "0 " + (-2 * rx), "0"];

            return d.join(',');
        } else {

            // The difference between the line width and the arrow width
            var dW = arrowheadWidth - lineWidth;
            // The angle of the line
            var angle = Math.atan2(dy, dx) * 180 / Math.PI;

            // Generate a path describing the arrow. For simplicity we define it as a
            // horizontal line of the right length, and starting at 0,0. Then we rotate
            // and move it into place with a transform attribute.

            var d = ['M', 0, -lineWidth / 2,
                'h', len - arrowheadLength,
                'v', -dW / 2,
                'L', len, 0,
                'L', len - arrowheadLength, arrowheadWidth / 2,
                'v', -dW / 2,
                'H', 0,
                'Z'
            ];
        }


        return d.join(' ');

    }

    function getTransform(from, to) {
        // rotate the arrow if it represent an ORF in reverse strand
        var dx = to.x - from.x;
        var dy = to.y - from.y;
        var angle = Math.atan2(dy, dx) * 180 / Math.PI;
        var scaleTxt = angle === 180 ? 'scale(-1, 1)' : '';

        return "translate(" + from.x + "," + from.y + ") " + scaleTxt;
    }



});