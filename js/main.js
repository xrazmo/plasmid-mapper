$(document).ready(function() {

    var columns = ["accession", "location", "organism", "nat_host", "isolation_source", "plasmid_name", "strain", "country", "collection_date", "sequence"]
    var freezedColumns = []
    var isolate_tablular_data = Object.keys(alignments_data).map(function(k) {
        table_row = {}
        $.each(columns, function(idx, col) {
            table_row[col] = alignments_data[k][col]
        });
        table_row['id'] = k;
        return table_row;
    });



    tabulate(isolate_tablular_data, columns)
    drawSequences();

    function tabulate(data, columns) {
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
            alignmentBox(alignID, data);
        });
    }

    function alignmentBox(alignID, data) {

        var margin = { top: 10, right: 10, bottom: 1, left: 10 };

        var width = 2000,
            height = 50;

        var alignmentLength = 21e3;
        var x = d3.scaleLinear().domain([1, alignmentLength]).range([0, width]),
            y = d3.scaleLinear().range([height, 0]);


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
        // console.log(data);
        focus.selectAll('.alignbox')
            .data(data)
            .enter()
            .append("rect")
            .attr('class', 'alignbox')
            .attr('x', function(d) { return x(d.qry_index[0]); })
            .attr('y', 0)
            .attr('width', function(d) { return x(Math.abs(d.qry_index[0] - d.qry_index[1])); })
            .attr('height', 15);

        var seq_comp = []
        $.each(data, function(k, arr) {
            seq_comp = seq_comp.concat(arr['line_annot']);
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


    }

});