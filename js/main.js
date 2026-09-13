$(document).ready(function () {
  tag = document.getElementsByTagName("title");
  var columns;
  var freezedColumns = [];
  var CONTEXT_DATA;
  var IS_MAINPAGE;
  var ALIGNMENT_LENGTH;

  $("#spinner").addClass("busy");

  $('[data-toggle="popover"]').popover();
  plotColorguid();

  var sample_name = {};

  $.searchParams = function (name) {
    var results = new RegExp("[?&]" + name + "=([^&#]*)").exec(
      window.location.href
    );
    if (!results) {
      return null;
    } else {
      return decodeURI(results[1]) || 0;
    }
  };

  catchParam();

  function catchParam() {
    var paramName = ["id"];
    paramName.forEach(function (param) {
      value = $.searchParams(param);
      if (value) {
        sample_name[param] = value;
      }
    });
  }

  if ($(tag).attr("id") != "main-page") {
    columns = [
      "#",
      "accession",
      "definition",
      "BioSample",
      "host",
      "isolation_source",
      "plasmid",
      "strain",
      "country",
      "collection_date",
      "sequence",
    ];
    CONTEXT_DATA = alignments_data[sample_name["id"]];
    if (typeof CONTEXT_DATA == "undefined") {
      $("#main-section").html(
        "<p class='h3 text-center text-danger'>No alignment to show!</p>"
      );
    }
    IS_MAINPAGE = false;
    ALIGNMENT_LENGTH = CONTEXT_DATA["qlen"];

    $(tag).text(sample_name["id"]);
    $("#header").html("Contig: <mark>" + sample_name["id"] + "</mark>");
  } else {
    columns = [
      "#",
      "alignments",
      "accession",
      "organism",
      "comments",
      "sequence",
    ];
    CONTEXT_DATA = reference_contexts;
    ALIGNMENT_LENGTH = 22000;
    IS_MAINPAGE = true;
  }

  tabulate(columns);
  drawSequences();
  $("#spinner").addClass("busy");

  function tabulate(columns) {
    var data = [];

    $.each(CONTEXT_DATA, function (k, cnt) {
      var row = {};

      $.each(columns, function (idx, col) {
        row[col] = CONTEXT_DATA[k][col];
      });
      row["id"] = k;
      if ((k.startsWith("aln") || IS_MAINPAGE) && !k.startsWith("orf")) {
        data.push(row);
      }
      if (k.startsWith("ref")) {
        data = [row].concat(data);
      }
    });

    var table = d3.select("#tbl-main");
    var thead = table.append("thead");
    var tbody = table.append("tbody");

    // append the header row
    thead
      .append("tr")
      .selectAll("th")
      .data(columns)
      .enter()
      .append("th")
      .append("button")
      .attr("class", "btn btn-light")
      .attr("data-toggle", "button")
      .text(function (column) {
        return column;
      })
      .on("click", function () {
        var col_name = $(this).text();
        var idx = freezedColumns.indexOf(col_name);

        if (idx > -1) {
          columnTh = $("table th:contains('" + col_name + "')");
          columnIndex = columnTh.index() + 1;
          $("table tr td:nth-child(" + columnIndex + ")")
            .css("position", "")
            .css("color", "")
            .css("background", "");
          columnTh.css("position", "");
          freezedColumns.splice(idx, 1);
        } else {
          freezedColumns.push(col_name);
          var leftMargin = 2;
          $.each(columns, function (i, col) {
            if (freezedColumns.indexOf(col) > -1) {
              columnTh = $("table th:contains('" + col + "')");
              columnIndex = columnTh.index() + 1;

              $("table tr td:nth-child(" + columnIndex + ")")
                .css("position", "sticky")
                .css("left", leftMargin + "px")
                .css("color", "#2c7fb8")
                .css("background", "#cccccc47");
              columnTh.css("position", "sticky").css("left", leftMargin + "px");
              leftMargin += columnTh.width() + 5;
            }
          });
        }
      });

    // create a row for each object in the data

    var rows = tbody.selectAll("tr").data(data).enter().append("tr");

    // create a cell in each row for each column
    var rownr = 0;
    rows
      .selectAll("td")
      .data(function (row) {
        return columns.map(function (column) {
          if (column == "sequence") {
            return { column: "sequence", id: row["id"] };
          }
          if (column == "alignments") {
            return { column: "alignments", acc: row["accession"] };
          }
          return { column: column, value: row[column] };
        });
      })
      .enter()
      .append("td")
      .attr("class", function (d) {
        return d.column == "sequence" ? "seqviewer" : undefined;
      })
      .attr("id", function (d) {
        return d.column == "sequence" ? d.id : undefined;
      })
      .html(function (d) {
        var colW = 21;
        if (d.column == "#") {
          rownr = rownr + 1;
          return rownr;
        } else if (d.column == "accession") {
          return rownr > 1 && !IS_MAINPAGE
            ? "<a href='https://www.ncbi.nlm.nih.gov/nuccore/" +
                d.value +
                "' target='_blank'>" +
                d.value +
                "</a>"
            : d.value;
        } else if (d.column == "BioSample") {
          return d.value
            ? "<a href='https://www.ncbi.nlm.nih.gov/biosample/" +
                d.value +
                "' target='_blank'>" +
                d.value +
                "</a>"
            : "";
        } else if (d.column == "alignments") {
          return (
            '<a class="btn btn-link" href="./html/aln.html?id=' +
            d.acc +
            '" target="_blank"><i class="fa fa-external-link"></i></a>'
          );
        } else if (d.column == "sequence") {
          return undefined;
        } else if (d.column == "definition" || d.column == "comments") {
          colW = 50;
        }
        return chunkSubstr(d.value, colW);
      });
  }

  function drawSequences() {
    $(".seqviewer").map(function () {
      var alignID = $(this).attr("id");
      var data = CONTEXT_DATA[alignID]["ranges"];
      // console.log(alignID, data);
      var isRef = (alignID == "ref" ? true : false) || IS_MAINPAGE;
      alignmentBox(alignID, data, CONTEXT_DATA[alignID]["accession"], isRef);
    });
  }

  function alignmentBox(alignID, data, accession, isRef) {
    var margin = { top: 6, right: 10, bottom: 1, left: 10 };

    var width = 2000,
      height = 60;

    var x = d3.scaleLinear().domain([1, ALIGNMENT_LENGTH]).range([0, width]),
      y = d3.scaleLinear().domain([1, 3]).range([height, 0]);

    var xAxis = d3
      .axisTop(x)
      .tickSize(height - 2)
      .tickFormat((t) => t + " bp");
    var dummy_xAxis = d3.axisTop(x).tickSize(0).tickValues([]);
    var yAxis = d3.axisLeft(y).tickSize(0).tickValues([]);

    var svg = d3
      .select("#" + alignID)
      .append("svg")
      .attr("width", width)
      .attr("height", height)
      .attr("xmlns", "http://www.w3.org/2000/svg")
      .attr("version", "1.1");

    var focus = svg
      .append("g")
      .attr("class", isRef ? "focus refcontext" : "focus")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
      .style("fill", "#c6dbef");

    var seq_comp = [],
      orflist = [];
    $.each(data, function (k, arr) {
      seq_comp = seq_comp.concat(arr["line_annot"]);
      orflist = orflist.concat(arr["orfs"]);
    });

    if (!isRef) {
      focus
        .selectAll(".alignbox")
        .data(data)
        .enter()
        .append("rect")
        .attr("class", "alignbox")
        .attr("accession", accession)
        .attr("sbjindex", function (d) {
          return d.sbj_index[0] + "-" + d.sbj_index[1];
        })
        .attr("x", function (d) {
          return x(d.qry_index[0]);
        })
        .attr("y", 5)
        .attr("width", function (d) {
          return x(Math.abs(d.qry_index[0] - d.qry_index[1]));
        })
        .attr("height", 15)
        .on("click", function (e) {
          var acc = $(this).attr("accession");
          var indexes = $(this).attr("sbjindex");
          var pp = indexes.split("-");
          pp = pp.map((x) => parseInt(x));

          var ncbiGraphicStr =
            "https://www.ncbi.nlm.nih.gov/nuccore/" +
            acc +
            "?report=graph&amp;from=" +
            (pp[0] - 1000) +
            "&amp;to=" +
            (pp[1] + 1000) +
            "&amp;mk=" +
            pp[0] +
            ":" +
            pp[1] +
            "|Aligned region|008000&amp";

          var win = window.open(ncbiGraphicStr, "_blank");
          if (win) {
            //Browser has allowed it to be opened
            win.focus();
          } else {
            //Browser has blocked it
            alert("Please allow popups for this website");
          }
        });

      // console.log(seq_comp);
      focus
        .selectAll(".alig-guid")
        .data(seq_comp)
        .enter()
        .append("line")
        .attr("class", function (d) {
          return "alig-guid " + d.t;
        })
        .attr("x1", function (d) {
          return x(d.v);
        })
        .attr("y1", 5)
        .attr("x2", function (d) {
          return x(d.v);
        })
        .attr("y2", 20)
        .style("stroke-width", 1);
    }

    focus
      .append("g")
      .attr("class", "axis x-axis")
      .attr("transform", "translate(0," + height + ")")
      .call(xAxis);

    focus
      .append("g")
      .attr("class", "axis x-axis")
      .attr("transform", "translate(0," + 0 + ")")
      .call(dummy_xAxis);

    focus.append("g").attr("class", "axis axis--y").call(yAxis);

    var xMax = 0;

    $.each(orflist, function (i, orf_dic) {
      var orfInfo = CONTEXT_DATA["orf_collection"][orf_dic["id"]];
      if (typeof orfInfo == "undefined" || !orfInfo) {
        orfInfo = {
          type: "unknown",
          dbname: "N/A",
          refprotien: "N/A",
          idty: "N/A",
          cov: "N/A",
          gap: "N/A",
          mismatch: "N/A",
          dscr: "unknown ",
        };
      }
      orflist[i] = Object.assign({}, orfInfo, orf_dic);

      xMax = Math.max(xMax, orf_dic["sidx"]);
      xMax = Math.max(xMax, orf_dic["eidx"]);
    });

    focus
      .selectAll(".orf")
      .data(orflist)
      .enter()
      .append("path")
      .attr("id", function (d) {
        return d.id;
      })
      .attr("class", function (d) {
        // console.log(d);
        return "orf " + d.type;
      })
      .attr("d", function (d) {
        return getPath(
          {
            x: x(d.sidx),
            y: y(isRef ? 2 : 1.7),
          },
          {
            x: x(d.eidx),
            y: y(isRef ? 2 : 1.7),
          },
          height / 3,
          height / 3,
          height / 3
        );
      })
      .attr("transform", function (d) {
        return getTransform(
          {
            x: x(d.sidx),
            y: y(isRef ? 2 : 1.7),
          },
          {
            x: x(d.eidx),
            y: y(isRef ? 2 : 1.7),
          }
        );
      })
      .on("mouseover", function (e) {
        // var newpopover = Mustache.render(BLASTX_POPOVER_TEMPLATE, d);
        orfId = $(this).attr("id");
        var d = CONTEXT_DATA["orf_collection"][orfId];
        if (!d) {
          d = {
            idty: "N/A",
            cov: "N/A",
            dscr: "unknown ",
          };
        }
        $(this).popover({
          placement: "auto",
          trigger: "hover",
          html: true,
          content: function () {
            return Mustache.render(BLASTX_POPOVER_TEMPLATE, d);
          },
        });
        $(this).popover("show");

        d3.select(this).transition().attr("style", "stroke-width:5px;");
      })
      .on("mouseleave", function (d, i) {
        d3.select(this)
          .transition()
          .delay(100)
          .attr("style", "stroke-width:0px;");
      });

    var xVisibleMax = x(xMax);
    var CHAR_SPACE = 4;

    focus
      .selectAll(".orfLbl")
      .data(orflist)
      .enter()
      .append("text")
      .attr("class", "orfLbl")
      .attr("transform", function (d) {
        return getTextTransform(
          d,
          isRef ? 1.9 : 1.6,
          d.dscr.length,
          xVisibleMax
        );
      })
      .attr("display", (d) => (textFits(d, xVisibleMax) ? null : "none"))
      .text((d) => trimText(d, xVisibleMax));

    function getTextTransform(d, ydt, textlen, maxAxis) {
      var leftPoint = Math.min(x(d.eidx), x(d.sidx));
      var rightPoint = Math.max(x(d.eidx), x(d.sidx));
      leftPoint = Math.max(0, leftPoint);
      rightPoint = Math.min(rightPoint, maxAxis);
      var xtr = leftPoint + Math.abs(rightPoint - leftPoint) / 2;
      return (
        "translate(" + (xtr - 0.5 * textlen * CHAR_SPACE) + "," + y(ydt) + ")"
      );
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
      var tt = visibleLength - d.dscr.length * CHAR_SPACE;
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
      var d = [
        "M" + -rx,
        "0a" + rx,
        ry + " 0 1",
        "0 " + 2 * rx,
        "0a" + rx,
        ry + " 0 1",
        "0 " + -2 * rx,
        "0",
      ];

      return d.join(",");
    } else {
      // The difference between the line width and the arrow width
      var dW = arrowheadWidth - lineWidth;
      // The angle of the line
      var angle = (Math.atan2(dy, dx) * 180) / Math.PI;

      // Generate a path describing the arrow. For simplicity we define it as a
      // horizontal line of the right length, and starting at 0,0. Then we rotate
      // and move it into place with a transform attribute.

      var d = [
        "M",
        0,
        -lineWidth / 2,
        "h",
        len - arrowheadLength,
        "v",
        -dW / 2,
        "L",
        len,
        0,
        "L",
        len - arrowheadLength,
        arrowheadWidth / 2,
        "v",
        -dW / 2,
        "H",
        0,
        "Z",
      ];
    }

    return d.join(" ");
  }

  function getTransform(from, to) {
    // rotate the arrow if it represent an ORF in reverse strand
    var dx = to.x - from.x;
    var dy = to.y - from.y;
    var angle = (Math.atan2(dy, dx) * 180) / Math.PI;
    var scaleTxt = angle === 180 ? "scale(-1, 1)" : "";

    return "translate(" + from.x + "," + from.y + ") " + scaleTxt;
  }

  function plotColorguid() {
    var svg = d3
      .select("#colorguid")
      .attr("viewBox", [0, 0, 100, 90])
      .attr("preserveAspectRatio", "xMidYMid meet");
    var orfWidth = 90,
      orfheight = 8;
    var margin = 5;
    var orfDic = {
      args: "Antibiotic resistance genes",
      biocidemetal: "Biocide and metal resistance genes",
      isel: "Insertion sequences",
      transposase: "Other transposases",
      integrase: "Integrase genes",
      virulence: "Virulance factors",
      other: "Other genes",
      unknown: "ORFs without annotation",
    };
    var row = 5;
    $.each(orfDic, function (klass, label) {
      var g = svg.append("g").attr("transform", "translate(5," + row + ")");
      g.append("path")
        .attr("class", klass)
        .attr(
          "d",
          getPath(
            { x: margin, y: 1 },
            { x: orfWidth + margin, y: 1 },
            orfheight,
            orfheight,
            orfheight
          )
        );

      g.append("text")
        .attr("class", "orfLbl")
        .attr("transform", "translate(" + (45 - label.length) + ",2)")
        .text(label)
        .style("font-size", "4px");

      row += orfheight + 2;
    });
  }

  function chunkSubstr(in_str, size) {
    if (typeof in_str == "undefined") {
      return in_str;
    }
    const numChunks = Math.ceil(in_str.length / size);
    //   const chunks = new Array(numChunks)
    var new_str = "";
    for (let i = 0, o = 0; i < numChunks; ++i, o += size) {
      // chunks[i] = str.substr(o, size)
      new_str = new_str + "<br>" + in_str.substr(o, size);
    }

    return new_str;
  }
});
