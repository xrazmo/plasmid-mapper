$(document).ready(function() {


    var controls, touched, ringNr, qryId, selected_alignments = [];
    const half_pi = Math.PI / 2.0;
    // User-adjustable rendering settings (diameter, per-ring thickness,
    // ring spacing), persisted across re-renders unlike the old
    // radius/radiusStep locals recreated fresh on every update_page() call.
    // Read by update_page()/plotBlastRings() instead of the previous
    // hardcoded literals; the Settings UI controls write into this object
    // and trigger rerenderCurrentPlasmid().
    var renderSettings = {
        radius: 200,
        radiusStep: -14,
        ringThickness: 10,
        orfLegendStyle: 'curved', // 'curved' | 'rect'
        blastLegendStyle: 'curved',
        // Per-element font settings, read by every text-drawing call site
        // below instead of a hardcoded literal. Seeded from the fonts each
        // element used before this became configurable. The Style modal
        // (js/pl_style_modal.js) writes into this object and calls
        // rerenderCurrentPlasmid(), the same mechanism the Settings sidebar
        // controls already use.
        fonts: {
            axisLabel: { family: 'monospace', size: 8 }, // ring tick / mb axis text
            orfLabel: { family: 'Helvetica', size: 8 }, // on-figure ORF short labels
            legendText: { family: 'Helvetica', size: 9 }, // ORF-category legend
            blastLegendText: { family: 'sans-serif', size: 7 } // BLAST-subject legend
        },
        // Per-ORF-category fill colors, seeded from ORF_COLOR (js/inline_style.js)
        // so nothing changes visually until the Style modal is used to
        // override one. bandHighlight is the new zoom-band background fill.
        colors: Object.assign({}, ORF_COLOR, { bandHighlight: '#ffe58f' })
    };
    window.PlasmidMapperSettings = renderSettings;
    // Shared string keys for PlasmidMapperEdits.setLegendPosition/
    // getLegendPosition, defined once to avoid typo drift between
    // plotLegend's and plotBlastLegend's call sites.
    var ORF_LEGEND_KEY = 'orfLegend';
    var BLAST_LEGEND_KEY = 'blastLegend';
    initForm()

    function initForm() {

        var qrySelect = document.getElementById('qryselect');
        $.each(Contig_ref, function(k, d) {
            qrySelect.options[qrySelect.options.length] = new Option(k, k);
        });
    }

    function update_page(qryId) {
        var size = 800;
        var radius = renderSettings.radius,
            radiusStep = renderSettings.radiusStep;

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
        renderAnnotationControls(qryId, effectiveData.annotations);
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
        plotLegend(qryId, effectiveData);

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

    // Scans candidate angles (every 15 degrees) and returns the first
    // whose [angle, angle + requiredArcWidth] window doesn't overlap any
    // entry in occupiedRanges (each {start, end}, radians, normalized into
    // [0, 2*PI) by the caller). Falls back to the candidate with the
    // smallest total overlap if every candidate collides with something --
    // this can happen with several user-added zoom bands plus both
    // legends already occupying angular space -- rather than ever leaving
    // an angle undefined (which would regress to the pre-Stage-0 NaN
    // legend bug this builds on top of).
    function findClearAngle(occupiedRanges, requiredArcWidth) {
        var pi2 = 2 * Math.PI;
        var normalize = function(a) { return ((a % pi2) + pi2) % pi2; };

        function windowOverlap(startA, endA) {
            var total = 0;
            occupiedRanges.forEach(function(range) {
                var s = normalize(range.start), e = normalize(range.end);
                if (e < s) e += pi2; // occupied range straddles the 0/2*PI seam
                // Test the candidate window against the occupied range,
                // and again shifted by +/- 2*PI so a candidate window that
                // itself straddles the seam is checked correctly too.
                [-pi2, 0, pi2].forEach(function(shift) {
                    var os = Math.max(startA, s + shift);
                    var oe = Math.min(endA, e + shift);
                    if (oe > os) total += (oe - os);
                });
            });
            return total;
        }

        var bestAngle = 0, bestOverlap = Infinity;
        for (var deg = 0; deg < 360; deg += 15) {
            var candidate = deg * Math.PI / 180;
            var overlap = windowOverlap(candidate, candidate + requiredArcWidth);
            if (overlap <= 0) return candidate;
            if (overlap < bestOverlap) {
                bestOverlap = overlap;
                bestAngle = candidate;
            }
        }
        return bestAngle; // starvation fallback: smallest-overlap candidate
    }

    // Annotation bands are drawn as a dashed highlight box near the main
    // query ring (recR in plotPlasmid, roughly radius-7 to radius+15) --
    // well inside where the legends are drawn (radius+15 to radius+37), so
    // only bands would actually visually intersect a legend at all; this
    // helper converts data.annotations into the {start, end} ranges
    // findClearAngle expects.
    function annotationOccupiedRanges(data, coord2Angle) {
        return (data.annotations || []).map(function(ann) {
            return { start: coord2Angle(ann.sidx), end: coord2Angle(ann.eidx) };
        });
    }

    // Makes a curved (arc-following) legend draggable AROUND THE RING only
    // -- angle, not free x/y -- since the legend's whole visual identity
    // is text following the ring at a fixed radius, and letting it drift
    // off that ring would look wrong. Persists the angle (radians) via
    // PlasmidMapperEdits.
    //
    // Unlike a simple rotate() transform on pre-built content, this calls
    // back into `drawAt(angle)` (the same function used for the initial
    // render) on every drag tick to fully rebuild the legend's geometry
    // at the live angle -- required because `drawAt` internally decides
    // whether to flip the text/swatch radii ("revert") so text reads
    // right-side-up, and that decision depends on which half of the
    // circle the angle falls in. A rotate()-only transform bakes that
    // decision in once at the pre-drag angle, so dragging into the
    // opposite hemisphere left the text upside-down (confirmed bug) --
    // recomputing via drawAt() on every tick keeps it correctly oriented
    // throughout the drag, not just after it settles.
    //
    // The grabbed point (e.g. one specific category's swatch) is generally
    // NOT at the legend's angle-0 origin -- it's offset by wherever that
    // category falls in the arc -- so the angle must track the CHANGE in
    // pointer angle from drag-start, added to the angle the legend already
    // had, not the pointer's absolute angle each tick (which would
    // snap/jump the legend's origin to the cursor on the very first tick).
    function makeLegendDraggable(legend, qryId, legendKey, initialAngle, drawAt) {
        var currentAngle = initialAngle;
        var angleAtDragStart, pointerAngleAtDragStart;

        legend.style('cursor', 'grab');

        legend.call(d3.drag()
            .container(d3.select('#focus').node())
            .on('start', function(event) {
                this.style.cursor = 'grabbing';
                angleAtDragStart = currentAngle;
                pointerAngleAtDragStart = Math.atan2(event.y, event.x);
            })
            .on('drag', function(event) {
                var pointerAngleNow = Math.atan2(event.y, event.x);
                currentAngle = angleAtDragStart + (pointerAngleNow - pointerAngleAtDragStart);
                drawAt(currentAngle);
            })
            .on('end', function(event) {
                this.style.cursor = 'grab';
                PlasmidMapperEdits.setLegendPosition(qryId, legendKey, null, null, currentAngle * 180 / Math.PI);
            }));
    }

    // Rectangular/"boxed" legend alternative to the curved-arc style:
    // plain <rect>+<text> rows, centered in the figure's empty middle by
    // default (draggable via makeLegendDraggable, since the center hole
    // shrinks as more BLAST comparison rings are selected and can
    // genuinely collide with a naively-centered box). rows is an array of
    // {swatchColor, label}; a saved manual position fully replaces the
    // default centered placement (no angle to preserve, unlike curved
    // mode's findClearAngle).
    function plotLegendRect(qryId, legendKey, rows, boxWidth, fontKey) {
        var font = renderSettings.fonts[fontKey] || renderSettings.fonts.legendText;
        var rowHeight = 14, swatchSize = 10, padding = 8;
        var boxHeight = rows.length * rowHeight + padding * 2;

        var legend = d3.select('#focus').append('g');

        legend.append('rect')
            .attr('x', 0).attr('y', 0)
            .attr('width', boxWidth).attr('height', boxHeight)
            .style('stroke', '#bdbdbd')
            .style('fill', '#cccccc2b')
            .style('stroke-width', '0.5');

        var row = legend.selectAll('.legend-rect-row')
            .data(rows)
            .enter()
            .append('g')
            .attr('transform', function(d, i) { return 'translate(' + padding + ',' + (padding + i * rowHeight) + ')'; });

        row.append('rect')
            .attr('width', swatchSize).attr('height', swatchSize)
            .style('fill', function(d) { return d.swatchColor; })
            .style('stroke', '#737373')
            .style('stroke-width', 0.5);

        row.append('text')
            .attr('x', swatchSize + 5)
            .attr('y', swatchSize - 1)
            .style('font-size', font.size + 'px')
            .style('font-family', font.family)
            .text(function(d) { return d.label; });

        var saved = PlasmidMapperEdits.getLegendPosition(qryId, legendKey);
        var defaultX = -boxWidth / 2, defaultY = -boxHeight / 2;
        var startX = saved ? saved.x : defaultX, startY = saved ? saved.y : defaultY;

        var currentTransform = { x: startX, y: startY };
        var dragOrigin, pointerOrigin;
        legend.attr('transform', 'translate(' + currentTransform.x + ',' + currentTransform.y + ')')
            .style('cursor', 'grab');

        legend.call(d3.drag()
            .container(d3.select('#focus').node())
            .on('start', function(event) {
                this.style.cursor = 'grabbing';
                dragOrigin = { x: currentTransform.x, y: currentTransform.y };
                pointerOrigin = { x: event.x, y: event.y };
            })
            .on('drag', function(event) {
                currentTransform.x = dragOrigin.x + (event.x - pointerOrigin.x);
                currentTransform.y = dragOrigin.y + (event.y - pointerOrigin.y);
                d3.select(this).attr('transform', 'translate(' + currentTransform.x + ',' + currentTransform.y + ')');
            })
            .on('end', function(event) {
                this.style.cursor = 'grab';
                PlasmidMapperEdits.setLegendPosition(qryId, legendKey, currentTransform.x, currentTransform.y);
            }));
    }

    function plotLegend(qryId, data) {
        var orf_labels_for_rect = {
            'ARGs': 'args', 'Insertion sequences': 'isel', 'Transposons': 'transposase',
            'Virulence factors': 'virulence', 'Biocide and metal resistance': 'biocidemetal',
            'Integron': 'integrase', 'Hypothetical proteins': 'hypothetical', 'Other': 'other'
        };
        if (renderSettings.orfLegendStyle === 'rect') {
            var rows = Object.keys(orf_labels_for_rect).map(function(txt) {
                return { swatchColor: renderSettings.colors[orf_labels_for_rect[txt]], label: txt };
            });
            plotLegendRect(qryId, ORF_LEGEND_KEY, rows, 160, 'legendText');
            return;
        }

        var deg = Math.PI / 180,
            pi2 = 2 * Math.PI;
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

        // Total angular width this legend's arc actually needs, computed
        // from the real (fixed) category list rather than assumed --
        // sum of every coef times step, matching how the loop below lays
        // categories out one after another starting at sAngle.
        var totalCoef = 0;
        $.each(orf_labels, function(txt, d) { totalCoef += d.coef; });
        // Widened from 2.7deg/unit: at the smaller default radius (200,
        // was 260) the legend's curved category text was visibly cramped
        // and overlapping ("AR Insertion sTransp..."). A bigger per-unit
        // step spreads the same 8 categories across more of the circle so
        // each curved label has room to read cleanly.
        var legendStepDeg = 6;
        var requiredArcWidth = totalCoef * (legendStepDeg * deg);
        var step = legendStepDeg * deg;
        var radius = controls.radius;

        // Default (for a plasmid ID not in the hand-curated lengendAngle
        // map above): avoid both any user-added zoom band (data.annotations)
        // and plotBlastLegend's own placement, computed dynamically via
        // findClearAngle rather than a fixed guess -- a fixed default used
        // to just avoid the OTHER legend's fixed default, which didn't
        // account for annotation bands at all (confirmed: selecting a band
        // under a legend's fixed position made both unreadable together).
        var coord2Angle = d3.scaleLinear().range([0, pi2]).domain([0, data.qlen]);
        var occupiedRanges = annotationOccupiedRanges(data, coord2Angle);
        var autoAngle = lengendAngle[qryId] ? lengendAngle[qryId] : findClearAngle(occupiedRanges, requiredArcWidth);

        // Draws (or redraws, clearing first) all of this legend's content
        // AT A GIVEN ANGLE -- factored out so a drag can call this
        // directly with the live pointer angle on every tick, instead of
        // building geometry once and rotating the whole group as a
        // transform. That matters specifically because of `revert` below:
        // it decides whether text/swatch radii need to flip so the text
        // reads right-side-up, and that decision depends on which half of
        // the circle sAngle falls in -- a decision baked in once at build
        // time goes stale the moment a transform-only drag moves the
        // group into the opposite hemisphere (confirmed: this is exactly
        // what caused upside-down legend text when dragged to the bottom).
        // Recomputing revert from the CURRENT angle on every redraw keeps
        // the text correctly oriented throughout the drag, not just after
        // it settles.
        function drawAt(sAngle) {
            legend.selectAll('*').remove();

            // Stash this legend's occupied window so plotBlastLegend
            // (called later, from a different handler once BLAST rings
            // are selected) can avoid it too. Only meaningful/updated
            // while this is the auto-placed (non-dragged) angle; see the
            // dragEndAngle handling below.
            controls.orfLegendOccupied = { start: sAngle, end: sAngle + requiredArcWidth };

            var lgTxt = legend.append('text');
            var k = 0, ta, tb, bias, sa, sb;

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
                    .style('fill', renderSettings.colors[d.kl])
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
                    .style('font-size', renderSettings.fonts.legendText.size + 'px')
                    .style('font-weight', 600)
                    .style('font-family', renderSettings.fonts.legendText.family);
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

        var legend = d3.select('#focus').append('g');
        var saved = PlasmidMapperEdits.getLegendPosition(qryId, ORF_LEGEND_KEY);
        var initialAngle = saved && typeof saved.angleDeg === 'number' ? saved.angleDeg * deg : autoAngle;
        drawAt(initialAngle);
        makeLegendDraggable(legend, qryId, ORF_LEGEND_KEY, initialAngle, drawAt);

    }

    function plotBlastLegend(qryId, data) {
        if (renderSettings.blastLegendStyle === 'rect') {
            var rows = selected_alignments.map(function(key, i) {
                var subjectName = key.split('$')[1] || key;
                // Truncate by character count for a predictable box width,
                // consistent with this file's chunkSubstr precedent rather
                // than a getBBox() measurement pass for a row set whose
                // count/content varies with the current ring selection.
                if (subjectName.length > 20) subjectName = subjectName.slice(0, 19) + '…';
                return { swatchColor: Color_collection[i % Color_collection.length], label: subjectName };
            });
            plotLegendRect(qryId, BLAST_LEGEND_KEY, rows, 180, 'blastLegendText');
            return;
        }

        var deg = Math.PI / 180,
            pi2 = 2 * Math.PI;

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
        // Widened slightly from 4.5deg/unit, alongside plotLegend's
        // ORF-category arc -- but NOT by the same ~2.2x factor: this
        // legend's total width already scales with the number of BLAST-
        // subject rings selected (requiredArcWidth below), so applying the
        // same multiplier here blew up into a huge arc for a many-subject
        // example that collided with plotLegend's own arc. A modest bump
        // keeps each plasmid-name label readable without that blowup.
        var blastLegendStepDeg = 5.5;
        var coefLocal = 2, stepLocal = blastLegendStepDeg * deg;
        var step = blastLegendStepDeg * deg;
        var radius = controls.radius;
        var r1 = radius + 18, r2 = radius + 28, coef = 2;

        // This legend's arc width varies with how many BLAST-subject rings
        // are checked (a text arc is shared per pair of rings) -- computed
        // from the real selection count rather than assumed fixed.
        var requiredArcWidth = Math.ceil(selected_alignments.length / 2) * coefLocal * stepLocal;

        // Default: avoid any user-added zoom band AND plotLegend's already-
        // chosen window (stashed on controls.orfLegendOccupied when it ran),
        // computed dynamically via findClearAngle rather than a fixed guess
        // that only ever avoided the other legend's own fixed default.
        var coord2Angle = d3.scaleLinear().range([0, pi2]).domain([0, data.qlen]);
        var occupiedRanges = annotationOccupiedRanges(data, coord2Angle);
        if (controls.orfLegendOccupied) occupiedRanges.push(controls.orfLegendOccupied);
        var autoAngle = lengendAngle[qryId] ? lengendAngle[qryId] : findClearAngle(occupiedRanges, requiredArcWidth);

        // See the matching comment in plotLegend()'s drawAt(): factored out
        // so a drag can rebuild geometry (including the text-orientation
        // "revert" decision) at the live angle on every tick, instead of
        // baking that decision in once and rotating pre-built content --
        // which left text upside-down when dragged into the opposite
        // hemisphere.
        function drawAt(sAngle) {
            legend.selectAll('*').remove();

            var k = 0, ta, tb, bias, sa, sb, r;
            var lgTxt = legend.append('text');

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
                    .style('font-size', renderSettings.fonts.blastLegendText.size + 'px')
                    .style('font-weight', 600)
                    .style('font-family', renderSettings.fonts.blastLegendText.family);
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

        var legend = d3.select('#focus').append('g');
        var saved = PlasmidMapperEdits.getLegendPosition(qryId, BLAST_LEGEND_KEY);
        var initialAngle = saved && typeof saved.angleDeg === 'number' ? saved.angleDeg * deg : autoAngle;
        drawAt(initialAngle);
        makeLegendDraggable(legend, qryId, BLAST_LEGEND_KEY, initialAngle, drawAt);

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
            .style('font-size', renderSettings.fonts.axisLabel.size + 'px')
            .style('font-family', renderSettings.fonts.axisLabel.family)
            .style('font-weight', 'bold');

    }

    function plotBlastRings(data, radius, color_index) {
        var qLen = data.qlen;
        var bl_focus = d3.select('#bl-focus');
        var coord2Angle = d3.scaleLinear().range([0, 2 * Math.PI]).domain([0, qLen])
        var arcW = renderSettings.ringThickness,
            panelW = renderSettings.ringThickness,
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
            .style('font-size', renderSettings.fonts.axisLabel.size + 'px')
            .style('font-weight', 600)
            .style('font-family', renderSettings.fonts.axisLabel.family);

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

        enableRegionSelection(qryfocus, qryId, qLen, coord2Angle, radius);

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
        var sColor = '#dd3497';
        var recR = [orfR[0] - 10, radius + 15];

        // Each zoom band gets its own stacked-outward secondRadius and its
        // own tcoord2Angle scale, stored per-band in bandContexts instead
        // of being reassigned into a single shared variable each loop
        // iteration. With only one band the old shared-variable approach
        // happened to work (there was nothing else to collide with); with
        // two or more, every band's ORFs/labels would render using
        // whichever band's values were left over after the loop finished.
        var sortedAnnotations = (data.annotations || []).slice().sort(function(a, b) { return a.sidx - b.sidx; });
        var bandCount = sortedAnnotations.length;
        var BAND_STEP = bandCount > 0 ? Math.min(24, Math.max(10, 200 / bandCount)) : 24;
        var bandContexts = [];

        $.each(sortedAnnotations, function(i, d) {

            var bandSecondRadius = recR[0] + 30 + i * BAND_STEP;

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
            var maxExpansion = Math.min(2 * qryLen, qLen / 15);
            // Clamp each side's expansion against the midpoint gap to the
            // neighboring band (in sorted order) so one band's "zoomed
            // canvas" never crosses into an adjacent band's own core
            // [sidx, eidx] range -- without this, two bands anywhere near
            // each other would have their tick marks/ORFs/labels overlap
            // in angle even though they're drawn at different radii.
            var prevBand = sortedAnnotations[i - 1];
            var nextBand = sortedAnnotations[i + 1];
            var leftExpansion = prevBand ? Math.max(0, Math.min(maxExpansion, (d.sidx - prevBand.eidx) / 2)) : maxExpansion;
            var rightExpansion = nextBand ? Math.max(0, Math.min(maxExpansion, (nextBand.sidx - d.eidx) / 2)) : maxExpansion;
            var arcSidx = d.sidx - leftExpansion,
                arcEidx = d.eidx + rightExpansion;

            // A light filled background behind the whole expanded band
            // (main-ring boundary out through the band's own tick axis and
            // duplicated ORF arcs) so the zoomed-in region reads as a
            // single visual unit at a glance, rather than only being
            // marked by the thin dashed boundary/connector lines. Drawn
            // first (before the boundary arc, ticks, and ORFs below) so it
            // paints underneath all of them.
            qryfocus.append('path')
                .attr('class', 'band-highlight')
                .attr('d', d3.arc()
                    .innerRadius(recR[1])
                    .outerRadius(bandSecondRadius + 10)
                    .startAngle(coord2Angle(arcSidx))
                    .endAngle(coord2Angle(arcEidx))
                )
                .style('fill', renderSettings.colors && renderSettings.colors.bandHighlight || '#ffe58f')
                .style('opacity', 0.25)
                .style('pointer-events', 'none');

            var ticks = d3.range(d.sidx, d.eidx, 1e3)
            var bandTcoord2Angle = d3.scaleLinear().range([coord2Angle(arcSidx) % pi2, coord2Angle(arcEidx) % pi2]).domain([d.sidx, d.eidx]);
            var tx = d3.scaleBand()
                .range([coord2Angle(arcSidx), coord2Angle(arcEidx)])
                .domain(d3.range(d.sidx, d.eidx));

            var ty = d3.scaleRadial()
                .range([bandSecondRadius - 2, bandSecondRadius]) // Domain will be define later.
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
                .attr("d", getUnaligned_deletion(recR[1], bandSecondRadius, coord2Angle(d.sidx + (qryLen / 2)), coord2Angle(arcSidx), coord2Angle(arcEidx)))
                .style('stroke', sColor)
                .style("stroke-dasharray", ("1,1"))
                .style('stroke-width', 1)
                .attr('fill', 'none');

            // Fall back to a positional id when the annotation has none --
            // hand-curated data predating the Python pipeline's id-stamping
            // (js/ref_data.js) has annotations with no `id` field at all,
            // which would otherwise produce DOM ids like "orf-123-undefined".
            var bandId = d.id != null ? d.id : ('band-idx-' + i);
            bandContexts.push({ ann: { sidx: d.sidx, eidx: d.eidx, id: bandId }, secondRadius: bandSecondRadius, tcoord2Angle: bandTcoord2Angle });

        });

        $.each(data.orfs, function(i, d) {

            // An ORF can fall inside more than one zoom band; render it
            // once per matching band (each band's ring shows what's
            // actually zoomed there -- dropping it from one would look
            // like missing data, and "first match" has no principled
            // meaning since band order is creation order, not spatial).
            var matchedBands = bandContexts.filter(function(bc) {
                return d.eidx >= bc.ann.sidx && d.sidx <= bc.ann.eidx;
            });

            var outerOrfArc = qryfocus.append("path")
                .attr('class', "orf " + d.type)
                .attr("d", getArrowedArc(orfR[0], orfR[1], coord2Angle(d.sidx),
                    coord2Angle(d.eidx), d.strand == 1))
                .style('fill', renderSettings.colors[d.type])
                .style('stroke', '#737373')
                .style('stroke-width', 0.3);
            attachOrfTooltip(outerOrfArc, d);

            $.each(matchedBands, function(bandIdx, bc) {
                // Suffix every id with this band's id so the same ORF
                // rendered in multiple bands doesn't collide in the DOM.
                var idSuffix = d.id + '-' + bc.ann.id;
                var innerOrfArc = qryfocus.append("path")
                    .attr('class', "orf " + d.type)
                    .attr('id', 'orf-' + idSuffix)
                    .attr('data-orf-id', d.id)
                    .attr('data-band-id', bc.ann.id)
                    .attr("d", getArrowedArc(bc.secondRadius + 2, bc.secondRadius + 8, bc.tcoord2Angle(d.sidx),
                        bc.tcoord2Angle(d.eidx), d.strand == 1))
                    .style('fill', renderSettings.colors[d.type])
                    .style('stroke', '#737373')
                    .style('stroke-width', 0.3).on('click', function(event) {
                        var orfId = d3.select(this).attr('data-orf-id');
                        var bandId = d3.select(this).attr('data-band-id');
                        var suffix = orfId + '-' + bandId;
                        var curStat = d3.select('#line-' + suffix).attr('display')
                        d3.select('#line-' + suffix).attr('display', curStat == 'none' ? 'block' : 'none');
                        d3.select('#txt-' + suffix).attr('display', curStat == 'none' ? 'block' : 'none');
                    });
                attachOrfTooltip(innerOrfArc, d);
            });

            if (d.type == 'hypothetical') return;

            $.each(matchedBands, function(bandIdx, bc) {
                var idSuffix = d.id + '-' + bc.ann.id;
                var secondRadius = bc.secondRadius,
                    tcoord2Angle = bc.tcoord2Angle;

                textg.append('path')
                    .attr('id', 'line-' + idSuffix)
                    // Visible by default: a zoom band is the region the
                    // user is actively curating for a figure, so its ORF
                    // labels should be visible immediately rather than
                    // requiring a click to reveal. Still click-to-hide
                    // via the ORF-arrow handler above if it gets crowded.
                    .attr("d", getORFLables(secondRadius + 8, secondRadius + 18,
                        tcoord2Angle(d.sidx), tcoord2Angle(d.eidx)))
                    .style('stroke', '#000')
                    .style("stroke-dasharray", ("1,1"))
                    .style('stroke-width', '0.1')
                    .style('fill', 'none');


                // Must match getORFLables' own midpoint formula (below) so
                // the label lands where its leader line actually points --
                // this used to omit the /2, which for eidx > sidx reduces
                // to just `eidx` (the ORF's end coordinate, not its middle),
                // scattering labels away from their ORFs for anything but
                // the narrowest features.
                var midPoint = (d.sidx + d.eidx) / 2,
                    x = (secondRadius + 18) * Math.cos(tcoord2Angle(midPoint) - half_pi),
                    y = (secondRadius + 18) * Math.sin(tcoord2Angle(midPoint) - half_pi);
                var labelOverride = d._labelOverride;
                var initRotation = 0;
                if (labelOverride && typeof labelOverride.x === 'number') {
                    x = labelOverride.x;
                    y = labelOverride.y;
                    initRotation = labelOverride.rotation || 0;
                }
                var labelText = textg.append('g')
                    .append('text')
                    .attr('id', 'txt-' + idSuffix)
                    .attr('x', x)
                    .attr('y', y)
                    .attr('transform', 'rotate(' + initRotation + ',' + x + ',' + y + ')')
                    .style("font-size", renderSettings.fonts.orfLabel.size + "px")
                    .style('font-weight', 600)
                    .style('font-family', renderSettings.fonts.orfLabel.family)
                    .style('cursor', 'grab')
                    .text((labelOverride && labelOverride.text) || d.dscr.replace('family transposase', ''))
                    .on("dblclick", function(event) {
                        event.preventDefault();
                        openLabelTextEditor(this, qryId, d.id, $(this).text());
                    })
                    .on('mousewheel', function(event) {

                        event.preventDefault();
                        var sig = event.wheelDelta > 0 ? 1 : -1;
                        var x1 = parseInt($(this).attr('x')),
                            y1 = parseInt($(this).attr('y'));

                        var line = qryfocus.select('#line-' + idSuffix);
                        var sp = line.attr("d").split(" ");

                        var tr = d3.select(this).attr("transform");
                        pp = tr.replace('rotate(', '').replace(');', '').split(',')
                        var newRotation = parseInt(pp[0]) + sig * 5;
                        $(this).attr('transform', 'rotate(' + newRotation +
                            ',' + x1 + ',' + y1 + ')');
                        PlasmidMapperEdits.setLabelPosition(qryId, 'orf-' + d.id, x1, y1, newRotation);
                    });

                // d3.drag() (not hand-rolled mousedown/mousemove/mouseup)
                // captures the pointer for the whole drag gesture
                // regardless of which element is under the cursor -- the
                // old per-element mousemove handler lost events whenever
                // a fast movement exited the label's small hit-area,
                // which is what produced the jerky/laggy dragging.
                // .container(qryfocus.node()) matches the reference-
                // element convention already used elsewhere in this
                // codebase (js/pl_editor.js's enableRegionSelection/
                // plotFreeformLabels both call d3.pointer(event,
                // qryfocus.node())), so drag coordinates land in the same
                // space the leader-line path was already drawn in.
                labelText.call(d3.drag()
                    .container(qryfocus.node())
                    .on('start', function(event) {
                        this.style.cursor = "grabbing";
                        d3.select(this).style('font-size', '10px');
                    })
                    .on('drag', function(event) {
                        var x1 = event.x, y1 = event.y;

                        var line = qryfocus.select('#line-' + idSuffix);
                        var sp = line.attr("d").split(" ");
                        var x, y, x0 = sp[1],
                            y0 = sp[2];
                        var bias = 0;
                        x = x1, y = y1;
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
                    })
                    .on('end', function(event) {
                        this.style.cursor = "grab";
                        d3.select(this).style('font-size', renderSettings.fonts.orfLabel.size + 'px');

                        var finalX = parseFloat($(this).attr('x')),
                            finalY = parseFloat($(this).attr('y'));
                        var tr = d3.select(this).attr("transform");
                        var rotation = parseFloat(tr.replace('rotate(', '').split(',')[0]) || 0;
                        PlasmidMapperEdits.setLabelPosition(qryId, 'orf-' + d.id, finalX, finalY, rotation);
                    }));
            });

        });

        plotFreeformLabels(qryfocus, textg, qryId, data._freeformLabels || []);
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

    function rerenderCurrentPlasmid() {
        d3.select("#main-svg").selectAll('*').remove();
        d3.select("#tbl-main").selectAll('*').remove();
        clear_canvas();
        selected_alignments = [];
        controls = update_page(qryId);
    }
    // Exposed globally so pl_editor.js can trigger a full re-render after an
    // annotation/label edit, reusing the exact same clear-and-redraw path as
    // the query-select dropdown.
    window.rerenderCurrentPlasmid = rerenderCurrentPlasmid;

    $("#qryselect").on('change', function() {
        qryId = this.value;
        rerenderCurrentPlasmid();
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

            var effectiveDataForLegend = PlasmidMapperEdits.mergeEdits(qryId, Contig_ref[qryId]);
            plotBlastLegend(qryId, effectiveDataForLegend);

        }

        // Labels/leader-lines (inside .qry-focus, appended by plotPlasmid()
        // before this handler creates #bl-focus/the BLAST legend) must stay
        // on top so they remain clickable/draggable -- otherwise later-drawn
        // BLAST ring/legend paths silently absorb pointer events meant for
        // double-click/drag on a label underneath them in SVG paint order.
        d3.select('#focus').select('.qry-focus').raise();

    });


    $('#genBtn').on('click', function(event) {
        var format = $('#export-format-select').val();
        exportFigure(qryId, format);
    });


    $('.inpt').on('change', function() {
        var id = $(this).attr('id');
        if (id == 'width-inpt') {
            $('#canvas').attr('width', $(this).val() + 'px')

        } else {
            $('#canvas').attr('height', $(this).val() + 'px')

        }
    });

    $('#addLabelBtn').on('click', function(event) {
        armFreeformLabelPlacement(qryId);
    });

    $('#exportEditsBtn').on('click', function(event) {
        PlasmidMapperEdits.exportEdits(qryId);
    });

    $('#importEditsInput').on('change', function(event) {
        var file = event.target.files[0];
        if (!file) return;
        PlasmidMapperEdits.importEditsFromFile(file, function(err) {
            if (err) {
                alert('Could not import edits file: ' + err.message);
                return;
            }
            rerenderCurrentPlasmid();
        });
        event.target.value = '';
    });

    $('#radius-inpt').on('change', function() {
        var val = parseFloat($(this).val());
        if (!isNaN(val) && val > 0) {
            renderSettings.radius = val;
            rerenderCurrentPlasmid();
        }
    });

    $('#ring-thickness-inpt').on('change', function() {
        var val = parseFloat($(this).val());
        if (!isNaN(val) && val > 0) {
            renderSettings.ringThickness = val;
            rerenderCurrentPlasmid();
        }
    });

    $('#ring-spacing-inpt').on('change', function() {
        var val = parseFloat($(this).val());
        if (!isNaN(val) && val < 0) {
            renderSettings.radiusStep = val;
            rerenderCurrentPlasmid();
        }
    });

    $('#orf-legend-style-select').on('change', function() {
        renderSettings.orfLegendStyle = $(this).val();
        rerenderCurrentPlasmid();
    });

    $('#blast-legend-style-select').on('change', function() {
        renderSettings.blastLegendStyle = $(this).val();
        rerenderCurrentPlasmid();
    });

    var qrySelectEl = document.getElementById('qryselect');
    var defaultQryId = qrySelectEl.options.length > 1 ? qrySelectEl.options[1].value : qrySelectEl.options[0].value;
    $('#qryselect').val(defaultQryId).change();
});