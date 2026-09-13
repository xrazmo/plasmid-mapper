// Interactive figure-editing layer for mapper.html: lets a user select a
// genomic region to render as a zoomed "band" (generalizing the
// previously-hardcoded Contig_ref[...].annotations), and add/edit/remove
// labels — all stored in a separate "edits overlay" kept out of
// ref_data.js/pl_data.js so re-running the Python pipeline never discards
// a user's manual figure polish.
//
// Loaded after pl_mapper.js. Exposes a single global, PlasmidMapperEdits,
// used by pl_mapper.js at the two integration points documented inline
// there (mergeEdits() call before plotPlasmid(), and persistence hooks on
// the existing label drag handlers).
var PlasmidMapperEdits = (function() {
    var STORAGE_KEY = 'plasmid-mapper-edits';
    var store = loadFromStorage();

    function loadFromStorage() {
        try {
            var raw = window.localStorage.getItem(STORAGE_KEY);
            return raw ? JSON.parse(raw) : {};
        } catch (e) {
            // Private browsing / blocked storage: degrade to in-memory only,
            // edits still work within the session, just don't persist.
            return {};
        }
    }

    function persist() {
        try {
            window.localStorage.setItem(STORAGE_KEY, JSON.stringify(store));
        } catch (e) {
            // Ignore quota/availability errors; edits remain usable in-memory.
        }
    }

    function ensurePlasmidEntry(qryId) {
        if (!store[qryId]) {
            store[qryId] = {
                annotations: { added: [], removed: [] },
                labels: {}
            };
        }
        return store[qryId];
    }

    function uid(prefix) {
        return prefix + '-' + Math.random().toString(36).slice(2, 10);
    }

    // --- Annotations (zoomed "band" regions) ---------------------------

    function addAnnotation(qryId, sidx, eidx) {
        var entry = ensurePlasmidEntry(qryId);
        var id = uid('ann');
        entry.annotations.added.push({ id: id, sidx: sidx, eidx: eidx });
        persist();
        return id;
    }

    function removeAnnotation(qryId, annotationId) {
        var entry = ensurePlasmidEntry(qryId);
        var addedIdx = entry.annotations.added.findIndex(function(a) { return a.id === annotationId; });
        if (addedIdx > -1) {
            // Session-added annotation never touched the base data; just drop it.
            entry.annotations.added.splice(addedIdx, 1);
        } else {
            entry.annotations.removed.push(annotationId);
        }
        persist();
    }

    // --- Labels ----------------------------------------------------------

    function setLabelPosition(qryId, labelKey, x, y, rotation) {
        var entry = ensurePlasmidEntry(qryId);
        var existing = entry.labels[labelKey] || {};
        entry.labels[labelKey] = Object.assign({}, existing, { x: x, y: y, rotation: rotation });
        persist();
    }

    function setLabelText(qryId, labelKey, text) {
        var entry = ensurePlasmidEntry(qryId);
        var existing = entry.labels[labelKey] || {};
        entry.labels[labelKey] = Object.assign({}, existing, { text: text });
        persist();
    }

    function addFreeformLabel(qryId, x, y, text, pointsTo) {
        var entry = ensurePlasmidEntry(qryId);
        var key = uid('freeform');
        entry.labels[key] = { x: x, y: y, rotation: 0, text: text, pointsTo: pointsTo, isFreeform: true };
        persist();
        return key;
    }

    function deleteLabel(qryId, labelKey) {
        var entry = ensurePlasmidEntry(qryId);
        if (entry.labels[labelKey] && entry.labels[labelKey].isFreeform) {
            delete entry.labels[labelKey];
        } else {
            // ORF-backed label: revert to default instead of leaving a hole,
            // since the underlying ORF arrow stays regardless of the override.
            delete entry.labels[labelKey];
        }
        persist();
    }

    // --- Merge into base Contig_ref data at render time -------------------

    function mergeEdits(qryId, contigRefEntry) {
        if (!contigRefEntry) return contigRefEntry;
        var entry = store[qryId];
        if (!entry) return contigRefEntry;

        var danglingWarnings = [];

        var baseAnnotations = (contigRefEntry.annotations || []);
        var removedIds = {};
        (entry.annotations.removed || []).forEach(function(id) { removedIds[id] = true; });

        var validRemovedIds = baseAnnotations.map(function(a) { return a.id; });
        (entry.annotations.removed || []).forEach(function(id) {
            if (validRemovedIds.indexOf(id) === -1) {
                danglingWarnings.push('a removed band (' + id + ') no longer exists in the regenerated data');
            }
        });

        var effectiveAnnotations = baseAnnotations
            .filter(function(a) { return !removedIds[a.id]; })
            .concat(entry.annotations.added || []);

        var effectiveOrfs = (contigRefEntry.orfs || []).map(function(orf) {
            var override = entry.labels['orf-' + orf.id];
            if (!override) return orf;
            var merged = Object.assign({}, orf);
            if (override.text !== undefined) merged.dscr = override.text;
            merged._labelOverride = override;
            return merged;
        });

        var validOrfIds = {};
        (contigRefEntry.orfs || []).forEach(function(o) { validOrfIds[o.id] = true; });
        Object.keys(entry.labels).forEach(function(key) {
            if (key.indexOf('orf-') === 0) {
                var orfId = key.slice(4);
                var matches = (contigRefEntry.orfs || []).some(function(o) { return String(o.id) === orfId; });
                if (!matches) {
                    danglingWarnings.push('a saved label edit (' + key + ') no longer matches an ORF in the regenerated data');
                }
            }
        });

        if (danglingWarnings.length > 0) {
            showDanglingEditsNotice(danglingWarnings);
        }

        var freeformLabels = Object.keys(entry.labels)
            .filter(function(key) { return key.indexOf('freeform-') === 0; })
            .map(function(key) {
                return Object.assign({ id: key }, entry.labels[key]);
            });

        var merged = Object.assign({}, contigRefEntry);
        merged.annotations = effectiveAnnotations;
        merged.orfs = effectiveOrfs;
        merged._freeformLabels = freeformLabels;
        return merged;
    }

    function showDanglingEditsNotice(messages) {
        if (document.getElementById('edits-dangling-notice')) return;
        var el = document.createElement('div');
        el.id = 'edits-dangling-notice';
        el.className = 'alert alert-warning alert-dismissible';
        el.style.margin = '8px 0';
        el.innerHTML = '<button type="button" class="close" data-dismiss="alert" aria-label="Close">' +
            '<span aria-hidden="true">&times;</span></button>' +
            'Some saved figure edits could not be reapplied: ' + messages.join('; ') + '.';
        var container = document.querySelector('#main-section') || document.body;
        container.parentNode.insertBefore(el, container);
        el.querySelector('.close').addEventListener('click', function() { el.remove(); });
    }

    function exportEdits(qryId) {
        var payload = qryId ? { qryId: qryId, edits: store[qryId] || {} } : { all: store };
        var blob = new Blob([JSON.stringify(payload, null, 2)], { type: 'application/json' });
        var name = (qryId || 'all') + '.edits.json';
        saveAs(blob, name);
    }

    function importEditsFromFile(file, callback) {
        var reader = new FileReader();
        reader.onload = function(e) {
            try {
                var payload = JSON.parse(e.target.result);
                if (payload.qryId) {
                    store[payload.qryId] = payload.edits;
                } else if (payload.all) {
                    store = payload.all;
                }
                persist();
                if (callback) callback(null);
            } catch (err) {
                if (callback) callback(err);
            }
        };
        reader.readAsText(file);
    }

    return {
        mergeEdits: mergeEdits,
        addAnnotation: addAnnotation,
        removeAnnotation: removeAnnotation,
        setLabelPosition: setLabelPosition,
        setLabelText: setLabelText,
        addFreeformLabel: addFreeformLabel,
        deleteLabel: deleteLabel,
        exportEdits: exportEdits,
        importEditsFromFile: importEditsFromFile,
        getEntry: ensurePlasmidEntry
    };
})();

// Inline text editing for a label: double-click opens a foreignObject-
// hosted <input> positioned over the label (SVG has no native editable
// text primitive), pre-filled with the current text. Enter/blur commits
// the change via PlasmidMapperEdits and re-renders; Escape cancels.
//
// labelKey is the full edits-overlay key (e.g. "orf-63908" or a
// freeform-<uuid> key) for freeform labels, or just the numeric ORF id
// for ORF-backed labels (setLabelText below normalizes this).
function openLabelTextEditor(textEl, qryId, orfOrLabelId, currentText, onFreeform) {
    var svg = document.getElementById('main-svg');
    var textNode = d3.select(textEl);
    var bbox = textEl.getBBox();
    var ctm = textEl.getScreenCTM();
    var svgRect = svg.getBoundingClientRect();

    var fo = document.createElementNS('http://www.w3.org/2000/svg', 'foreignObject');
    var padding = 4;
    fo.setAttribute('x', bbox.x - padding);
    fo.setAttribute('y', bbox.y - padding);
    fo.setAttribute('width', Math.max(120, bbox.width + padding * 2));
    fo.setAttribute('height', bbox.height + padding * 2 + 4);

    var input = document.createElement('input');
    input.type = 'text';
    input.value = currentText;
    input.style.width = '100%';
    input.style.fontSize = '10px';
    input.style.border = '1px solid #2c7fb8';
    input.style.boxSizing = 'border-box';

    fo.appendChild(input);
    textEl.parentNode.appendChild(fo);
    input.focus();
    input.select();

    var settled = false;

    function commit() {
        if (settled) return;
        settled = true;
        var newText = input.value;
        fo.remove();
        if (newText === currentText) return;
        var key = onFreeform ? orfOrLabelId : ('orf-' + orfOrLabelId);
        PlasmidMapperEdits.setLabelText(qryId, key, newText);
        textNode.text(newText);
    }

    function cancel() {
        if (settled) return;
        settled = true;
        fo.remove();
    }

    input.addEventListener('keydown', function(e) {
        if (e.key === 'Enter') { e.preventDefault(); commit(); }
        else if (e.key === 'Escape') { e.preventDefault(); cancel(); }
    });
    input.addEventListener('blur', commit);
}

// Region ("band") selection: click-drag directly on the query ring to pick
// a start/end coordinate and create a new zoomed-region annotation,
// generalizing what used to be a hardcoded Contig_ref[...].annotations
// entry. Only armed while "Add band" mode is toggled on (see
// renderAnnotationControls) so normal label-dragging elsewhere on the
// figure isn't accidentally hijacked.
var regionSelectionArmed = false;

function enableRegionSelection(qryfocus, qryId, qLen, coord2Angle, radius) {
    var half_pi = Math.PI / 2;
    // The visible ring path is only 0.1 radius units thick — far too thin
    // to reliably grab with a pointer. A separate, wider, transparent arc
    // layered on top is the actual drag target.
    var hitArc = qryfocus.append('path')
        .attr('class', 'region-select-hit')
        .attr('d', d3.arc()
            .innerRadius(radius - 10)
            .outerRadius(radius + 10)
            .startAngle(0)
            .endAngle(2 * Math.PI))
        .style('fill', 'transparent')
        .style('cursor', 'crosshair')
        .style('pointer-events', function() { return regionSelectionArmed ? 'all' : 'none'; });

    var dragStartCoord = null;
    var previewPath = null;

    function angleFromPointer(event) {
        var p = d3.pointer(event, qryfocus.node());
        var angle = Math.atan2(p[1], p[0]) + half_pi;
        // Snap floating-point noise near the 0/2*PI seam (the plasmid's 0kb
        // origin tick) to exactly 0 before wrapping negatives into [0, 2*PI).
        // Without this, a pointer angle computed as e.g. -5.68e-14 (meant to
        // be exactly 0) wraps all the way around to ~2*PI, i.e. a coordinate
        // near qLen instead of near 0 -- exactly the seam a real drag
        // starting at the plasmid's origin tick would hit.
        var epsilon = 1e-9;
        if (Math.abs(angle) < epsilon || Math.abs(angle - 2 * Math.PI) < epsilon) {
            angle = 0;
        } else if (angle < 0) {
            angle += 2 * Math.PI;
        }
        return angle;
    }

    hitArc.on('mousedown', function(event) {
        if (!regionSelectionArmed) return;
        event.preventDefault();
        var angle = angleFromPointer(event);
        dragStartCoord = Math.round(coord2Angle.invert(angle));
        previewPath = qryfocus.append('path')
            .attr('class', 'region-select-preview')
            .style('fill', '#dd3497')
            .style('opacity', 0.3);
    });

    hitArc.on('mousemove', function(event) {
        if (!regionSelectionArmed || dragStartCoord === null) return;
        var angle = angleFromPointer(event);
        var currentCoord = Math.round(coord2Angle.invert(angle));
        var sidx = Math.min(dragStartCoord, currentCoord);
        var eidx = Math.max(dragStartCoord, currentCoord);
        previewPath.attr('d', d3.arc()
            .innerRadius(radius - 15)
            .outerRadius(radius + 15)
            .startAngle(coord2Angle(sidx))
            .endAngle(coord2Angle(eidx)));
    });

    function finishDrag(event) {
        if (!regionSelectionArmed || dragStartCoord === null) return;
        var angle = angleFromPointer(event);
        var currentCoord = Math.round(coord2Angle.invert(angle));
        var sidx = Math.min(dragStartCoord, currentCoord);
        var eidx = Math.max(dragStartCoord, currentCoord);
        dragStartCoord = null;
        if (previewPath) { previewPath.remove(); previewPath = null; }
        if (eidx - sidx < 10) return; // too small to be an intentional selection
        PlasmidMapperEdits.addAnnotation(qryId, sidx, eidx);
        window.rerenderCurrentPlasmid();
    }

    hitArc.on('mouseup', finishDrag);
    hitArc.on('mouseleave', function(event) {
        // Don't cancel on leave; the user may drag slightly outside the
        // hit-arc's radial bounds while still intending to complete the
        // selection. Only mouseup (wherever it happens) finalizes it.
    });
    d3.select(window).on('mouseup.region-select-' + qryId, function(event) {
        if (dragStartCoord !== null) finishDrag(event);
    });
}

function setRegionSelectionArmed(armed) {
    regionSelectionArmed = armed;
    d3.selectAll('.region-select-hit').style('pointer-events', armed ? 'all' : 'none');
}

// Bootstrap list-group UI (outside the SVG) listing the current plasmid's
// zoomed-region bands with a delete button each, plus an "Add band" toggle
// that arms enableRegionSelection's drag gesture on the query ring.
function renderAnnotationControls(qryId, effectiveAnnotations) {
    var containerId = 'annotation-controls';
    var existing = document.getElementById(containerId);
    if (existing) existing.remove();

    var container = document.createElement('div');
    container.id = containerId;
    container.className = 'mt-2 mb-2';

    // Sectioned layout (mapper.html, Step 5 reorg) provides a dedicated
    // mount point inside the "Zoom bands" card; older/simpler layouts
    // (e.g. examples/KPC33_p1/kpc33_viewer.html, not yet reorganized)
    // fall back to inserting the controls directly above #main-section.
    var mount = document.getElementById('annotation-controls-mount');

    var toggleBtn = document.createElement('button');
    toggleBtn.type = 'button';
    toggleBtn.className = 'btn btn-sm btn-outline-primary mr-2';
    toggleBtn.textContent = 'Add band';
    toggleBtn.addEventListener('click', function() {
        var nowArmed = !regionSelectionArmed;
        setRegionSelectionArmed(nowArmed);
        toggleBtn.textContent = nowArmed ? 'Click-drag on the ring, or click here to cancel' : 'Add band';
        toggleBtn.className = nowArmed
            ? 'btn btn-sm btn-primary mr-2'
            : 'btn btn-sm btn-outline-primary mr-2';
    });
    container.appendChild(toggleBtn);

    var list = document.createElement('div');
    list.className = 'list-group list-group-horizontal flex-wrap d-inline-flex';
    (effectiveAnnotations || []).forEach(function(ann) {
        var item = document.createElement('span');
        item.className = 'badge badge-light border mr-1 mb-1 p-2';
        item.textContent = ann.sidx + '–' + ann.eidx + ' bp ';

        var del = document.createElement('button');
        del.type = 'button';
        del.className = 'btn btn-sm btn-link text-danger p-0 ml-1';
        del.textContent = '×';
        del.title = 'Remove this band';
        del.addEventListener('click', function() {
            PlasmidMapperEdits.removeAnnotation(qryId, ann.id);
            window.rerenderCurrentPlasmid();
        });
        item.appendChild(del);
        list.appendChild(item);
    });
    container.appendChild(list);

    if (mount) {
        mount.appendChild(container);
    } else {
        var mainSection = document.getElementById('main-section');
        mainSection.parentNode.insertBefore(container, mainSection);
    }
}

// Renders freeform labels (added via "Add label", not tied to any ORF) as
// text plus a straight leader line from the clicked coordinate to the
// label's (possibly dragged-away) position. getORFLables/getArrowedArc in
// pl_mapper.js are not reusable here: they're declared inside that file's
// $(document).ready(...) closure and never exposed on window, and in any
// case expect a single (innerRadius, outerRadius, angle) radial geometry,
// not two arbitrary (x,y) points — a plain two-point line is both correct
// and simpler for this case. Each freeform label supports the same
// double-click-to-edit-text and drag interactions as ORF-backed labels,
// plus a delete ("×", removes it outright rather than reverting to a
// default, since there is no underlying ORF arrow).
function plotFreeformLabels(qryfocus, textg, qryId, freeformLabels) {
    freeformLabels.forEach(function(label) {
        if (label.pointsTo) {
            textg.append('line')
                .attr('class', 'freeform-leader')
                .attr('x1', label.pointsTo.x)
                .attr('y1', label.pointsTo.y)
                .attr('x2', label.x)
                .attr('y2', label.y)
                .style('stroke', '#000')
                .style('stroke-dasharray', '1,1')
                .style('stroke-width', '0.3');
        }

        var textEl = textg.append('g')
            .append('text')
            .attr('id', 'lbl-' + label.id.replace(/[^a-zA-Z0-9_-]/g, ''))
            .attr('x', label.x)
            .attr('y', label.y)
            .attr('transform', 'rotate(' + (label.rotation || 0) + ',' + label.x + ',' + label.y + ')')
            .style('font-size', '8px')
            .style('font-weight', 600)
            .style('font-family', 'Helvetica')
            .style('cursor', 'grab')
            .text(label.text)
            .on('dblclick', function(event) {
                event.preventDefault();
                openLabelTextEditor(this, qryId, label.id, $(this).text(), true);
            });

        var dragging = false;
        textEl.on('mousedown', function(event) {
            event.preventDefault();
            dragging = true;
            d3.select(this).style('font-size', '10px');
        }).on('mouseleave mouseup', function(event) {
            if (!dragging) return;
            dragging = false;
            d3.select(this).style('font-size', '8px');
            var finalX = parseFloat($(this).attr('x')),
                finalY = parseFloat($(this).attr('y'));
            PlasmidMapperEdits.setLabelPosition(qryId, label.id, finalX, finalY, label.rotation || 0);
        }).on('mousemove', function(event) {
            if (!dragging) return;
            var t = d3.pointer(event, qryfocus.node());
            $(this).attr('x', t[0]).attr('y', t[1]);
            d3.select(this).attr('transform', 'rotate(' + (label.rotation || 0) + ',' + t[0] + ',' + t[1] + ')');
        });
    });
}

// "Add label" click-to-place mode: armed by the toolbar button in
// mapper.html, disarmed after one placement (one-shot, re-click the
// button to add another). On click inside the figure, records the
// SVG-local coordinate, prompts for the label text inline via the same
// foreignObject editor used for existing labels, then creates the
// freeform label and re-renders.
var freeformPlacementArmed = false;

function armFreeformLabelPlacement(qryId) {
    freeformPlacementArmed = true;
    var svg = d3.select('#main-svg');
    svg.style('cursor', 'crosshair');

    function placementHandler(event) {
        if (!freeformPlacementArmed) return;
        freeformPlacementArmed = false;
        svg.style('cursor', null);
        svg.on('click.freeform-placement', null);

        var focus = document.getElementById('focus');
        var pt = d3.pointer(event, focus);
        var key = PlasmidMapperEdits.addFreeformLabel(qryId, pt[0], pt[1], 'New label');
        window.rerenderCurrentPlasmid();
    }

    svg.on('click.freeform-placement', placementHandler);
}

// Hover tooltip shown on any ORF arrow (every ORF everywhere on the
// plasmid, not just ones inside a zoomed band) with the full detail a
// permanent label doesn't show by default: name, identity%, coverage%,
// source database, and category. One reusable <div id="orf-tooltip">
// (declared once in mapper.html/kpc33_viewer.html) is filled and
// positioned on mouseenter, hidden on mouseleave -- cheaper than a
// per-element tooltip and there's only ever one visible at a time anyway.
var ORF_CATEGORY_LABELS = {
    args: 'ARG', isel: 'Insertion sequence', transposase: 'Transposon',
    integrase: 'Integron', virulence: 'Virulence factor',
    biocidemetal: 'Biocide/metal resistance', hypothetical: 'Hypothetical protein',
    other: 'Other', unknown: 'Unknown'
};

function attachOrfTooltip(selection, d) {
    selection.on('mouseenter', function(event) {
        var tooltip = document.getElementById('orf-tooltip');
        if (!tooltip) return;
        var category = ORF_CATEGORY_LABELS[d.type] || d.type;
        var idty = typeof d.idty === 'number' ? d.idty.toFixed(1) + '%' : '—';
        var cov = typeof d.cov === 'number' ? d.cov.toFixed(1) + '%' : '—';
        var dbname = d.dbname || '—';
        tooltip.innerHTML =
            '<span class="orf-tooltip-name">' + (d.dscr || 'Unknown') + '</span>' +
            'Category: ' + category + '<br>' +
            'Identity: ' + idty + ' &nbsp; Coverage: ' + cov + '<br>' +
            'Database: ' + dbname;
        tooltip.style.display = 'block';
        tooltip.style.left = (event.clientX + 12) + 'px';
        tooltip.style.top = (event.clientY + 12) + 'px';
    }).on('mousemove', function(event) {
        var tooltip = document.getElementById('orf-tooltip');
        if (!tooltip) return;
        tooltip.style.left = (event.clientX + 12) + 'px';
        tooltip.style.top = (event.clientY + 12) + 'px';
    }).on('mouseleave', function() {
        var tooltip = document.getElementById('orf-tooltip');
        if (tooltip) tooltip.style.display = 'none';
    });
}

// Exports the current figure in the requested format. Replaces the old
// #genBtn handler, which claimed "Save png" but actually called
// canvas.toDataURL("image/tiff") -- not a real canvas MIME type, so it
// silently produced PNG bytes saved with a .tiff extension.
//
// SVG is the recommended/best-quality option: the already-serialized SVG
// string is saved directly, no rasterization involved. PNG/JPEG go
// through the existing canvas-based rasterization path. PDF is
// deliberately not offered (would need a new external dependency this
// project doesn't otherwise have); the SVG export can be opened in any
// vector editor or the browser's own print-to-PDF when a PDF is needed.
function exportFigure(qryId, format) {
    var svgString = new XMLSerializer().serializeToString(document.querySelector('svg'));

    if (format === 'svg') {
        var blob = new Blob([svgString], { type: 'image/svg+xml;charset=utf-8' });
        saveAs(blob, qryId + '.svg');
        return;
    }

    var canvas = document.getElementById('canvas');
    var ctx = canvas.getContext('2d');
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    var DOMURL = self.URL || self.webkitURL || self;
    var img = new Image();
    var svgBlob = new Blob([svgString], { type: 'image/svg+xml;charset=utf-8' });
    var url = DOMURL.createObjectURL(svgBlob);
    img.onload = function() {
        ctx.drawImage(img, 0, 0);
        if (format === 'jpeg') {
            var jpegDataUrl = canvas.toDataURL('image/jpeg', 0.92);
            saveAs(jpegDataUrl, qryId + '.jpg');
        } else {
            var pngDataUrl = canvas.toDataURL('image/png');
            saveAs(pngDataUrl, qryId + '.png');
        }
        DOMURL.revokeObjectURL(url);
    };
    img.src = url;
}
