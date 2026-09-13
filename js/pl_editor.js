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
