// Builds and wires the "Style" modal (font-family/size per text element,
// color per ORF category + the zoom-band highlight fill). Kept in its own
// file rather than growing js/pl_mapper.js further.
//
// Reads/writes window.PlasmidMapperSettings (the same renderSettings
// object pl_mapper.js's Settings sidebar controls already write into) and
// calls window.rerenderCurrentPlasmid() on Apply -- identical mechanism to
// how e.g. #radius-inpt's change handler works, just batched behind one
// button instead of firing on every keystroke (a modal with this many
// fields re-rendering on every keystroke would be janky).
$(function() {

    var FONT_CHOICES = ['Helvetica', 'Arial', 'Tahoma', 'Times New Roman', 'monospace', 'sans-serif'];

    var FONT_FIELD_LABELS = {
        axisLabel: 'Axis labels (ring ticks)',
        orfLabel: 'ORF labels',
        legendText: 'ORF legend text',
        blastLegendText: 'BLAST legend text'
    };

    var COLOR_FIELD_LABELS = {
        args: 'ARGs',
        isel: 'Insertion sequences',
        transposase: 'Transposons',
        integrase: 'Integron',
        virulence: 'Virulence factors',
        biocidemetal: 'Biocide/metal resistance',
        other: 'Other',
        hypothetical: 'Hypothetical proteins',
        bandHighlight: 'Zoom band highlight'
    };

    function buildModal() {
        var settings = window.PlasmidMapperSettings;
        if (!settings) return; // pl_mapper.js hasn't initialized yet

        var $fonts = $('#style-modal-fonts').empty();
        $.each(settings.fonts, function(key, val) {
            var label = FONT_FIELD_LABELS[key] || key;
            var $row = $('<div class="style-font-row"></div>');
            $row.append('<label for="style-font-family-' + key + '">' + label + '</label>');

            var $select = $('<select class="custom-select custom-select-sm" id="style-font-family-' + key + '"></select>');
            $.each(FONT_CHOICES, function(i, fam) {
                var $opt = $('<option></option>').attr('value', fam).text(fam);
                if (fam.toLowerCase() === String(val.family).toLowerCase()) $opt.attr('selected', 'selected');
                $select.append($opt);
            });
            $row.append($select);

            var $size = $('<input type="number" class="form-control form-control-sm" min="3" max="40" step="1">')
                .attr('id', 'style-font-size-' + key)
                .val(val.size);
            $row.append($size);

            $fonts.append($row);
        });

        var $colors = $('#style-modal-colors').empty();
        $.each(settings.colors, function(key, val) {
            var label = COLOR_FIELD_LABELS[key] || key;
            var $swatch = $('<div class="style-color-swatch"></div>');
            $swatch.append('<label for="style-color-' + key + '" class="mb-0 small">' + label + '</label>');
            $swatch.append(
                $('<input type="color" id="style-color-' + key + '">').val(val)
            );
            $colors.append($swatch);
        });
    }

    $('#openStyleModalBtn').on('click', function() {
        buildModal();
        $('#styleModal').modal('show');
    });

    $('#applyStyleModalBtn').on('click', function() {
        var settings = window.PlasmidMapperSettings;
        if (!settings) return;

        $.each(settings.fonts, function(key) {
            settings.fonts[key] = {
                family: $('#style-font-family-' + key).val(),
                size: parseFloat($('#style-font-size-' + key).val()) || settings.fonts[key].size
            };
        });

        $.each(settings.colors, function(key) {
            settings.colors[key] = $('#style-color-' + key).val();
        });

        $('#styleModal').modal('hide');
        if (typeof window.rerenderCurrentPlasmid === 'function') {
            window.rerenderCurrentPlasmid();
        }
    });

});
