/* MultiMapR GUI -- light client-side helpers (no framework).
   - <input type="color"> Shiny binding
   - live language switch (texts come from the R dictionary)
   - light / dark theme switch
   - character list: client-side search + click-to-toggle
   - click-to-sort statistics table
*/
(function () {
  'use strict';

  var MM = window.MultiMapR = window.MultiMapR || {};
  MM.dict = null;
  MM.lang = document.documentElement.getAttribute('lang') || 'en';

  /* ---- Color input binding ------------------------------------------------ */
  var colorBinding = new Shiny.InputBinding();
  $.extend(colorBinding, {
    find: function (scope) { return $(scope).find('input.mm-color'); },
    getValue: function (el) { return el.value; },
    setValue: function (el, v) { el.value = v; },
    subscribe: function (el, cb) {
      $(el).on('input.mmColor change.mmColor', function () { MM.syncHex(el); cb(true); });
    },
    unsubscribe: function (el) { $(el).off('.mmColor'); },
    receiveMessage: function (el, data) {
      if (data.value) { el.value = data.value; MM.syncHex(el); $(el).trigger('change'); }
    },
    getRatePolicy: function () { return { policy: 'debounce', delay: 400 }; }
  });
  Shiny.inputBindings.register(colorBinding, 'multimapr.color');

  MM.syncHex = function (el) {
    var target = el.getAttribute('data-hex-target');
    if (target) { var t = document.getElementById(target); if (t) t.textContent = el.value.toUpperCase(); }
  };

  /* ---- i18n ------------------------------------------------------------------ */
  MM.t = function (key) {
    if (!MM.dict || !MM.dict[key]) return null;
    return MM.dict[key][MM.lang] || MM.dict[key].en || null;
  };
  MM.applyLang = function (root) {
    if (!MM.dict) return;
    root = root || document;
    $(root).find('[data-i18n]').each(function () {
      var v = MM.t(this.getAttribute('data-i18n')); if (v !== null) this.textContent = v;
    });
    $(root).find('[data-i18n-ph]').each(function () {
      var v = MM.t(this.getAttribute('data-i18n-ph')); if (v !== null) this.setAttribute('placeholder', v);
    });
    $(root).find('[data-i18n-title]').each(function () {
      var v = MM.t(this.getAttribute('data-i18n-title'));
      if (v !== null) { this.setAttribute('title', v); this.setAttribute('aria-label', v); }
    });
    document.documentElement.setAttribute('lang', MM.lang);
    $('.mm-lang button').each(function () {
      var on = this.getAttribute('data-lang') === MM.lang;
      this.classList.toggle('is-active', on);
      this.setAttribute('aria-pressed', on ? 'true' : 'false');
    });
    MM.syncThemeLabel();
  };
  MM.setLang = function (lang) {
    MM.lang = lang;
    try { localStorage.setItem('mm-lang', lang); } catch (e) {}
    MM.applyLang(document);
    Shiny.setInputValue('mmr_lang', lang);
  };
  Shiny.addCustomMessageHandler('mm-dict', function (msg) {
    MM.dict = msg.dict;
    MM.lang = msg.lang;
    MM.applyLang(document);
    Shiny.setInputValue('mmr_lang', MM.lang);
  });
  // Server-rendered fragments already come translated; static pieces inside them
  // (and shiny-bound placeholders) are refreshed here.
  $(document).on('shiny:value', function (e) {
    setTimeout(function () { MM.applyLang(e.target); }, 0);
  });

  /* ---- Theme ------------------------------------------------------------------ */
  MM.syncThemeLabel = function () {
    var dark = document.documentElement.getAttribute('data-bs-theme') === 'dark';
    $('.mm-theme-btn .mm-theme-label').each(function () {
      var key = dark ? 'theme_dark' : 'theme_light';
      this.setAttribute('data-i18n', key);
      var v = MM.t(key); if (v !== null) this.textContent = v;
    });
    $('.mm-theme-btn').attr('aria-pressed', dark ? 'true' : 'false');
  };
  MM.setTheme = function (theme) {
    document.documentElement.setAttribute('data-bs-theme', theme);
    try { localStorage.setItem('mm-theme', theme); } catch (e) {}
    MM.syncThemeLabel();
  };
  try {
    var saved = localStorage.getItem('mm-theme');
    if (saved === 'dark' || saved === 'light') document.documentElement.setAttribute('data-bs-theme', saved);
  } catch (e) {}

  $(document).on('click', '.mm-lang button', function () { MM.setLang(this.getAttribute('data-lang')); });
  $(document).on('click', '.mm-theme-btn', function () {
    var dark = document.documentElement.getAttribute('data-bs-theme') === 'dark';
    MM.setTheme(dark ? 'light' : 'dark');
  });

  /* ---- Character list ------------------------------------------------------------ */
  MM.selectedChars = [];
  MM.markSelected = function () {
    var sel = MM.selectedChars;
    $('.mm-charrow').each(function () {
      var on = sel.indexOf(this.getAttribute('data-char')) >= 0;
      this.classList.toggle('is-selected', on);
      this.setAttribute('aria-pressed', on ? 'true' : 'false');
    });
  };
  Shiny.addCustomMessageHandler('mm-chars', function (msg) {
    MM.selectedChars = [].concat(msg.sel || []);
    MM.markSelected();
  });
  $(document).on('shiny:value', function (e) {
    if (e.name === 'char_list') setTimeout(function () { MM.markSelected(); MM.filterChars(); }, 0);
  });
  $(document).on('click', '.mm-charrow, .mm-chip-x', function () {
    Shiny.setInputValue('char_toggle', this.getAttribute('data-char'), { priority: 'event' });
  });
  MM.filterChars = function () {
    var box = document.getElementById('char_search');
    var q = box ? box.value.trim().toLowerCase() : '';
    $('.mm-charrow').each(function () {
      var name = (this.getAttribute('data-char') || '').toLowerCase();
      this.hidden = q.length > 0 && name.indexOf(q) < 0;
    });
  };
  $(document).on('input', '#char_search', MM.filterChars);

  /* ---- Sortable statistics table ------------------------------------------------- */
  $(document).on('click', '.mm-table thead th[data-col]', function () {
    var th = this, table = th.closest('table'), tbody = table.tBodies[0];
    var col = parseInt(th.getAttribute('data-col'), 10);
    var numeric = th.getAttribute('data-type') === 'num';
    var dir = th.getAttribute('data-dir') === 'asc' ? 'desc' : 'asc';
    $(table).find('thead th').removeClass('is-sorted').removeAttr('data-dir').attr('aria-sort', 'none');
    th.classList.add('is-sorted'); th.setAttribute('data-dir', dir);
    th.setAttribute('aria-sort', dir === 'asc' ? 'ascending' : 'descending');
    $(th).find('.mm-sort').text(dir === 'asc' ? '\u25B2' : '\u25BC');
    var rows = Array.prototype.slice.call(tbody.rows);
    rows.sort(function (a, b) {
      var x = a.cells[col].getAttribute('data-v'), y = b.cells[col].getAttribute('data-v');
      var r;
      if (numeric) r = parseFloat(x) - parseFloat(y);
      else r = x.localeCompare(y, undefined, { numeric: true, sensitivity: 'base' });
      return dir === 'asc' ? r : -r;
    });
    rows.forEach(function (r) { tbody.appendChild(r); });
  });

  /* ---- Data-category colors (scored / missing / inapplicable) ------------------- */
  // Matrix cells, table bars and KPI accents read these CSS variables.
  Shiny.addCustomMessageHandler('mm-stats-colors', function (msg) {
    var st = document.documentElement.style;
    [['s', 'scored'], ['m', 'missing'], ['i', 'inapplicable']].forEach(function (p) {
      st.setProperty('--mm-c-' + p[0], msg[p[1]]);
      st.setProperty('--mm-c-' + p[0] + '-ink', msg[p[1] + '_ink']);
    });
    document.documentElement.classList.toggle('mm-no-borders', !msg.borders);
  });

  /* ---- Taxon x character matrix: row/column highlight + readout ----------------- */
  MM.mxCategory = { 'mm-c-s': 'col_scored', 'mm-c-m': 'col_missing', 'mm-c-i': 'col_inapp' };
  MM.mxClear = function (table) {
    $(table).find('.is-hl, .is-hl-cell').removeClass('is-hl is-hl-cell');
  };
  $(document).on('mouseover', '.mm-matrix td', function () {
    var td = this, row = td.parentNode, table = td.closest('table');
    var ci = td.cellIndex;
    MM.mxClear(table);
    row.classList.add('is-hl');
    var head = table.tHead.rows[0].cells[ci];
    head.classList.add('is-hl');
    Array.prototype.forEach.call(table.tBodies[0].rows, function (r) {
      if (r.cells[ci]) r.cells[ci].classList.add('is-hl');
    });
    td.classList.add('is-hl-cell');
    var info = document.getElementById('mx_info');
    if (info) {
      var cat = MM.t(MM.mxCategory[td.className.split(' ')[0]]) || '';
      info.textContent = row.cells[0].textContent + '  ·  ' + head.textContent +
        '  =  ' + td.textContent + (cat ? '   (' + cat + ')' : '');
    }
  });
  $(document).on('mouseleave', '.mm-matrix', function () { MM.mxClear(this); });

  /* ---- Preview canvas height ---------------------------------------------------- */
  MM.setCanvasHeight = function (h) {
    var c = document.getElementById('mm-canvas'); if (c && h) c.style.height = h + 'px';
  };
  $(document).on('shiny:inputchanged', function (e) {
    if (e.name === 'preview_height') MM.setCanvasHeight(e.value);
  });

  /* ---- Misc server messages -------------------------------------------------------- */
  // Shows the source of preloaded data inside a fileInput's text box
  Shiny.addCustomMessageHandler('mm-file-label', function (msg) {
    var el = document.getElementById(msg.id);
    if (!el) return;
    var box = $(el).closest('.input-group').find('input[type=text]');
    box.val(msg.text || '');
    if (!msg.text) {
      $(el).val('');
      $(el).closest('.shiny-input-container').find('.progress').css('visibility', 'hidden');
    }
  });

  Shiny.addCustomMessageHandler('mm-class', function (msg) {
    var el = document.getElementById(msg.id);
    if (el) { el.classList.toggle(msg.cls, !!msg.on); el.setAttribute('aria-pressed', msg.on ? 'true' : 'false'); }
  });

  // Logs stay scrolled to the newest line
  $(document).on('shiny:value', function (e) {
    var el = e.target;
    if (el && el.classList && el.classList.contains('mm-log')) {
      setTimeout(function () { el.scrollTop = el.scrollHeight; }, 0);
    }
  });

  $(document).on('shiny:connected', function () {
    $('input.mm-color').each(function () { MM.syncHex(this); });
    MM.syncThemeLabel();
  });
})();
