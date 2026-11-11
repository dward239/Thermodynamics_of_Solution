(() => {
  const model = globalThis.VesselEditorModel;
  const byId = id => document.getElementById(id);
  const make = (tag, value, className) => {
    const node = document.createElement(tag);
    if (value) node.textContent = value;
    if (className) node.className = className;
    return node;
  };
  const button = (label, action, className) => {
    const node = make('button', label, className);
    node.type = 'button';
    node.addEventListener('click', action);
    return node;
  };

  let data = model.empty();
  let comments = [];
  let section = model.groups[0].id;
  let dirty = false;
  let issues = [];
  const fieldId = (key, row, column) => `field-${key}-${row}-${column}`;

  function notify(message, error = false) {
    const notice = byId('notice');
    notice.hidden = !message;
    notice.classList.toggle('error', error);
    notice.textContent = message;
  }

  function renderSection() {
    for (const groupButton of byId('sections').children)
      groupButton.setAttribute('aria-pressed', String(groupButton.dataset.section === section));
    const cards = model.entries.filter(([, spec]) => spec.group === section).map(([key, spec]) => {
      const card = make('section', '', 'record');
      const head = make('div', '', 'record-head');
      const copy = make('div');
      copy.append(make('h2', spec.label), make('p', spec.hint));
      head.append(copy);
      if (!spec.singleton) {
        const add = button('Add row', () => {
          data[key].push(model.blankRow(spec));
          dirty = true;
          renderSection();
          byId(fieldId(key, data[key].length - 1, 0)).focus();
        });
        add.setAttribute('aria-label', `Add row to ${spec.label.toLowerCase()}`);
        head.append(add);
      }
      card.append(head);
      const rows = make('div', '', 'rows');
      data[key].forEach((row, index) => {
        const container = make('div', '', 'record-row');
        spec.fields.forEach((field, column) => {
          const group = make('div', '', 'field');
          const input = make('input');
          input.id = fieldId(key, index, column);
          input.type = 'text';
          input.spellcheck = false;
          input.autocomplete = 'off';
          if (field.numeric) input.inputMode = field.integer ? 'numeric' : 'decimal';
          if (field.source) input.setAttribute('list', `${field.source}-options`);
          input.value = row[column];
          input.addEventListener('input', () => {
            row[column] = input.value;
            dirty = true;
            refresh();
          });
          const label = make('label', field.label);
          label.htmlFor = input.id;
          if (field.unit) label.append(make('span', ` · ${field.unit}`, 'unit'));
          group.append(label, input);
          container.append(group);
        });
        if (!spec.singleton && data[key].length > (spec.min ?? 0)) {
          const remove = button('Remove', () => {
            data[key].splice(index, 1);
            dirty = true;
            renderSection();
            const nextRow = Math.min(index, data[key].length - 1);
            (byId(fieldId(key, nextRow, 0)) ?? cardAddButton(key)).focus();
          }, 'remove');
          remove.setAttribute('aria-label', `Remove ${spec.label.toLowerCase()} row ${index + 1}`);
          container.append(remove);
        }
        rows.append(container);
      });
      card.dataset.record = key;
      card.append(rows);
      return card;
    });
    byId('fields').replaceChildren(...cards);
    refresh();
  }

  function cardAddButton(key) {
    return byId('fields').querySelector(`[data-record="${key}"] button`);
  }

  function refresh() {
    for (const source of model.referenceSources) {
      const list = byId(`${source}-options`);
      const values = [...model.identifiers(data, source)];
      if (list.options.length !== values.length || values.some((value, i) => list.options[i].value !== value)) {
        list.replaceChildren(...values.map(value => {
          const option = make('option');
          option.value = value;
          return option;
        }));
      }
    }
    issues = model.validate(data);
    for (const input of byId('fields').querySelectorAll('input')) input.removeAttribute('aria-invalid');
    for (const issue of issues) byId(fieldId(issue.key, issue.row, issue.field))?.setAttribute('aria-invalid', 'true');
    const issueList = byId('issues');
    issueList.replaceChildren(...issues.map(issue => {
      const li = make('li');
      const spec = model.records[issue.key];
      li.append(button(`${spec.label}: ${issue.message}`, () => {
        if (section !== spec.group) {
          section = spec.group;
          renderSection();
        }
        byId(fieldId(issue.key, issue.row, issue.field))?.focus();
      }));
      return li;
    }));
    const status = byId('status');
    const statusText = issues.length ? `${issues.length} to resolve` : 'Ready to export';
    if (status.textContent !== statusText) status.textContent = statusText;
    status.classList.toggle('ready', !issues.length);
    const [step, count] = data.simulate[0].map(Number);
    byId('duration').textContent = step > 0 && count > 0 && Number.isFinite(step * count)
      ? `Run duration: ${(step * count).toLocaleString(undefined, { maximumSignificantDigits: 8 })} s` : 'Enter a time step and step count.';
    byId('preview').textContent = model.serialize(data, comments).trimEnd() || 'Your input records will appear here.';
    byId('export').disabled = Boolean(issues.length);
  }

  function replace(parsed, filename, message) {
    if (parsed.errors.length) {
      notify(parsed.errors.join('\n'), true);
      return;
    }
    if (dirty && !window.confirm('Replace the current case and discard its edits?')) return;
    data = parsed.data;
    comments = parsed.comments;
    dirty = false;
    section = model.groups[0].id;
    byId('filename').value = filename;
    notify(message);
    renderSection();
  }

  for (const source of model.referenceSources) {
    const list = make('datalist');
    list.id = `${source}-options`;
    byId('fields').after(list);
  }
  for (const group of model.groups) {
    const control = button(group.label, () => {
      section = group.id;
      renderSection();
    });
    control.dataset.section = group.id;
    byId('sections').append(control);
  }
  byId('new').addEventListener('click', () =>
    replace({ data: model.empty(), comments: [], errors: [] }, 'vessel_model.txt', 'New case.'));
  byId('example').addEventListener('click', () =>
    replace(model.parse(model.example), 'synthetic_vessel.txt', 'Loaded the synthetic example.'));
  byId('open').addEventListener('click', () => byId('file').click());
  byId('file').addEventListener('change', async event => {
    const file = event.target.files[0];
    event.target.value = '';
    if (!file) return;
    try {
      replace(model.parse(await file.text()), file.name, `Imported ${file.name}.`);
    } catch {
      notify('Unable to read the selected file.', true);
    }
  });
  byId('export').addEventListener('click', () => {
    if (issues.length) return;
    let name = byId('filename').value.trim() || 'vessel_model.txt';
    if (!/\.txt$/i.test(name)) name += '.txt';
    const url = URL.createObjectURL(new Blob([model.serialize(data, comments)], { type: 'text/plain;charset=utf-8' }));
    const link = make('a');
    link.href = url;
    link.download = name;
    link.click();
    setTimeout(() => URL.revokeObjectURL(url), 30000);
    dirty = false;
    notify(`Downloaded ${name}.`);
  });
  renderSection();
})();
