// The same definitions drive the form, file format, and field validation.
(() => {
  const token = (label, source) => ({ label, source });
  const number = (label, unit, bounds = {}) => ({
    label, unit, numeric: true, min: 0, ...bounds,
  });
  const fraction = (label, bounds) => number(label, '', { max: 1, ...bounds });
  const positive = { exclusiveMin: true };
  const nuclide = label => token(label, 'nuclide');
  const groups = [
    { id: 'decay', label: 'Decay & inventory' },
    { id: 'sources', label: 'Energy & radiolysis' },
    { id: 'vessel', label: 'Vessel & run' },
  ];
  const records = {
    nuclide: { group: 'decay', label: 'Nuclides', min: 1, unique: [0], provides: 0,
      hint: 'Use a zero half-life for a stable nuclide.',
      fields: [token('Nuclide ID'), token('Element'), number('Half-life', 's'), number('Emitted energy', 'J/decay')] },
    daughter: { group: 'decay', label: 'Daughter branches', unique: [0, 1],
      hint: 'Branch fractions for each radioactive parent must total 1.',
      fields: [nuclide('Parent'), nuclide('Daughter'), fraction('Branch fraction', positive)] },
    atoms: { group: 'decay', label: 'Initial inventory', min: 1, unique: [0],
      hint: 'Omitted nuclides start with zero atoms.',
      fields: [nuclide('Nuclide'), number('Initial inventory', 'atoms')] },
    deposit: { group: 'sources', label: 'Energy deposition', unique: [0, 1], provides: 1,
      hint: 'Deposition fractions for each nuclide may total up to 1.',
      fields: [nuclide('Source nuclide'), token('Phase ID'), fraction('Deposited fraction')] },
    gvalue: { group: 'sources', label: 'Radiolysis yields',
      hint: 'These yields contribute to the gross products column in the CSV.',
      fields: [token('Phase', 'deposit'), token('Product ID'), number('G-value', 'molecules/100 eV')] },
    thermal: { group: 'vessel', label: 'Thermal properties', min: 1, singleton: true,
      hint: 'Set heat loss conductance to zero for an insulated vessel.',
      fields: [number('Thermal mass', 'kg', positive), number('Specific heat', 'J/(kg·K)', positive),
        number('Heat loss conductance', 'W/K'), number('Ambient temperature', 'K', positive), fraction('Sensible heat fraction')] },
    headspace: { group: 'vessel', label: 'Headspace', min: 1, singleton: true,
      hint: 'Gas inventory and volume are held constant in the vessel CLI.',
      fields: [number('Volume', 'm³', positive), number('Initial gas', 'mol'), number('Starting temperature', 'K', positive)] },
    simulate: { group: 'vessel', label: 'Time grid', min: 1, singleton: true,
      hint: 'The run duration is the time step multiplied by the number of steps.',
      fields: [number('Fixed time step', 's', positive), number('Number of steps', '', { min: 1, max: 100000, integer: true })] },
  };
  const entries = Object.entries(records);
  const referenceSources = entries.filter(([, spec]) => spec.provides !== undefined).map(([key]) => key);
  const example = `# SYNTHETIC TEST ONLY. X and Y are fictitious and 2 J/decay is not physical.
nuclide X SyntheticElementX 10 2
nuclide Y SyntheticElementY 0 0
daughter X Y 1
deposit X water 0.5
atoms X 8
gvalue water H2 1
thermal 2 500 0 298.15 0.5
headspace 0.01 1 298.15
simulate 10 3
`;
  const blankRow = spec => spec.fields.map(() => '');
  const empty = () => Object.fromEntries(entries.map(([key, spec]) =>
    [key, Array.from({ length: spec.min ?? 0 }, () => blankRow(spec))]));

  function parse(text) {
    const data = Object.fromEntries(entries.map(([key]) => [key, []]));
    const comments = [];
    const errors = [];
    text.split(/\r?\n/).forEach((line, index) => {
      const [record, ...comment] = line.split('#');
      if (comment.length) comments.push('#' + comment.join('#'));
      if (!record.trim()) return;
      const [key, ...values] = record.trim().split(/\s+/);
      const spec = Object.hasOwn(records, key) ? records[key] : undefined;
      if (!spec) errors.push(`Line ${index + 1}: unknown keyword "${key}".`);
      else if (values.length !== spec.fields.length) errors.push(`Line ${index + 1}: ${key} needs ${spec.fields.length} values.`);
      else if (spec.singleton && data[key].length) errors.push(`Line ${index + 1}: duplicate ${key} record.`);
      else data[key].push(values);
    });
    if (errors.length) return { errors };
    for (const [key, spec] of entries) {
      while (data[key].length < (spec.min ?? 0)) data[key].push(blankRow(spec));
    }
    return { data, comments, errors };
  }

  const identifiers = (data, source) => new Set(
    data[source].map(row => row[records[source].provides].trim()).filter(Boolean)
  );
  const decimal = /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i;
  function validate(data) {
    const issues = [];
    const add = (key, row, field, message) => issues.push({ key, row, field, message });
    const sources = Object.fromEntries(referenceSources.map(id => [id, identifiers(data, id)]));
    for (const [key, spec] of entries) {
      const seen = new Set();
      data[key].forEach((row, index) => {
        spec.fields.forEach((field, column) => {
          const value = row[column].trim();
          const report = message => add(key, index, column, `${field.label} ${message}`);
          if (!value) report('is required.');
          else if (!field.numeric) {
            if (/[\s#]/.test(value)) report('must be a single token without spaces or #.');
            else if (field.source && !sources[field.source].has(value)) report('must reference a declared entry.');
          } else {
            const n = Number(value);
            if (!(field.integer ? /^\+?\d+$/.test(value) : decimal.test(value)) || !Number.isFinite(n))
              report(field.integer ? 'must be a whole number.' : 'must be a finite decimal number.');
            else if (field.exclusiveMin ? n <= field.min : n < field.min)
              report(`must be ${field.exclusiveMin ? 'greater than' : 'at least'} ${field.min}.`);
            else if (field.max !== undefined && n > field.max) report(`must be at most ${field.max}.`);
          }
        });
        if (spec.unique && spec.unique.every(column => row[column].trim())) {
          const id = JSON.stringify(spec.unique.map(column => row[column].trim()));
          if (seen.has(id)) add(key, index, spec.unique[0], 'This record duplicates an earlier entry.');
          seen.add(id);
        }
      });
    }
    const branchTotals = new Map();
    const depositTotals = new Map();
    for (const [id, , fraction] of data.daughter)
      branchTotals.set(id.trim(), (branchTotals.get(id.trim()) ?? 0) + Number(fraction));
    for (const [id, , fraction] of data.deposit)
      depositTotals.set(id.trim(), (depositTotals.get(id.trim()) ?? 0) + Number(fraction));
    const step = Number(data.simulate[0][0]);
    const count = Number(data.simulate[0][1]);
    if (step > 0 && count > 0 && !Number.isFinite(step * count))
      add('simulate', 0, 0, 'The total simulation duration exceeds the numeric range.');
    data.nuclide.forEach(([rawId, , halfLife], index) => {
      const id = rawId.trim();
      if (!id || !halfLife.trim()) return;
      const life = Number(halfLife);
      if (life === 0 && branchTotals.has(id)) add('nuclide', index, 2, `${id} is stable and has daughter branches.`);
      if (life > 0 && Math.abs((branchTotals.get(id) ?? 0) - 1) > 1e-10)
        add('nuclide', index, 2, `Daughter branch fractions for ${id} must sum to 1.`);
      if (life > 0) {
        const rate = Math.LN2 / life;
        if (!Number.isFinite(rate)) add('nuclide', index, 2, `${id} has a half-life below the solver’s numeric range.`);
        else if (step * rate > 120000)
          add('simulate', 0, 0, `Shorten the time step for the half-life of ${id}.`);
      }
      if ((depositTotals.get(id) ?? 0) > 1 + 1e-10)
        add('deposit', data.deposit.findIndex(row => row[0].trim() === id), 2,
          `Energy deposition fractions for ${id} exceed 1.`);
    });
    return issues;
  }

  function serialize(data, comments = []) {
    const lines = entries.flatMap(([key]) => data[key].filter(row => row.every(value => value.trim()))
      .map(row => `${key} ${row.map(value => value.trim()).join(' ')}`));
    return [...comments, ...lines].join('\n') + '\n';
  }

  globalThis.VesselEditorModel = { groups, records, entries, referenceSources, example, blankRow, empty, parse, identifiers, validate, serialize };
})();
