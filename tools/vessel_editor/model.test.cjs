const assert = require('node:assert/strict');
const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const { test } = require('node:test');
require('./model.js');

const model = globalThis.VesselEditorModel;
const example = () => model.parse(model.example);

test('the browser preset round-trips the bundled C++ vessel example', () => {
  const original = model.parse(readFileSync(join(__dirname, '../../examples/synthetic_vessel.txt'), 'utf8'));
  assert.deepEqual(original.errors, []);
  assert.deepEqual(example(), original);
  assert.deepEqual(model.validate(original.data), []);
  assert.deepEqual(model.parse(model.serialize(original.data, original.comments)), original);
});

test('bad syntax rejects the import and identifies the source line', () => {
  assert.deepEqual(model.parse('# heading\nunknown X\nthermal 2\n').errors, [
    'Line 2: unknown keyword "unknown".',
    'Line 3: thermal needs 5 values.',
  ]);
  assert.match(model.parse(model.example + 'simulate 1 2\n').errors[0], /duplicate simulate/);
});

test('references and cross-record constraints prevent invalid export', () => {
  const { data } = example();
  data.daughter[0][2] = '.8';
  data.deposit.push(['X', 'organic', '.6']);
  data.gvalue[0][0] = 'missing';
  data.atoms[0][0] = 'unknown';
  data.simulate[0][1] = '1.5';
  const issues = model.validate(data);
  assert.ok(issues.some(issue => issue.key === 'nuclide' && /sum to 1/.test(issue.message)));
  assert.ok(issues.some(issue => issue.key === 'deposit' && /exceed 1/.test(issue.message)));
  for (const key of ['gvalue', 'atoms'])
    assert.ok(issues.some(issue => issue.key === key && /reference/.test(issue.message)));
  assert.ok(issues.some(issue => issue.key === 'simulate' && /whole number/.test(issue.message)));
});

test('a missing required record imports as a field the user can repair', () => {
  const parsed = model.parse(model.example.replace(/^headspace .*\n/m, ''));
  assert.deepEqual(parsed.errors, []);
  assert.ok(model.validate(parsed.data).some(issue => issue.key === 'headspace'));
  parsed.data.headspace[0] = ['1e-2', '1', '298.15'];
  assert.deepEqual(model.validate(parsed.data), []);
});
