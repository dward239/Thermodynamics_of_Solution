# Vessel input editor

Open `index.html` from a local checkout in a modern desktop browser.
Use **Load example** to try the synthetic case or import a vessel text
file. The validation list links to the corresponding form fields.
Download the finished case and pass it to `vessel_cli`, for example:

```sh
./build/vessel_cli vessel_model.txt output.csv
```

On Windows with a Visual Studio build, the executable is commonly
`build\Debug\vessel_cli.exe`. The GUI edits the fixed-step vessel CLI format.
Adaptive cases can be prepared using the [vessel input reference](../../docs/vessel-input-format.md).

Files are read and generated in the browser. Case data stays in the open
page until downloaded. Export groups records by keyword and moves imported
comments to the top of the file.

## Extending the editor

- `model.js` defines record order, fields, units, bounds, cardinality,
  references, and rules that span records. The form, parser and exporter
  read these definitions.
- `editor.js` renders records from those definitions. It handles file
  selection, validation feedback, and downloads.
- When adding a supported record, update the C++ parser and the record
  definition together. Put domain rules that involve multiple records in
  the model validation function and exercise the exported file with
  `vessel_cli`.

Run the model checks with Node.js 18 or newer:

```sh
node --test tools/vessel_editor/model.test.cjs
```
