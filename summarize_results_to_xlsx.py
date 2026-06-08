#!/usr/bin/env python3
import csv
import glob
import math
import os
import re
import sys
import zipfile
from datetime import datetime, timezone
from xml.sax.saxutils import escape


HEADERS = [
    "instance",
    "B",
    "N",
    "R",
    "time_to_best_s",
    "best",
    "avg",
    "std",
    "best_seed",
    "runs",
    "verified_runs",
]


def find_latest_results_csv():
    candidates = sorted(glob.glob("batch_results_*/results.csv"), key=os.path.getmtime, reverse=True)
    if not candidates:
        raise FileNotFoundError("No batch_results_*/results.csv found in current directory")
    return candidates[0]


def instance_key(name):
    m = re.search(r"B(\d+)_N(\d+)_R(\d+)", name)
    if not m:
        return (10**9, 10**9, 10**9, name)
    return (int(m.group(1)), int(m.group(2)), int(m.group(3)), name)


def parse_instance(name):
    m = re.search(r"B(\d+)_N(\d+)_R(\d+)", name)
    if not m:
        return "", "", ""
    return m.group(1), m.group(2), m.group(3)


def fmt_time(value):
    return f"{value:.4f}".rstrip("0").rstrip(".")


def sample_std(values):
    n = len(values)
    if n <= 1:
        return 0.0
    avg = sum(values) / n
    return math.sqrt(sum((x - avg) ** 2 for x in values) / (n - 1))


def summarize(results_csv):
    groups = {}
    with open(results_csv, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        required = {"instance", "seed", "best_profit", "best_time", "verified"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{results_csv} missing required columns: {sorted(missing)}")
        for row in reader:
            groups.setdefault(row["instance"], []).append(row)

    rows = []
    for inst, items in groups.items():
        vals = []
        for item in items:
            vals.append(
                {
                    "seed": int(item["seed"]),
                    "profit": int(item["best_profit"]),
                    "time": float(item["best_time"]),
                    "verified": item["verified"],
                }
            )

        profits = [v["profit"] for v in vals]
        best = max(profits)
        avg = sum(profits) / len(profits)
        std = sample_std(profits)
        best_runs = [v for v in vals if v["profit"] == best]
        best_time = min(v["time"] for v in best_runs)
        best_seed = min(v["seed"] for v in best_runs if v["time"] == best_time)
        b, n, r = parse_instance(inst)

        rows.append(
            {
                "instance": inst,
                "B": int(b) if b else "",
                "N": int(n) if n else "",
                "R": int(r) if r else "",
                "time_to_best_s": float(fmt_time(best_time)),
                "best": best,
                "avg": round(avg, 2),
                "std": round(std, 2),
                "best_seed": best_seed,
                "runs": len(vals),
                "verified_runs": sum(1 for v in vals if v["verified"] == "YES"),
            }
        )

    rows.sort(key=lambda row: instance_key(row["instance"]))
    return rows


def col_name(index):
    name = ""
    while index:
        index, rem = divmod(index - 1, 26)
        name = chr(65 + rem) + name
    return name


def cell_xml(row_idx, col_idx, value, style=0):
    ref = f"{col_name(col_idx)}{row_idx}"
    style_attr = f' s="{style}"' if style else ""
    if isinstance(value, (int, float)) and value != "":
        return f'<c r="{ref}"{style_attr}><v>{value}</v></c>'
    return f'<c r="{ref}" t="inlineStr"{style_attr}><is><t>{escape(str(value))}</t></is></c>'


def build_sheet_xml(rows):
    xml_rows = []
    header_cells = [cell_xml(1, i + 1, h, 1) for i, h in enumerate(HEADERS)]
    xml_rows.append(f'<row r="1">{"".join(header_cells)}</row>')

    for r_idx, row in enumerate(rows, start=2):
        cells = []
        for c_idx, header in enumerate(HEADERS, start=1):
            style = 2 if header in {"avg", "std", "time_to_best_s"} else 0
            cells.append(cell_xml(r_idx, c_idx, row[header], style))
        xml_rows.append(f'<row r="{r_idx}">{"".join(cells)}</row>')

    last_row = len(rows) + 1
    last_col = col_name(len(HEADERS))
    cols = (
        '<cols>'
        '<col min="1" max="1" width="28" customWidth="1"/>'
        '<col min="2" max="4" width="8" customWidth="1"/>'
        '<col min="5" max="5" width="15" customWidth="1"/>'
        '<col min="6" max="8" width="13" customWidth="1"/>'
        '<col min="9" max="11" width="13" customWidth="1"/>'
        '</cols>'
    )
    return f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
  <sheetViews>
    <sheetView workbookViewId="0">
      <pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/>
    </sheetView>
  </sheetViews>
  {cols}
  <sheetData>
    {''.join(xml_rows)}
  </sheetData>
  <autoFilter ref="A1:{last_col}{last_row}"/>
</worksheet>'''


def write_xlsx(rows, output_xlsx):
    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    sheet_xml = build_sheet_xml(rows)
    files = {
        "[Content_Types].xml": '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>
  <Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>
  <Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
  <Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>
  <Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>
</Types>''',
        "_rels/.rels": '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
  <Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>
</Relationships>''',
        "docProps/app.xml": '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">
  <Application>Python</Application>
</Properties>''',
        "docProps/core.xml": f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <dc:title>GEO batch summary</dc:title>
  <dc:creator>summarize_results_to_xlsx.py</dc:creator>
  <dcterms:created xsi:type="dcterms:W3CDTF">{now}</dcterms:created>
  <dcterms:modified xsi:type="dcterms:W3CDTF">{now}</dcterms:modified>
</cp:coreProperties>''',
        "xl/workbook.xml": '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <sheets>
    <sheet name="summary" sheetId="1" r:id="rId1"/>
  </sheets>
</workbook>''',
        "xl/_rels/workbook.xml.rels": '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>
</Relationships>''',
        "xl/styles.xml": '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
  <fonts count="2">
    <font><sz val="11"/><name val="Calibri"/></font>
    <font><b/><sz val="11"/><name val="Calibri"/></font>
  </fonts>
  <fills count="2">
    <fill><patternFill patternType="none"/></fill>
    <fill><patternFill patternType="gray125"/></fill>
  </fills>
  <borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>
  <cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>
  <cellXfs count="3">
    <xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>
    <xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/>
    <xf numFmtId="2" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/>
  </cellXfs>
  <cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>
</styleSheet>''',
        "xl/worksheets/sheet1.xml": sheet_xml,
    }

    os.makedirs(os.path.dirname(os.path.abspath(output_xlsx)) or ".", exist_ok=True)
    with zipfile.ZipFile(output_xlsx, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for name, content in files.items():
            zf.writestr(name, content)


def main():
    if len(sys.argv) > 3:
        raise SystemExit("Usage: python3 summarize_results_to_xlsx.py [results.csv] [output.xlsx]")

    results_csv = sys.argv[1] if len(sys.argv) >= 2 else find_latest_results_csv()
    if len(sys.argv) >= 3:
        output_xlsx = sys.argv[2]
    else:
        output_xlsx = os.path.join(os.path.dirname(results_csv), "summary_best_avg_std.xlsx")

    rows = summarize(results_csv)
    write_xlsx(rows, output_xlsx)

    run_counts = sorted({row["runs"] for row in rows})
    verified_bad = [row for row in rows if row["runs"] != row["verified_runs"]]

    print(f"Input: {results_csv}")
    print(f"Output: {output_xlsx}")
    print(f"Instances: {len(rows)}")
    print("Runs per instance: " + ", ".join(str(x) for x in run_counts))
    if verified_bad:
        print(f"WARNING: {len(verified_bad)} instances have unverified runs")


if __name__ == "__main__":
    main()
