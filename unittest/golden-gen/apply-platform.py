#!/usr/bin/env python3
"""Losslessly expand frozen value overrides; never execute model/generator code."""
import csv
import sys

base, overlay, output = sys.argv[1:]
with open(base) as f:
    rows = list(csv.reader(line for line in f if not line.startswith("#")))
with open(overlay) as f:
    replacements = list(csv.reader(line for line in f if not line.startswith("#")))
if any(len(row) != 4 for row in replacements):
    raise ValueError("platform overrides may contain only dist,params,draw_index,value")
values = {tuple(row[:3]): row[3] for row in replacements}
if len(values) != len(replacements):
    raise ValueError("duplicate platform override")
for row in rows:
    row[3] = values.pop(tuple(row[:3]), row[3])
if values:
    raise ValueError("platform override is absent from the frozen schedule")
with open(output, "w", newline="") as f:
    csv.writer(f, lineterminator="\n").writerows(rows)
