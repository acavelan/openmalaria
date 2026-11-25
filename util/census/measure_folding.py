#!/usr/bin/env python3
"""Pair-folding LOC measurement for the modules spike (header-DELETION upside).

A .h+.cpp -> single .cppm fold deletes outright, per category:
  1. in-class member prototypes in .h whose body lives in a .cpp (Class::m)
  2. include-guard scaffold (#ifndef + #define + closing #endif) / #pragma once
  3. #include directives (net of ~1 import line per consumer file)
  4. forward-declaration lines (class X; / struct X;)
Handles OM brace style ({ on own line). xml/ subtotal reported separately.
"""
import re, sys, pathlib

root = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "model")
headers = sorted(root.rglob("*.h"))
sources = sorted(root.rglob("*.cpp"))

def is_xml(p): return "xml/" in str(p).replace("\\", "/")
def code(s):
    s = s.strip()
    return s and not s.startswith("//") and not s.startswith("/*") and not s.startswith("*")

defs = set()
def_re = re.compile(r"^\s*(?:template\s*<[^>]*>\s*)?[A-Za-z_][\w:<>,\s\*&\[\]]*?\b(\w+)::(\w+)\s*\(")
for s in sources:
    for line in s.read_text(errors="replace").splitlines():
        m = def_re.match(line)
        if m and not line.lstrip().startswith("//"):
            defs.add((m.group(1), m.group(2)))

kw = {"if","for","while","switch","return","catch","do","else","delete","new","throw"}
proto = re.compile(r"^\s*(?:virtual\s+|static\s+|explicit\s+|constexpr\s+|friend\s+|inline\s+)*"
                   r"[\w:<>,\s\*&\[\]]*?\b(\w+)\s*\([^;{)]*\)\s*"
                   r"(?:const\s*)?(?:override\s*)?(?:noexcept\s*)?(?:=\s*(?:0|delete|default)\s*)?;\s*(?://.*)?$")
fwd_re = re.compile(r"^\s*(class|struct|union|enum(?:\s+class)?)\s+\w+\s*(?::\s*[\w,: ]+)?\s*;\s*(?://.*)?$")
cls_re = re.compile(r"^\s*(?:template\s*<[^>]*>\s*)?(?:class|struct|union)\s+(\w+)(?:\s+(?:final|override)\b|\s*:\s*[^{;]+)?\s*$")

cat1 = cat1x = 0
ex1 = []
n_guard = n_guardx = n_inc = n_incx = n_fwd = n_fwdx = 0

for h in headers:
    lines = h.read_text(errors="replace").splitlines()
    lead = [l for l in lines if code(l)]
    # guard scaffold: first code line #ifndef, last code line #endif, optional #define after it
    if lead and lead[0].lstrip().startswith("#ifndef"):
        n_guard += 1; n_guardx += is_xml(h)
        if len(lead) > 1 and lead[1].lstrip().startswith("#define"):
            n_guard += 1; n_guardx += is_xml(h)
        if lead[-1].lstrip().startswith("#endif"):
            n_guard += 1; n_guardx += is_xml(h)
    n_guard += sum(1 for l in lines if "#pragma once" in l)
    n_guardx += sum(is_xml(h) for l in lines if "#pragma once" in l)

    stack = []            # (class_name, depth_when_open)
    depth = 0
    pending = None
    for line in lines:
        if line.strip().startswith("#include"):
            n_inc += 1; n_incx += is_xml(h); continue
        if fwd_re.match(line):
            n_fwd += 1; n_fwdx += is_xml(h); continue
        cm = cls_re.match(line)
        if cm and ";" not in line.split("//")[0]:
            pending = cm.group(1)
        if "{" in line and pending:
            stack.append((pending, depth)); pending = None
        m = proto.match(line)
        if m and stack and "{" not in line:
            name = m.group(1)
            cname = stack[-1][0]
            if name not in kw and name != cname and (cname, name) in defs:
                cat1 += 1; cat1x += is_xml(h)
                if len(ex1) < 3: ex1.append(f"{h}:{line.strip()[:70]}")
        depth += line.count("{") - line.count("}")
        while stack and depth <= stack[-1][1]:
            stack.pop()

cpp_inc = sum(1 for s in sources for l in s.read_text(errors="replace").splitlines()
              if l.strip().startswith("#include"))
cpp_incx = sum(is_xml(s) for s in sources for l in s.read_text(errors="replace").splitlines()
               if l.strip().startswith("#include"))

print(f"files: {len(headers)} .h + {len(sources)} .cpp")
print(f"1 header prototypes w/ body in paired cpp : {cat1:5d}  (xml/: {cat1x})")
print(f"2 guard/#pragma scaffold lines            : {n_guard:5d}  (xml/: {n_guardx})")
print(f"3a #include lines in .h                   : {n_inc:5d}  (xml/: {n_incx})")
print(f"3b #include lines in .cpp                 : {cpp_inc:5d}  (xml/: {cpp_incx})")
print(f"4 forward-decl lines in .h                : {n_fwd:5d}  (xml/: {n_fwdx})")
print("examples of cat-1:")
for e in ex1: print("  ", e)
