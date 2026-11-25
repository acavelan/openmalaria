#!/bin/sh
# draws.sh — regenerate the draw-site census and re-apply `// parity:` comments.
# Safe to re-run on an annotated tree (annotate.py is idempotent). From repo root:
#   ./util/census/draws.sh
set -e
cd "$(dirname "$0")/../.."
python3 util/census/draws.py
python3 util/census/annotate.py "$@"
# comment-only proof: every diff hunk added line must be `// parity:` or blank,
# and nothing may be removed
git diff -U0 -- model/ | awk '
/^-/   && !/^---/ { rem++ }
/^\+/  && !/^\+\+\+/ { add++;
    line=$0; sub(/^\+/,"",line); gsub(/^[ \t]+|[ \t]+$/,"",line);
    if (line != "" && line !~ /^\/\/ parity:/) { print "NON-COMMENT:" line; bad++ } }
END    { printf "added=%d removed=%d bad=%d %s\n", add, rem, bad,
         (rem==0 && bad==0) ? "COMMENT-ONLY:PASS" : "COMMENT-ONLY:FAIL"; exit bad+rem }'
