#!/usr/bin/env bash
# =============================================================================
# scan_private.sh - Block private material from the public repository
# =============================================================================
# Usage: bash scripts/qa/scan_private.sh      (also run by .githooks/pre-commit)
# Exit 0 = clean, 1 = finding.
#
# Scans every file git would publish for credentials, e-mail addresses, local
# or HPC paths, cluster account IDs, reviewer names and manuscript-handling
# terms, and for files larger than 50 MB.
set -uo pipefail
cd "$(git rev-parse --show-toplevel)" || exit 1

SECRETS='BEGIN [A-Z ]*PRIVATE KEY|api[_-]?key[[:space:]]*[:=]|password[[:space:]]*[:=]|secret[[:space:]]*[:=]|token[[:space:]]*[:=][[:space:]]*[A-Za-z0-9]{16,}|ghp_[A-Za-z0-9]{20,}|AKIA[0-9A-Z]{16}|WEBIN_(USERNAME|PASSWORD)'
EMAILS='[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}'
LOCALPATH='/Users/[a-zA-Z]|/home/[a-zA-Z]|C:\\\\Users|/crex/|/gorilla/|/proj/[a-z]|/sw/|OneDrive|Dropbox|iCloud'
ACCOUNTS='uppmax20[0-9]{2}-[0-9]+-[0-9]+|naiss20[0-9]{2}-[0-9]+-[0-9]+|snic20[0-9]{2}-[0-9]+-[0-9]+|fdube|farukmeddie'
PEOPLE='[Rr]eviewer|Doyle|Laing|Rockman|[Rr]ebuttal|[Cc]over letter'

# The scanner and the hook hold these patterns as literals. The UPPMAX profile
# binds the cluster's top-level mount points (/proj, /gorilla), which carry no
# project or user identifiers.
SELF_EXEMPT='^scripts/qa/scan_private\.sh$|^\.githooks/'
PATH_EXEMPT='^pipeline/conf/uppmax\.config$'

fail=0
report() { printf '  %-8s %s\n' "$1" "$2"; fail=1; }
n=0
while IFS= read -r f; do
  [ -f "$f" ] || continue
  size=$(wc -c < "$f")
  [ "$size" -gt 52428800 ] && report "SIZE" "$f ($size bytes > 50 MB)"
  printf '%s' "$f" | grep -qE "$SELF_EXEMPT" && continue
  case "$f" in
    *.gz) content() { gzip -dc "$f" 2>/dev/null; } ;;
    *.tif|*.png|*.pdf) content() { strings "$f" 2>/dev/null; } ;;
    *) content() { cat "$f"; } ;;
  esac
  n=$((n + 1))
  content | grep -qiE "$SECRETS"  && report "SECRET" "$f"
  # Random byte runs in images can look like addresses; check those for the rest only.
  case "$f" in *.tif|*.png|*.pdf) ;; *)
  content | grep -qE  "$EMAILS"   && report "EMAIL" "$f -> $(content | grep -oE "$EMAILS" | head -1)" ;; esac
  content | grep -qE  "$ACCOUNTS" && report "ACCOUNT" "$f -> $(content | grep -oE "$ACCOUNTS" | head -1)"
  content | grep -qE  "$PEOPLE"   && report "PEOPLE" "$f -> $(content | grep -oE "$PEOPLE" | head -1)"
  if ! printf '%s' "$f" | grep -qE "$PATH_EXEMPT"; then
    content | grep -qE "$LOCALPATH" && report "PATH" "$f -> $(content | grep -oE "$LOCALPATH" | head -1)"
  fi
done < <(git ls-files --cached --others --exclude-standard)

echo "scanned $n file(s) git would publish"
[ "$fail" -ne 0 ] && { echo "FAIL - see findings above."; exit 1; }
echo "PASS - no secrets, addresses, private paths, accounts or reviewer material."
