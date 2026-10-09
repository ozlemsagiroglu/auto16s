#!/usr/bin/env bash
# Turns failed Nextflow tasks into GitHub annotations, so the cause of a failed CI run is visible
# on the run page (and through the API) without downloading logs.
shopt -s nullglob
esc() { local s="${1//'%'/'%25'}"; s="${s//$'\r'/'%0D'}"; printf '%s' "${s//$'\n'/'%0A'}"; }
escp() { local s; s=$(esc "$1"); s="${s//':'/'%3A'}"; printf '%s' "${s//','/'%2C'}"; }
n=0
for f in work/*/*/.exitcode; do
  code=$(cat "$f")
  [ "$code" = "0" ] && continue
  d=$(dirname "$f")
  task=$(grep -m1 -oP "(?<=^### name: ').*(?=')|(?<=# NEXTFLOW TASK: ).*" "$d/.command.run" 2>/dev/null || echo "$d")
  err=$(tail -c 3500 "$d/.command.err" 2>/dev/null; echo; tail -c 1500 "$d/.command.out" 2>/dev/null)
  echo "::error title=$(escp "Failed task $task (exit $code)")::$(esc "$err")"
  n=$((n + 1))
  [ $n -ge 8 ] && break
done
if [ $n -eq 0 ] && [ -f .nextflow.log ]; then
  msg=$(grep -E "ERROR|Exception|Error" .nextflow.log | tail -n 25)
  echo "::error title=$(escp "Nextflow error (no failed task found)")::$(esc "$msg")"
fi
exit 0
