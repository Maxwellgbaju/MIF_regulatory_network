#!/bin/bash
# Quick check that all table scripts are present and importable. Run: bash check_setup.sh
for f in tab_common.py table1_integrated_evidence.py table2_coexpression_sensitivity.py table3_survival_sensitivity.py immune_bridge_check.py; do
  [ -f "$f" ] && echo "OK   $f ($(wc -l < $f) lines)" || echo "MISSING $f"
done
python3 -c "import tab_common; print('OK   tab_common imports')" 2>&1 | tail -1
python3 -c "import lifelines; print('OK   lifelines installed')" 2>&1 | tail -1
