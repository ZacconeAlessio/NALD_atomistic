#!/usr/bin/env python3
"""Validate the transcribed Kobayashi et al. JCP 2026 benchmark table.

For rows containing K_inf-K_0 and tau_zeta, the paper's Eq. (4) requires

    zeta = (K_inf-K_0) tau_zeta.

With K in GPa and tau in ps, K*tau is in mPa s.  The published zeta column is
in units of 10^-1 mPa s, so the table value should be approximately
10*K_GPa*tau_ps. Small differences are expected from printed rounding.
"""

from __future__ import annotations

import csv
from pathlib import Path


def main() -> None:
    path = Path(__file__).with_name("kobayashi_reference.csv")
    worst = (0.0, None)
    checked = 0
    with path.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row["model"] == "experiment" or not row["Kinf_minus_K0_GPa"]:
                continue
            zeta = float(row["zeta_1e-1_mPa_s"])
            kdiff = float(row["Kinf_minus_K0_GPa"])
            tau = float(row["tau_zeta_ps"])
            reconstructed = 10.0 * kdiff * tau
            rel = abs(reconstructed - zeta) / zeta
            checked += 1
            if rel > worst[0]:
                worst = (rel, row)
    if checked == 0:
        raise SystemExit("no benchmark rows checked")
    print(f"checked {checked} MD rows")
    print(f"worst relative discrepancy from zeta=(K_inf-K_0)tau: {worst[0]:.4%}")
    if worst[0] > 0.05:
        raise SystemExit("FAIL: benchmark transcription/formula mismatch exceeds 5%")
    print("PASS: benchmark table is internally consistent within printed rounding.")


if __name__ == "__main__":
    main()
