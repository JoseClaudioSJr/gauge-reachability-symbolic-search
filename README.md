# Hypergeometric Sums Lab (WZ / Zeilberger-style)

Goal: turn numerical conjectures into identities with recurrence-based certificates.

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Available Families

- `binom-sum`: \(\sum_{k=0}^{n} \binom{n}{k} = 2^n\)
- `binom-square`: \(\sum_{k=0}^{n} \binom{n}{k}^2 = \binom{2n}{n}\)
- `weighted-binom`: \(\sum_{k=0}^{n} k\binom{n}{k} = n2^{n-1}\)
- `franel-cube`: \(\sum_{k=0}^{n} \binom{n}{k}^3\) (Franel, order-2 recurrence)

## 1) WZ Certificate for One Family

Run:

```bash
python3 src/wz_lab.py prove-family --family binom-square
```

The script:

- automatically searches for a rational certificate `R(n,k)` in an ansatz space,
- symbolically verifies the residual of the WZ equation.

For `franel-cube`, the current certificate searcher (first-order normalized mode) returns `n/a`.
Now `prove-family` automatically tries a higher-order telescoping mode using the family's known recurrence.

## 2) Conjecture -> Recurrence

Run:

```bash
python3 src/wz_lab.py guess-family --family franel-cube --n-start 0 --n-end 30 --order -1 --degree 2
```

This command:

- generates exact terms of the sum,
- attempts to guess a linear recurrence with polynomial coefficients,
- validates on additional exact terms.

Note: `--order -1` uses the family-recommended order.

## 3) Automatic Report (3 Families)

Run:

```bash
python3 src/wz_lab.py report --n-end 30 --degree 2 --max-degree 4
```

Output:

- status of the conjectured recurrence,
- recurrence order used per family,
- status of the found WZ certificate,
- symbolic verification of the residual.

## Notes

- The `prove-family` command generates a symbolic WZ-style certificate for the selected family.
- The `guess-family` command is useful for discovering candidate recurrences (conjecture pipeline).
- The `report` command automates the cycle across all internal families.
- Use `--no-auto-order` in `report` to force a single order for all families.
- For families with recurrence order > 1 (e.g., `franel-cube`), the system attempts a higher-order telescoping certificate via rational ansatz.
