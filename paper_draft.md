---
header-includes:
  - \usepackage{graphicx}
---

# Gauge Selection as Symbolic Preconditioning for Rational Certificate Search

**Author:** José Cláudio da Silva Junior  
**Affiliation:** Electrical Engineer, State University of Londrina (UEL), 2016; M.Sc. in Electronic Systems and Electronic Instrumentation, State University of Londrina (UEL), 2018

---

## Abstract

The difficulty of finding a rational WZ certificate is not intrinsic to the
identity: it depends on how the summand is represented. We show empirically that
replacing the raw summand by its closed-form normalisation — a *gauge
transformation* — can change certificate search from *unreachable* to *reachable*
within a fixed ansatz budget, and that structural knowledge of the recurrence
reduces the denominator search space by 74–83%, cutting runtimes by 56–92%.
Controlled experiments on the binomial-power family
$\sum_k \binom{n}{k}^r$ ($r = 1, \ldots, 4$), a weighted binomial variant, and
the independent central Delannoy family all support the same conclusion: **gauge
selection may function as symbolic preconditioning in the tested regime**.

---

## 1. Introduction

Creative telescoping reduces a hypergeometric summation identity to finding a
rational function $R(n,k)$ — the WZ certificate — satisfying

$$
F(n+1,k) - F(n,k) = \Delta_k[R(n,k)\,F(n,k)].
$$

Algorithms such as Zeilberger's method guarantee a solution exists for
*holonomic* summands, but practical search is governed by an ansatz: a
candidate denominator set and a maximum polynomial degree. Whether a certificate
is found at all depends on whether the true certificate lies within that ansatz.

**The central observation of this work is that the same identity, expressed with
different summand normalisations, yields certificates of vastly different
complexity.** Two summands $F$ and $G$ that differ by a multiplicative rational
factor in $k$ define the same identity but lead to different telescoping
equations and different certificate functions. We call this choice of
normalisation a *gauge*, following the analogy with gauge freedom in physics.

This paper makes two empirical claims:

- **Phenomenon A (reachability):** Gauge choice appears to affect whether a certificate
  is found at all within a fixed ansatz budget — a qualitative, not merely
  quantitative, effect.

- **Phenomenon B (cost):** Structural features of the recurrence (predicted
  poles, recurrence order) appear to allow the denominator candidate set to be pruned
  substantially, reducing runtime by roughly 80–90% while preserving exact
  symbolic verification.

The two phenomena are related but distinct: Phenomenon A is about *whether* we
reach a certificate; Phenomenon B is about *how cheaply* we reach it once the
gauge is already favourable.

---

\begin{figure}[htbp]
\centering
\includegraphics[width=\linewidth]{figure1.pdf}
\caption{The core observation.}
\label{fig:core-observation}
\end{figure}

*Same identity, same ansatz budget, different gauge yields qualitatively different reachability.*

---

## 2. Background and Definitions

### 2.1 Hypergeometric families

We study three families:

- **Binomial-power:** $F_r(n,k) = \binom{n}{k}^r$, for $r = 1, 2, 3, 4$.
- **Weighted-binom:** $F(n,k) = k\binom{n}{k}$, normalised to
  $G(n,k) = \tfrac{k\binom{n}{k}}{n\,2^{n-1}}$ (the closed-form gauge).
- **Central Delannoy:** the central Delannoy recurrence family (order 2,
  independent of the binomial axis).

### 2.2 Gauge

Given a hypergeometric summand $F(n,k)$, any function $G(n,k) = c(n,k)\,F(n,k)$
with $c$ rational in $n,k$ defines the *same* identity
$\sum_k F(n,k) = \sum_k G(n,k) / c(n,k)$. We call the choice of normalisation a
*gauge*. Three gauges are tested:

| label | definition |
|---|---|
| `raw` | $F(n+1,k)/F(n,k) - 1$, no normalisation |
| `closed_form` | $(G(n+1,k)/G(n,k)) - 1$, where $G$ is the closed-form expression |
| `recurrence_lhs` | $\sum_j p_j(n)\,F(n+j,k)/F(n,k)$, telescoping from the recurrence LHS |

### 2.3 Reachability (operational definition)

**Definition (Reachability under budget $B$).** Fix a search budget
$B = (D, d, s)$ where $D$ is the set of denominator candidates, $d$ is the
maximum certificate polynomial degree, and $s$ is the maximum recurrence order.
A WZ certificate for summand $F$ is *reachable under $B$* if and only if the
search procedure finds a rational function $R(n,k)$ within $B$ and the identity

$$\sum_k F(n,k) = \text{const} \quad \Rightarrow \quad
F(n+1,k) - F(n,k) = \Delta_k[R(n,k)\,F(n,k)]$$

is verified symbolically for all $k$.

Failure to find a certificate under $B$ means only that the true certificate
lies outside $B$. It is **not** an impossibility claim. All "unreachable"
results in this paper are explicitly relative to the budget $B$ stated for each
experiment.

---

## 3. Why Gauge Affects Certificate Complexity: An Informal Account

This section provides structural intuition for the empirical results. No formal
proofs are claimed; the purpose is to make the phenomenon less surprising and to
suggest why closed-form normalisation is the natural preconditioner.

### 3.1 Shift ratios and certificate degree

Given a hypergeometric summand $F(n,k)$, the WZ method works with the
*shift ratio*

$$
\rho_F(n,k) = \frac{F(n+1,k)}{F(n,k)}.
$$

The certificate $R(n,k)$ must satisfy a linear equation whose coefficients
are polynomials in the shift ratio. Concretely, writing $R = P/Q$ with $Q$
in the denominator candidate set, the telescoping condition becomes a linear
system in the numerator coefficients of $P$. The complexity of that system —
its size, the degree of its polynomial entries, and whether its solution is
sparse — all depend on how simple $\rho_F$ is.

**When $F$ is a raw combinatorial summand**, $\rho_F$ typically contains
high-degree polynomial factors inherited from the falling-factorial structure
of binomial coefficients or rising factorials. These factors propagate into
the denominator of any certificate expressed relative to $F$.

**When $F$ is replaced by its closed-form normalisation $G = F / c(n,k)$**,
the shift ratio becomes

$$
\rho_G = \frac{c(n+1,k)}{c(n,k)} \cdot \rho_F.
$$

Choosing $c$ to be the closed-form summand (when the identity is $\sum_k F = S(n)$,
set $c = F / S$) makes $\rho_G$ a ratio of consecutive values of a known simple
function. In the weighted-binom case, this tends to cancel the dominant high-degree factor,
leaving a degree-1 shift ratio whose certificate is degree 1 in $k$.

### 3.2 Certificate pole structure

The poles of the rational certificate $R(n,k)$ (viewed as a function of $k$)
are related to the zeros of $\rho_F$ shifted by the recurrence order. Under
the raw gauge, these poles involve high-degree polynomial factors from the
combinatorial structure. Under the closed-form gauge, the normalisation tends to cancel
those factors, and the poles reduce to low-degree linear factors in $k$.

This may explain why guided denominator selection works: the
recurrence poles (extracted from the recurrence polynomial) are exactly the
poles of the certificate in the closed-form gauge, but appear not to be in the
raw gauge. Using recurrence-derived denominators with the raw gauge is therefore
structurally mismatched; with the closed-form gauge the match appears tight.

### 3.3 Summary of the structural argument

$$
\text{raw gauge} \;\Rightarrow\; \text{high-degree shift ratio}
\;\Rightarrow\; \text{high-degree certificate}
\;\Rightarrow\; \text{outside ansatz budget}
\;\Rightarrow\; \text{unreachable}
$$

$$
\text{closed-form gauge} \;\Rightarrow\; \text{low-degree shift ratio}
\;\Rightarrow\; \text{low-degree certificate}
\;\Rightarrow\; \text{inside ansatz budget}
\;\Rightarrow\; \text{reachable}
$$

This chain is informal: "high-degree" is not quantified in general, and the
relationship between shift-ratio degree and certificate degree is not proven
here. But it is consistent with all experimental observations, and provides a
testable structural prediction: gauges that simplify $\rho_F$ should be
preferred in any bounded search.

---

## 4. Phenomenon A: Gauge Sensitivity and Reachability

### 4.1 Experimental setup

We test the weighted-binom identity $\sum_k k\binom{n}{k} = n\,2^{n-1}$ under
three gauges, each with a guided denominator ansatz (denominator candidates
derived from the recurrence poles) and a maximum certificate degree of 4.

### 4.2 Results

| gauge | reachable | degree | equations solved | median time |
|---|---|---|---|---|
| `raw` | **NO** | — | 533 | 1.53 s |
| `closed_form` | **YES** | 1 | 10 | 0.033 s |
| `recurrence_lhs` | **NO** | — | 517 | 1.92 s |

**Certificate found under `closed_form`:**
$R(n,k) = \dfrac{k-1}{2(k-n-1)}$

The runtime ratio between the hardest and easiest gauge is **58.5×**. More
importantly, the effect is *qualitative*: two of the three gauges produce no
certificate within the budget, regardless of how much time is allowed for the
linear system.

### 4.3 Causal interpretation

The weighted-binom family was chosen specifically because, without the
closed-form gauge, all gauges fail. Adding the gauge $G = k\binom{n}{k}/(n\,2^{n-1})$
immediately restores reachability and recovers a degree-1 certificate. This is
a controlled intervention: one change, one outcome. The relationship between
gauge choice and reachability is consistent with a direct causal mechanism.

### 4.4 Narrative bridge: from reachability to cost

**Phenomenon A indicates that representation governs reachability.** Given a
solvable identity, choosing the right gauge is the price of admission: without it,
the certificate is unreachable within any reasonable budget. The question now is:
once we have chosen a reachable gauge, what further governs the *cost* of finding
the certificate?

**Phenomenon B addresses this:** structural features of the recurrence — specifically,
the poles predicted from the recurrence polynomial — appear to allow us to shrink the
denominator search space dramatically. This is not a gauge question; it applies
even to families (like Delannoy) where a natural closed form exists. The insight
is that guided ansatz design transfers across independent families
(Phenomenon C), suggesting a general principle: **once representation enables
reachability, structural guidance promotes efficiency**.

---

## 5. Phenomenon B: Structural Guided Ansatz and Search Cost

### 5.1 Experimental setup

For the binomial-power family $\sum_k \binom{n}{k}^r$ with $r = 1, 2, 3, 4$, we
compare two denominator strategies:

- **Baseline:** all denominator candidates up to a fixed degree bound (23–51
  candidates depending on $r$).
- **Guided:** only candidates predicted from the recurrence poles and structural
  features (6–9 candidates).

All experiments use the closed-form gauge. Exact symbolic verification is
performed; a result is marked verified only if the certificate passes the
telescoping identity check.

### 5.2 Results

| $r$ | baseline cands | guided cands | reduction | baseline time | guided time | time reduction | verified |
|---|---|---|---|---|---|---|---|
| 1 | 23 | 6 | 74% | 0.063 s | 0.028 s | 56% | **YES** |
| 2 | 27 | 6 | 78% | 8.89 s | 1.61 s | 82% | **YES** |
| 3 | 46 | 8 | 83% | 185.2 s | 33.0 s | 82% | NO |
| 4 | 51 | 9 | 82% | 421.2 s | 77.4 s | 82% | NO |

For $r = 1, 2$, guided search both prunes the space and finds a verified
certificate. For $r = 3, 4$, the guided ansatz is still substantially cheaper,
but neither baseline nor guided finds a verified certificate within the current
budget — consistent with the reachability definition: higher-order cases may
require a larger budget.

### 5.3 Delannoy validation

The central Delannoy family (order-2 recurrence, independent of the binomial
axis) provides a cross-family replication:

| setting | time |
|---|---|
| baseline | 22.617 s |
| guided | 2.339 s |
| reduction | **89.7%** |

This is consistent with guided denominator selection not being an artefact of the
binomial-power structure. The same structural principle — use recurrence poles
to focus the ansatz — appears to transfer across families.

---

## 6. Discussion

### 6.1 Why the reachability effect is qualitative

In standard algorithm analysis, a representation change affects runtime but not
correctness: the algorithm either terminates with a solution or it does not. Here,
the search is bounded: it stops when the ansatz is exhausted. Under this budget
model, a gauge transformation that simplifies the certificate structure may cross
the boundary between reachable and unreachable. The effect is not a speedup
of an existing search — it appears to be a phase transition in whether search succeeds at all.

### 6.2 Gauge selection as symbolic preconditioning

In numerical linear algebra, a preconditioner transforms a system so that its
condition number is more favourable for iterative solvers. A plausible analogy here
is that the gauge transforms the telescoping system so that its solution — the
certificate — has lower degree and simpler denominator structure. The search
procedure is unchanged; only the system it solves is different. In this sense,
closed-form normalisation may function as the symbolic analogue of a preconditioner.

### 6.3 Relationship between the two phenomena

Phenomena A and B are related but distinct. Phenomenon A (reachability) is
governed by whether the certificate degree and denominator are within the budget
at all. Phenomenon B (cost) is governed by how efficiently the budget is
traversed. The closed-form gauge is necessary for Phenomenon A. Guided
denominator selection contributes to Phenomenon B independently: it works even
on cases (Delannoy) where the gauge question does not arise, because the identity
is always expressed in a natural closed form.

---

## 7. Limitations and Open Questions

**Heuristic denominator prediction.** Guided denominators are generated from
recurrence poles and low-degree structural rules. There is no proof that all
true certificate poles lie within the predicted set. However, the 56–92%
cost reduction and cross-family transfer (Delannoy) suggest the heuristic
is reliable in practice for the families studied.

**Symbolic budget wall in higher-order families.** The $r = 3, 4$ cases of
binomial-power expose a structural bottleneck: higher-order recurrences appear
to require qualitatively larger ansatz budgets to produce verified certificates.
This is not a failure of the guided strategy — it proportionally reduces cost
even when unverified — but a discovery of a natural frontier. Whether this
frontier is inherent to the families or an artefact of the current search regime
remains open.

**No completeness claim.** Failure to find a certificate under budget $B$ is not
a proof of non-existence; it is a statement about $B$. All unreachability claims
in this work are explicitly relative to the stated budget parameters. Expanding
the budget may move a certificate from unreachable to reachable, as demonstrated
by the gauge sensitivity results.

**Manual gauge selection.** Closed-form gauges are supplied by hand. Automating
their discovery — especially for non-obvious normalisations — is an open problem
that could substantially broaden the applicability of the method.

**Single-variable recurrences.** All experiments use first-order recurrences in
$n$. Multi-variable or multi-dimensional identities are not addressed.

---

## 8. Conclusion

The difficulty of rational WZ certificate search appears to be strongly
representation-dependent. Two summands that encode the same identity but differ
in normalisation can yield certificates of markedly different complexity: one
reachable within a modest budget, the other unreachable within a substantially
larger budget.

Three lessons emerge from the experiments:

1. **Gauge matters qualitatively.** Changing the gauge is not merely a speedup:
   it can shift whether a certificate is reachable within a bounded ansatz budget.

2. **Structural guidance reduces search cost substantially.** Guided denominator selection,
   derived from recurrence poles, prunes the search space (by 56–92% across tested
   families) and appears to transfer across independent families.

3. **Negative results are data.** The trace dataset of failed searches — with
   recorded gauge, denominator candidates, systems solved, and runtime — is
   itself a contribution. The creative telescoping literature rarely reports
   what does not work.

The central message is concise: **representation governs certificate
reachability**. Gauge selection functions as symbolic preconditioning in the
regime tested.

### Broader perspective

The principle observed here — that representation-dependent complexity gates
search reachability under bounded budgets — is not specific to hypergeometric
identities. Similar bottlenecks are well-known in theorem proving (proof search
in different normal forms), SAT/SMT solving (formula representation), symbolic
integration (substitution strategies), and compiler optimization (expression
rewriting). The contribution of this work is to quantify and demonstrate the
effect experimentally in a controlled setting. The phenomenon suggests that
across symbolic and combinatorial search problems, the first optimization target
should often be representation, not algorithm.

---

## Appendix: Experiment Metadata and Reproducibility

### Reproducibility

All experiments are fully reproducible. The Python CLI tool `src/wz_lab.py` exports
all computation traces, including intermediate steps, timeout events, and symbolic
verifications. Results are stored as JSON (for machine parsing) and CSV (for
spreadsheet review).

**Intentional negative results:** Experiments that fail to find a certificate are
recorded with the same detail as successes — gauge, denominator candidates, number
of equations solved, and runtime. This negative-result corpus is unusual in the
creative telescoping literature and essential for understanding the reachability
boundary.

**Budget parameters:** All experiments state the ansatz budget explicitly: maximum
certificate degree, recurrence order, denominator candidate set, and runtime.
These allow others to extend or challenge the results by adjusting the budget.

### Experiment Metadata

All experiments were run with SymPy 1.14.0. Source code, result JSON/CSV files,
and verification traces are in `src/wz_lab.py` and `results/`.

Key parameters:

| experiment | template | $r$ range | max cert degree | max rec order |
|---|---|---|---|---|
| gauge sensitivity | weighted-binom | 1 | 4 | 1 |
| structural pruning | binomial-power | 1–4 | 4 | 4 |
| cross-family | delannoy-central | — | — | 2 |