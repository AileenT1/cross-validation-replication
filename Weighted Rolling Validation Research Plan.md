I am studying the toy example of rolling-validation model-selection consistency.

The original setup is:

- \(Y_1,\ldots,Y_n \overset{iid}{\sim} N(0,1)\)
- Model 0 uses the sequential sample mean:
  \[
  \hat\mu_{i-1}=\frac{1}{i-1}\sum_{j=1}^{i-1}Y_j
  \]
- Model 1 is the oracle predictor 0.
- The original rolling-validation statistic is
  \[
  \Xi_n
  =
  \sum_{i=2}^n
  \left[
  -2Y_i\hat\mu_{i-1}
  +
  \hat\mu_{i-1}^2
  \right].
  \]
- The incorrect-selection event is
  \[
  \Xi_n<0.
  \]
- In the unweighted case, the theory gives
  \[
  E[\Xi_n]\sim \log n,
  \qquad
  \operatorname{Var}(\Xi_n)=O(1),
  \]
  so
  \[
  P(\Xi_n<0)\to0.
  \]

My professor wants me to investigate whether this probability still goes to zero after introducing weights into the rolling-validation sum.

Please complete the project in the following order.

## Step 1: Reproduce the original unweighted simulation

Write clean modular R code that simulates the original statistic

\[
\Xi_n
=
\sum_{i=2}^n
\left[
-2Y_i\hat\mu_{i-1}
+
\hat\mu_{i-1}^2
\right].
\]

For each sample size \(n\), run many Monte Carlo replications and estimate

\[
P(\Xi_n<0).
\]

Use sample sizes such as

\[
n=100,200,500,1000,2000,5000,10000.
\]

Start with approximately 1000 Monte Carlo replications, but make this an adjustable parameter.

Verify that the estimated incorrect-selection probability decreases toward zero as \(n\) grows.

Also record:

- empirical mean of \(\Xi_n\)
- empirical variance of \(\Xi_n\)
- theoretical signal \(\log n\)
- empirical ratio
  \[
  \frac{\operatorname{Var}(\Xi_n)}
       {E[\Xi_n]^2}.
  \]

Do not move to weighted RV until the original simulation behaves correctly.

## Step 2: Introduce a general weighted statistic

Define

\[
\Xi_n^{(w)}
=
\sum_{i=2}^n
w_{i,n}
\left[
-2Y_i\hat\mu_{i-1}
+
\hat\mu_{i-1}^2
\right].
\]

Write the R functions so that the weight function can be passed as an argument.

Do not hard-code one specific weight.

The function should allow weights depending on \(i\), \(n\), or both.

## Step 3: Start with power weights

First investigate

\[
w_i=i^\alpha
\]

for

\[
\alpha=-1,-0.5,0,0.5,1.
\]

Remember that \(\alpha=0\) is the original unweighted rolling validation and should serve as the baseline.

For every combination of \(n\) and \(\alpha\), estimate:

1. \(P(\Xi_n^{(w)}<0)\)
2. \(E[\Xi_n^{(w)}]\)
3. \(\operatorname{Var}(\Xi_n^{(w)})\)
4. \(\operatorname{Var}(\Xi_n^{(w)})/[E(\Xi_n^{(w)})]^2\)

Store all results in one tidy data frame.

## Step 4: Derive the expected weighted signal

Show mathematically that

\[
E[\Xi_n^{(w)}]
=
\sum_{i=2}^n
\frac{w_{i,n}}{i-1}.
\]

Explain this from

\[
E[-2Y_i\hat\mu_{i-1}]=0
\]

and

\[
E[\hat\mu_{i-1}^2]
=
\frac{1}{i-1}.
\]

For \(w_i=i^\alpha\), determine the asymptotic order of

\[
\sum_{i=2}^n\frac{i^\alpha}{i-1}.
\]

In particular, compare:

- \(\alpha=-1\)
- \(\alpha=-0.5\)
- \(\alpha=0\)
- \(\alpha=0.5\)
- \(\alpha=1\)

Identify which choices produce a diverging expected signal and which produce only a bounded signal.

## Step 5: Investigate the variance

Do not assume that a diverging mean is enough for consistency.

For each weight, study how

\[
\operatorname{Var}(\Xi_n^{(w)})
\]

changes with \(n\).

Use simulation first.

Then try to derive or approximate its asymptotic order analytically.

The main quantity of interest is

\[
R_n
=
\frac{\operatorname{Var}(\Xi_n^{(w)})}
{\left(E[\Xi_n^{(w)}]\right)^2}.
\]

The goal is to determine whether

\[
R_n\to0.
\]

Explain that if this ratio goes to zero, Chebyshev's inequality gives

\[
P(\Xi_n^{(w)}<0)
\le
\frac{\operatorname{Var}(\Xi_n^{(w)})}
{\left(E[\Xi_n^{(w)}]\right)^2}
\to0.
\]

## Step 6: Make diagnostic plots

Create separate ggplot2 figures for:

1. estimated incorrect-selection probability versus \(n\)
2. empirical mean versus \(n\)
3. empirical variance versus \(n\)
4. variance divided by squared mean versus \(n\)

Use one line for each value of \(\alpha\).

A log scale for \(n\) is appropriate.

The most important plot is

\[
n
\quad\text{versus}\quad
P(\Xi_n^{(w)}<0).
\]

## Step 7: Check larger sample sizes if necessary

If some curves do not clearly stabilize, extend the sample sizes, for example to

\[
n=20000,50000,100000.
\]

Reduce the Monte Carlo replication count if runtime becomes excessive.

The code should be efficient enough that changing the \(n\)-grid and Monte Carlo count is easy.

## Step 8: Compare normalized weights

After studying \(w_i=i^\alpha\), also investigate normalized versions such as

\[
w_{i,n}
=
\left(\frac{i}{n}\right)^\alpha.
\]

Explain whether multiplying every weight by a common positive constant changes the sign of the RV statistic and therefore whether it changes the model-selection event.

Use this observation to avoid running redundant simulations.

## Step 9: Summarize the result

Create a final summary table with columns similar to:

- weight
- \(\alpha\)
- expected-signal order
- empirical variance behavior
- empirical error probability behavior
- whether simulations suggest \(P(\Xi_n^{(w)}<0)\to0\)

Do not claim a formal theorem unless the analytical argument actually proves it.

Clearly distinguish:

- simulation evidence
- mathematical derivation
- conjecture

## Coding requirements

Use R.

Keep the code modular. Ideally use functions such as:

- `simulate_weighted_rv_once()`
- `estimate_weighted_rv()`
- `run_weight_grid()`

Avoid unnecessary packages beyond tidyverse/ggplot2 if possible.

Set a random seed.

Do not overwrite results unnecessarily.

Make Monte Carlo repetitions, sample sizes, and weight functions easy to modify.

Return both the R code and a short mathematical explanation of the results.