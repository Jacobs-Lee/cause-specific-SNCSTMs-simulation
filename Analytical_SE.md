# Analytical standard errors

This note describes the **recursive influence-function (IF) calculations** for
one target cause. It is organized into two parts: **SI** and **IV**. All subject-level
IFs are on the full-sample scale, i.e.,
$\widehat\vartheta-\vartheta_0\approx n^{-1}\sum_{i=1}^n\psi_i$;
they have **not** been divided by $n$.

## Part I. Sequential Ignorability (SI)

### 1. Initialize each stage and estimate its increments

For stage $k$, let $s_1<\cdots<s_M$ be the distinct target-event times in
$[t_k,t_{k+1})$. Initialize the stage effect and its $n$ subject-level IFs as
$\widehat A_k(t_k-)=0$ and $\psi_{ik}^A(t_k-)=0$.

At time $s_j$, write $D_r=D_{rk}$, $\Delta N_{rj}=\Delta N_r^1(s_j)$, and
$a_{rj}=\widehat{\check Y}_r^1(s_j)$ (`risk_weight`). This weight is zero if
subject $r$ is not eligible at $t_k$, one if eligible and still observed at
$s_j$, and $\widehat G(s_j\mid L_r)/\widehat G(T_r\mid L_r)$ if eligible and
an observed competing event occurred before $s_j$; it is zero otherwise.
All sums below use the **full sample size** $n$ as denominator.

Using the stage-specific fitted logistic treatment model, compute

```math
\begin{aligned}
m_{rj}&=\mathrm{expit}\{X_{rk}^{\mathsf T}\widehat\alpha_k-\widehat A_k(s_j-)\},
& q_{rj}&=D_r-m_{rj},
& v_{rj}&=m_{rj}(1-m_{rj}),\\
C_j&=\frac1n\sum_r q_{rj}a_{rj}D_r,
& b_j:=\Delta\widehat A_k(s_j)&=\frac{\sum_r q_{rj}a_{rj}\Delta N_{rj}}{nC_j},
& r_{rj}&=\Delta N_{rj}-D_rb_j.
\end{aligned}
```

### 2. Calculate the direct IF, feedback, and nuisance correction

The **direct contribution** and **feedback coefficient** at $s_j$ are

```math
\boxed{
 h_{ij}=\frac{q_{ij}a_{ij}r_{ij}}{C_j},
 \qquad
 \kappa_j=\frac{1}{nC_j}\sum_r v_{rj}a_{rj}r_{rj}.
}
```

The feedback has a **positive sign**, since
$\partial(D_r-m_{rj})/\partial A_k=m_{rj}(1-m_{rj})=v_{rj}$.

Let $\phi_{ik}^{\alpha}$ be subject $i$'s full-sample-scale IF for
$\widehat\alpha_k$. The direct contribution from estimating the treatment
model is

```math
\boxed{
 f_{ij}^{\alpha}
 =-(\phi_{ik}^{\alpha})^{\mathsf T}
 \frac{\sum_r X_{rk}v_{rj}a_{rj}r_{rj}}{nC_j}.
}
```

For censoring, let $\phi_i^\gamma$ and $\phi_i^{\Lambda_C}(\cdot)$ be the IFs
for the Cox coefficient and baseline cumulative hazard (the latter already
includes the effect of estimating $\gamma$). For a subject $r$ with a
**competing event before** $s_j$, define

```math
\begin{aligned}
Q_{rij}=e^{X_{C,r}^{\mathsf T}\widehat\gamma}[&
\{\widehat\Lambda_C(s_j)-\widehat\Lambda_C(T_r)\}
 X_{C,r}^{\mathsf T}\phi_i^\gamma\\
&+\phi_i^{\Lambda_C}(s_j)-\phi_i^{\Lambda_C}(T_r)].
\end{aligned}
```

Set $Q_{rij}=0$ for all other subjects. The **censoring-nuisance contribution** is

```math
\boxed{
 f_{ij}^C=\frac{b_j}{nC_j}\sum_r q_{rj}a_{rj}D_rQ_{rij}.
}
```

This is the `risk_if` correction in the code: for prior competing events,
$\dot a_{rj}[\phi_i^C]=-a_{rj}Q_{rij}$, while
$\Delta N_{rj}=0$. Consequently, differentiating the weighted ratio gives
$+b_j\sum_r q_{rj}a_{rj}D_rQ_{rij}/(nC_j)$.

### 3. Update the IF recursively

At each event time, **compute everything above using $\widehat A_k(s_j-)$ and
$\psi_{ik}^A(s_j-)$, then update**

```math
\boxed{
\begin{aligned}
\psi_{ik}^A(s_j)
  &=(1+\kappa_j)\psi_{ik}^A(s_j-)
    +h_{ij}+f_{ij}^{\alpha}+f_{ij}^C,\\
\widehat A_k(s_j)&=\widehat A_k(s_j-)+b_j.
\end{aligned}
}
```

Repeat through the stage's event times; between event times the IF and effect
remain constant. **Tied target events are processed together**, not sequentially.
The implementation stores the fixed-nuisance and nuisance IFs as `E` and `U`:
$E_i\leftarrow(1+\kappa_j)E_i+h_{ij}$ and
$U_i\leftarrow(1+\kappa_j)U_i+f_{ij}^{\alpha}+f_{ij}^C$.
Their sum is exactly $\psi_{ik}^A=E_i+U_i$ above.

### 4. Combine stages and calculate standard errors

At the end of stage $k$, accumulate the stage estimates and **subject-level** IFs:

```math
\boxed{
\begin{aligned}
\widehat B_D(t_{k+1}-)&=\sum_{\ell=1}^k\widehat A_\ell(t_{\ell+1}-),\\
\psi_i^{B_D}(t_{k+1}-)&=\sum_{\ell=1}^k\psi_{i\ell}^A(t_{\ell+1}-),\\
\widehat{\mathrm{Var}}\{\widehat B_D(t_{k+1}-)\}
 &=\frac1{n^2}\sum_{i=1}^n(\psi_i^{B_D}(t_{k+1}-))^2.
\end{aligned}
}
```

**Add IFs before squaring**, not stage variances or standard errors. In the code,
these are `epsilon_A` (stage IF) and `epsilon_B` (cumulative IF).
The same **uncentered** convention, $\sum_i\psi_i^2/n^2$, applies to every SE.

For the pooled constant-effect coefficient, let $\rho_k(s)$ be the implemented
**unweighted** projection weight (the full-sample fraction eligible at $t_k$
and in the cause-specific subdistribution risk set at $s$), and let
$Q=\sum_k\int_{t_k}^{t_{k+1}}\rho_k(s)\,ds$ (approximated in the code by
midpoint integration). Then

```math
\widehat\beta
 =\frac{\sum_{k,j}\rho_k(s_j)b_{kj}}{Q},
\qquad
\psi_i^\beta
 =\frac{\sum_{k,j}\rho_k(s_j)
   \{\psi_{ik}^A(s_j)-\psi_{ik}^A(s_j-)\}}{Q},
\qquad
\widehat{\mathrm{Var}}(\widehat\beta)
 =\frac1{n^2}\sum_i(\psi_i^\beta)^2.
```

**Convention for pooled $\beta$:** the current calculation holds the empirical
$\rho_k$ and $Q$ fixed when computing its IF. If $\beta$ is interpreted as a
projection under a nonconstant true effect, the IF of these estimated weights
would also be needed. Under an exact constant-effect model, that extra
first-order contribution cancels.

## Part II. Instrumental Variables (IV)

### 1. Estimate the effect increment at each target-event time

Let $s_1<\cdots<s_M$ be all distinct target-event times, and let $k=k(j)$
denote the stage containing $s_j$. Write $\mathcal R_k$ for subjects
**eligible at $t_k$**. For $r\in\mathcal R_k$, the treatment is constant,
$D_r(t)=D_{rk}$ on $[t_k,t_{k+1})$, even if an event or censoring occurs
during that stage. If this happens, $D_{r,k+1}$ and later treatments are
**undefined**, not set to zero. All treatment-dependent sums at $s_j$ are
therefore restricted to $r\in\mathcal R_k$, for whom $D_{r1},\ldots,D_{rk}$
are defined. The normalization remains the **full sample size** $n$.

Let $a_{rj}=\widehat{\check Y}_r^1(s_j)$ be the Fine--Gray `risk_weight`
defined in Part I, and $\Delta N_{rj}=\Delta N_r^1(s_j)$. With the fitted
linear instrument model $E(Z\mid L;\alpha)=X^{\mathsf T}\alpha$, compute

```math
\begin{aligned}
 z_r^c&=Z_r-X_r^{\mathsf T}\widehat\alpha,
 &\eta_{rj}&=\sum_{h<j}D_{r,k(h)}b_h,
 &w_{rj}&=z_r^c e^{\eta_{rj}},\\
 C_j&=\frac1n\sum_{r\in\mathcal R_k}w_{rj}a_{rj}D_{rk},
 &b_j:=\Delta\widehat B_D(s_j)
 &=\frac{\sum_{r\in\mathcal R_k}w_{rj}\Delta N_{rj}}{nC_j},\\
 \Delta\widehat R_{rj}
 &=\Delta N_{rj}-a_{rj}D_{rk}b_j.
\end{aligned}
```

Here $\eta_{rj}$ uses **only earlier effect increments** $b_h$; ties at
$s_j$ are handled together. At a target event, $a_{rj}=1$, so the numerator
could equivalently include $a_{rj}$. **Do not** simplify $a_{rj}D_{rk}$
to $D_{rk}$: subjects with earlier competing events may have nonunit weights.

### 2. Direct IF and historical feedback

The fixed-nuisance **direct contribution** is

```math
\boxed{
 h_{ij}=\frac{w_{ij}\Delta\widehat R_{ij}}{C_j}
 \quad (i\in\mathcal R_k),
 \qquad h_{ij}=0\quad(i\notin\mathcal R_k).
}
```

For each previous or current stage $\ell\le k$, define the feedback
coefficient

```math
\boxed{
 g_{\ell j}=\frac1{nC_j}
 \sum_{r\in\mathcal R_k}
 w_{rj}\Delta\widehat R_{rj}D_{r\ell}.
}
```

Indeed, perturbing an earlier effect increment changes the current
exponential weight by $\delta w_{rj}=w_{rj}\sum_{h<j}D_{r,k(h)}\delta b_h$.
Thus the feedback into the current increment is
$\sum_{\ell\le k}g_{\ell j}\psi_{i\ell}(s_j-)$, using the **previous**
stage-specific IFs. The $g_{\ell j}$ are not additionally multiplied by
$b_j$.

### 3. Direct nuisance contributions

Let $\phi_i^\alpha$ be the full-sample-scale IF for the instrument regression,
and $\phi_i^\gamma$, $\phi_i^{\Lambda_C}(\cdot)$ the IFs for the Cox
censoring coefficients and baseline cumulative hazard, respectively.
The baseline-hazard IF already includes estimation of $\gamma$.

The **instrument-model contribution** is

```math
\boxed{
 f_{ij}^\alpha
 =-\frac{(\phi_i^\alpha)^{\mathsf T}}{nC_j}
 \sum_{r\in\mathcal R_k}X_r e^{\eta_{rj}}\Delta\widehat R_{rj}.
}
```

For a subject $r\in\mathcal R_k$ whose **competing event precedes** $s_j$,
write

```math
\begin{aligned}
 Q_{rij}=e^{X_{C,r}^{\mathsf T}\widehat\gamma}[&
 \{\widehat\Lambda_C(s_j)-\widehat\Lambda_C(T_r)\}
 X_{C,r}^{\mathsf T}\phi_i^\gamma\\
 &+\phi_i^{\Lambda_C}(s_j)-\phi_i^{\Lambda_C}(T_r)].
\end{aligned}
```

Set $Q_{rij}=0$ otherwise. Since the derivative of the competing-event
risk weight is $\dot a_{rj}[\phi_i^C]=-a_{rj}Q_{rij}$ and that subject
has $\Delta N_{rj}=0$, the **censoring contribution** is

```math
\boxed{
 f_{ij}^{C}=\frac{b_j}{nC_j}
 \sum_{r\in\mathcal R_k}w_{rj}a_{rj}D_{rk}Q_{rij}.
}
```

Together, $f_{ij}^{\theta}=f_{ij}^\alpha+f_{ij}^C$ is the full nuisance
`forcing` term (`if_alpha`, `if_gamma`, `if_Lambda`, and `risk_if` in the code).

### 4. Update the IF recursively

Initialize the $n\times K$ matrices $E=U=0$ **once**, retaining prior-stage
columns as time advances. For every event time $s_j$ in stage $k$, evaluate
all terms using the pre-update values $E_{i\ell}(s_j-)$ and
$U_{i\ell}(s_j-)$, then update **only column $k$**:

```math
\boxed{
\begin{aligned}
 \Delta E_{ik}(s_j)
 &=h_{ij}+\sum_{\ell=1}^k g_{\ell j}E_{i\ell}(s_j-),\\
 \Delta U_{ik}(s_j)
 &=f_{ij}^\alpha+f_{ij}^C
   +\sum_{\ell=1}^k g_{\ell j}U_{i\ell}(s_j-),\\
 E_{ik}(s_j)&=E_{ik}(s_j-)+\Delta E_{ik}(s_j),\\
 U_{ik}(s_j)&=U_{ik}(s_j-)+\Delta U_{ik}(s_j).
\end{aligned}
}
```

Equivalently, for $\psi_{i\ell}=E_{i\ell}+U_{i\ell}$,

```math
\boxed{
 \Delta\psi_{ik}(s_j)
 =h_{ij}+f_{ij}^\theta
  +\sum_{\ell=1}^k g_{\ell j}\psi_{i\ell}(s_j-).
}
```

Columns $\ell\ne k$ do not change at $s_j$. **The subject index $i$
still ranges over all $n$**, including individuals no longer stage-eligible:
although their direct $h_{ij}$ is zero, their previously accumulated IF
or nuisance IF can contribute through feedback.

### 5. Cumulative IF, standard errors, and pooled beta

For $t\in[t_k,t_{k+1})$, the effect and its **full** subject-level IF are

```math
\widehat B_D(t)=\sum_{j:s_j\le t}b_j,
\qquad
\psi_i^{B_D}(t)=\sum_{\ell=1}^k\{E_{i\ell}(t)+U_{i\ell}(t)\}.
```

As in Part I, combine all stage columns **before squaring** and use the
uncentered variance convention:

```math
\boxed{
\widehat{\mathrm{Var}}\{\widehat B_D(t)\}
 =\frac1{n^2}\sum_{i=1}^n\{\psi_i^{B_D}(t)\}^2,
\qquad \widehat{\mathrm{SE}}=\sqrt{\widehat{\mathrm{Var}}}.
}
```

For the pooled constant-effect coefficient, use the same projection weights
$\rho_k(s)$ and denominator $Q$ defined in Part I:

```math
\begin{aligned}
 \widehat\beta
 &=\frac{\sum_j\rho_{k(j)}(s_j)b_j}{Q},\\
 \psi_i^\beta
 &=\frac{\sum_j\rho_{k(j)}(s_j)
    \{\Delta E_{i,k(j)}(s_j)+\Delta U_{i,k(j)}(s_j)\}}{Q},\\
 \widehat{\mathrm{Var}}(\widehat\beta)
 &=\frac1{n^2}\sum_i(\psi_i^\beta)^2.
\end{aligned}
```

This pooled-$\beta$ IF follows the **same fixed-projection-weight convention**
as SI: no additional derivative of $\rho_k$ or $Q$ is included.
