# Analytical standard errors: oracle recursive IF

The purpose of this note is to explain where the recursive influence function (IF) comes from, for SI and IV. We work at the population truth; implementation simply substitutes estimated quantities and updates at event times using left limits.

Fix one cause and write $`P_0`$ for the true distribution. Suppose a cumulative effect $`B`$ solves

```math
dB(t)=\Psi(P,B,\theta;dt).
```

Perturb both the distribution and the nuisance parameters in subject $`i`$'s direction:
$`P_a=(1-a)P_0+a\delta_{O_i}`$ and $`\theta_a=\theta_0+a\phi_i^\theta`$.
Differentiating the estimating equation gives

```math
\boxed{
 d\psi_i(t)
 =D_P\Psi[\delta_{O_i}-P_0]
 +D_B\Psi[\psi_i]
 +D_\theta\Psi[\phi_i^\theta].
}
```

The three terms are the direct contribution, feedback from earlier effect estimates, and nuisance-estimation contribution. All derivatives below are evaluated at the truth.

We write $`\check Y^1(t)`$ for the stage-eligible Fine--Gray risk weight. For censoring nuisance, let $`Q(t;\phi_i^C)`$ satisfy

```math
\dot{\check Y}^{1}(t;\phi_i^C)=-\check Y^1(t)Q(t;\phi_i^C).
```

Under the Cox censoring model, this is the directional derivative of the censoring cumulative-hazard difference; it is zero unless the generic subject had a competing event before $`t`$. In particular, $`\dot{\check Y}^{1}\,dN^1=0`$.

For subjects eligible at stage start $`t_k`$, $`D(t)=D_k`$ throughout $`[t_k,t_{k+1})`$, even if an event or censoring occurs inside the interval. Treatments in later stages are undefined after exit, not set to zero. All stage-$`k`$ expressions involving $`D_k`$ are restricted to stage-eligible subjects.

## Part I. Sequential Ignorability (SI)

Fix stage $`k`$ and let $`A_{k0}(t)`$ be its true cumulative effect. Define

```math
\begin{aligned}
m(t)&=\mathrm{expit}\{X_k^\top\alpha_{k0}-A_{k0}(t^-)\},
& q(t)&=D_k-m(t),\\
v(t)&=m(t)\{1-m(t)\},
& C_k(t)&=E_0[q(t)\check Y^1(t)D_k],\\
dR(t)&=\check Y^1(t)\{dN^1(t)-D_k\,dA_{k0}(t)\}.
\end{aligned}
```

The estimator solves the ratio equation

```math
dA_k(t)=\frac{E[q(t)\check Y^1(t)dN^1(t)]}
{E[q(t)\check Y^1(t)D_k]}.
```

The key derivatives are $`\partial q/\partial A_k=v`$ and
$`\partial q/\partial\alpha_k=-vX_k`$. Applying the quotient rule, together with the censoring-weight derivative above, yields the full oracle IF recursion:

```math
\boxed{
\begin{aligned}
d\psi_{ik}^A(t)={}&\frac{q_i(t)}{C_k(t)}\,dR_i(t)
+\frac{E_0[v(t)dR(t)]}{C_k(t)}\,\psi_{ik}^A(t^-)\\
&-\frac{E_0[v(t)X_k^\top\phi_i^{\alpha_k}\,dR(t)]}{C_k(t)}
+\frac{E_0[q(t)\check Y^1(t)D_k Q(t;\phi_i^C)]}{C_k(t)}\,dA_{k0}(t),
\end{aligned}
\qquad \psi_{ik}^A(t_k^-)=0.
}
```

The first term is direct, the second is historical feedback, and the last two propagate estimation of treatment and censoring nuisance parameters. For $`t\in[t_k,t_{k+1})`$, the cumulative-effect IF is

```math
\boxed{\psi_i^{B_D}(t)=\sum_{\ell\lt k}\psi_{i\ell}^A(t_{\ell+1}^-)+\psi_{ik}^A(t).}
```

## Part II. Instrumental Variables (IV)

Let $`B_0(t)`$ be the true cumulative effect, and at $`t\in[t_k,t_{k+1})`$ define

```math
\begin{aligned}
J(t)&=\int_{[0,t)}D(u)\,dB_0(u),
& W(t)&=(Z-X^\top\alpha_0)e^{J(t)},\\
C(t)&=E_0[W(t)\check Y^1(t)D_k],
& dR(t)&=dN^1(t)-\check Y^1(t)D_k\,dB_0(t).
\end{aligned}
```

Again, the estimator is defined by a ratio:

```math
dB(t)=\frac{E[W(t)dN^1(t)]}{E[W(t)\check Y^1(t)D_k]}.
```

The source of the recursion is the derivative of the exponential weight:

```math
\delta W(t)
=-e^{J(t)}X^\top\delta\alpha
+W(t)\int_{[0,t)}D(u)\,d(\delta B)(u).
```

For a subject eligible at stage $`k`$, every $`D_1,\ldots,D_k`$ is defined and constant within its stage. Splitting the history by stages gives

```math
\int_{[0,t)}D(u)\,d\psi_i^{B_D}(u)
=\sum_{\ell=1}^k D_\ell\,\psi_{i\ell}(t^-),
```

where $`\psi_{i\ell}`$ is the IF contribution accumulated during stage $`\ell`$ (previous-stage contributions remain fixed). Thus the quotient rule gives

```math
\boxed{
\begin{aligned}
d\psi_{ik}(t)={}&\frac{W_i(t)}{C(t)}\,dR_i(t)
+\sum_{\ell=1}^k\frac{E_0[W(t)D_\ell\,dR(t)]}{C(t)}\,\psi_{i\ell}(t^-)\\
&-\frac{E_0[e^{J(t)}X^\top\phi_i^\alpha\,dR(t)]}{C(t)}
+\frac{E_0[W(t)\check Y^1(t)D_k Q(t;\phi_i^C)]}{C(t)}\,dB_0(t),
\end{aligned}
\qquad \psi_{ik}(t_k^-)=0.
}
```

The four contributions are direct, feedback from all earlier and current treatment stages, instrument nuisance, and censoring nuisance. At time $`t`$, only the current-stage contribution $`\psi_{ik}`$ is updated. The complete IF is

```math
\boxed{\psi_i^{B_D}(t)=\sum_{\ell=1}^k\psi_{i\ell}(t),\qquad t\in[t_k,t_{k+1}).}
```

## Standard error

The recursion is evaluated by replacing oracle quantities with estimates at each distinct cause-1 event time, always using the pre-update IF on the right-hand side. After combining all stage and nuisance contributions within each subject, use the uncentered variance convention

```math
\boxed{
\widehat{\mathrm{Var}}\{\widehat B_D(t)\}
=\frac1{n^2}\sum_{i=1}^n\{\widehat\psi_i^{B_D}(t)\}^2,
\qquad \widehat{\mathrm{SE}}=\sqrt{\widehat{\mathrm{Var}}}.
}
```
