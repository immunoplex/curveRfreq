# Comparison of five nonlinear regression models for standard curve fitting

Consider the five nonlinear models used in the `curveRfreq` package with
concentration as the predictor. Writing $x$ for concentration:

1.  **logistic4: 4-parameter logistic in $x$**

$$y \sim a + \frac{d - a}{1 + \exp\left( \frac{x - c}{b} \right)}$$

2.  **log-logistic4: 4-parameter log-logistic (Hill-type)**

$$y \sim a + \frac{d - a}{1 + \left( \frac{x}{c} \right)^{b}}$$

3.  **gompertz4: 4-parameter Gompertz-type in $x$**

$$y \sim a + (d - a)\exp\left( - \exp\left( - b(x - c) \right) \right)$$

4.  **logistic5: 5-parameter generalized logistic (asymmetric)**

$$y \sim d + \frac{a - d}{\left( 1 + \exp\left( \frac{x - c}{b} \right) \right)^{g}}$$

5.  **loglogistic5: 5-parameter generalized log-logistic (asymmetric)**

$$y \sim a + (d - a)\left( 1 + g\exp\left( - b(x - c) \right) \right)^{- 1/g}$$

We compare their forms, how they fit assay response vs. concentration,
and their performance and numerical stability.

## 1. logistic4: 4-Parameter logistic model in concentration (x)

## Model form

### logistic4

$$y(x) = a + \frac{d - a}{1 + \exp\left( \frac{x - c}{b} \right)}$$

### Parameters

- $a$: lower asymptote (response at very low $x$)
- $d$: upper asymptote (response at very high $x$)
- $c$: inflection point in $x$ (the $x$ value where the response is
  halfway between $a$ and $d$)
- $b$: horizontal scale (slope) parameter on the linear $x$ scale

### Transform and symmetry

Rearranging:

$$\frac{y - a}{d - y} = \exp\left( \frac{x - c}{b} \right)$$

Taking logarithms:

$$\ln\left( \frac{y - a}{d - y} \right) = \frac{x - c}{b}$$

So this transformed response is linear in $x$, and the curve is:

- Symmetric (S-shaped) around $x = c$ on the original $x$ axis.
- $b$ controls how quickly the response transitions around $c$.

The response at the inflection point is:

$$y_{\text{inflection}} = \frac{a + d}{2}$$

### Fit behavior, performance, and stability

- Works well when the S-shape of $y$ vs. $x$ is roughly symmetric in
  $x$.
- Easy to interpret and typically numerically stable with reasonable
  starting values.
- May be less appropriate if concentration spans orders of magnitude and
  the natural symmetry is on a log scale, or if the empirical curve is
  clearly asymmetric.

## 2. log-logistic4: 4-parmameter log-logistic (Hill) model

## Model form

## loglogistic4

$$y(x) = a + \frac{d - a}{1 + \left( \frac{x}{c} \right)^{b}}$$

## Parameters

- $a$: lower asymptote
- $d$: upper asymptote
- $c$: inflection point (the $x$ value where $y$ is halfway between $a$
  and $d$)
- $b$: Hill slope parameter acting on a logarithmic $x$ scale##
  Transform and symmetry

Rearranging: $$\frac{d - y}{y - a} = \left( \frac{x}{c} \right)^{b}$$

Taking logarithms:

$$\ln\left( \frac{d - y}{y - a} \right) = b\left\lbrack \ln(x) - \ln(c) \right\rbrack$$

So the transformed response is linear in $\ln(x)$. The curve is:

- Symmetric when plotted against $\ln(x)$ (or $\log_{10}(x)$).
- $b$ controls how steeply the response changes per unit change in
  $\ln(x)$.
- The response at the inflection point
  is:$$y_{\text{inflection}} = \frac{a + d}{2}$$

## Fit behavior, performance, and stability

- Well suited for dose–response relationships where $x$ spans orders of
  magnitude and the curve looks symmetric on a log(concentration) axis.

- Parameters have clear interpretations:

  - $a$ and $d$ are asymptotes,  
  - $c$ is the EC50-type midpoint,
  - $b$ is the Hill slope.

- Numerically robust in typical dose–response applications with
  multiplicative concentration steps (for example, $0.1,1,10,100$).

- The intended log-scale interpretation breaks down if $x$ is near zero
  or negative, though the function can still be fitted algebraically
  when $x > 0$.

## 3. gompertz4: 4-parameter Gompertz-type model in x

## Model Form

### gompertz4

$$y(x) = a + (d - a)\exp\left( - \exp\left( - b(x - c) \right) \right)$$

## Parameters

- $a$: lower asymptote
- $d$: upper asymptote as $x$ becomes large
- $c$: location parameter associated with the inflection region
- $b$: growth/slope parameter on the linear $x$ scale

## Shape and asymmetry

Define:

$$G(x) = \exp\left( - \exp\left( - b(x - c) \right) \right)$$

Then:

$$y(x) = a + (d - a)G(x)$$

where:

- $G(x)$ ranges from $0$ to $1$.
- For very low $x$,

$$\exp\left( - b(x - c) \right)$$

is large, so:

$$G(x) \approx \exp\left( - \text{large} \right) \approx 0$$

- For very high $x$,

$$- b(x - c)$$

is very negative, giving:

$$\exp\left( - b(x - c) \right) \approx 0$$

and therefore:

$$G(x) \approx \exp(0) = 1$$

The response at the inflection point is:

$$y_{\text{inflection}} = a + (d - a)e^{- 1}$$

## Curve characteristics

The Gompertz curve is:

- Inherently asymmetric in $x$.
- Typically rises relatively quickly from the lower tail and then
  approaches the upper asymptote more gradually (or vice versa depending
  on sign conventions).

## Fit behavior, performance, and stability

- Useful when data clearly show asymmetric sigmoidal behavior that
  logistic or log-logistic forms do not capture well.
- Provides flexibility particularly in how the upper tail approaches
  $d$.
- Parameter interpretation (especially $c$ and $b$) is generally less
  intuitive than for logistic/Hill-type models.
- Numerically more delicate because of nested exponentials; very large
  values of

$$\left| b(x - c) \right|$$

can cause underflow or overflow in exponential terms and may require
constrained fitting or careful starting values.

## 4. logistic5: 5-parameter generalized logistic model in $x$ (asymmetric logistic)

## Model form

### logistic5

$$y(x) = d + \frac{a - d}{\left( 1 + \exp\left( \frac{x - c}{b} \right) \right)^{g}}$$

It is convenient to rewrite the model as:

$$y(x) = a + (d - a)\left( 1 + \exp\left( \frac{x - c}{b} \right) \right)^{- g}$$

## Parameters

- $a$: lower asymptote
- $d$: upper asymptote
- $c$: location parameter analogous to a center or inflection-like
  position in $x$
- $b$: scale/slope parameter on the linear $x$ scale
- $g$: asymmetry (shape) parameter

The inflection location is:

$$x_{\text{inflection}} = c + b\log(g)$$

The response at the inflection point is:

$$y_{\text{inflection}} = d + \frac{a - d}{(1 + g)^{g}}$$

When:

$$g = 1$$

the `logistic4` model reduces to the standard 4-parameter logistic model
(up to the $a$–$d$ sign convention). Thus, `logistic5` generalizes
`logistic4` by introducing an additional shape parameter $g$ that allows
asymmetric sigmoidal forms on the linear $x$ scale.

## Shape and asymmetry

- For $g \neq 1$, the left and right shoulders of the sigmoid differ in
  steepness and extent; the curve is no longer symmetric about $x = c$.
- Values of $g > 1$ or $g < 1$ stretch or compress one side of the curve
  relative to the other.
- As: $\left. g\rightarrow 1 \right.$ the model approaches the symmetric
  logistic-in-$x$ behavior.

## Fit behavior

- The `logistic5` moodel can capture asymmetric sigmoidal relationships
  in $x$ while retaining logistic-like tails and a logistic-like central
  transition region.
- This is useful when the symmetric `logistic4` model misfits one tail
  but otherwise captures the overall S-shaped trend.

## Performance and stability

- More flexible than `logistic4` because of the additional asymmetry
  parameter $g$, but also more complex to fit.
- Increased flexibility can improve empirical fit quality but may also
  introduce parameter trade-offs and stronger parameter correlation
  (particularly among $b$, $c$, and $g$).
- These dependencies can reduce numerical stability or lead to
  near-nonidentifiability in some datasets.
- Good starting values (often initializing with $g \approx 1$) and
  parameter bounds are often important for stable estimation and to
  avoid degenerate curve shapes.

## loglogistic5: 5-parameter generalized log-logistic model (asymmetric)

## Model form

### loglogistic5

$$y(x) = a + (d - a)\left( 1 + g\exp\left( - b(x - c) \right) \right)^{- 1/g}$$

## Parameters

- $a$: lower asymptote
- $d$: upper asymptote
- $c$: location parameter on the $x$ scale
- $b$: scale/slope parameter, similar to a logistic slope factor in $x$
- $g$: shape/asymmetry parameter

The inflection location is:

$$x_{\text{inflection}} = c + \frac{\log(g)}{b}$$

The response at the inflection point is:

$$y_{\text{inflection}} = a + (d - a)(1 + g)^{- 1/g}$$

For small values of $g$,

$$\left. g\rightarrow 0 \right.$$

the model approaches an exponential-type logistic form. For nonzero $g$,
the model defines a generalized family that includes several asymmetric
sigmoidal shapes and is closely related to generalized logistic families
and certain reparameterizations of Gompertz and log-logistic models.

## Shape and asymmetry

- The curve is sigmoidal between $a$ and $d$, while allowing flexible
  asymmetric transition behavior.
- The parameter $g$ explicitly modulates curvature and asymmetry of the
  S-curve.
- Depending on the sign and magnitude of $g$, one tail may become
  longer, flatter, or steeper relative to the other.

## Relation to log-scale interpretations

Although the model contains the exponential term:

$$\exp\left( - b(x - c) \right)$$

the resulting model is not the standard symmetric log-logistic form.
Instead, it represents a generalized asymmetric family capable of
approximating:

- Logistic-type curves,
- Gompertz-like curves,
- Other asymmetric sigmoidal shapes.

Thus, `loglogistic5` provides substantially greater flexibility than the
symmetric `loglogistic4` model.

## Fit behavior, performance, and stability

- Potentially the most flexible of the five models for describing
  asymmetric S-shaped dose–response relationships in $x$.
- The additional flexibility comes with increased risk of overfitting
  and stronger parameter correlations, particularly among $b$, $c$, and
  $g$, and in interaction with asymptotes $a$ and $d$.
- Numerical fitting may be less stable when data do not strongly
  constrain both tails and the central transition region.
- Good starting values and reasonable parameter bounds for $g$ are often
  important for stable optimization.
- Interpretation of $g$ is more abstract than the primary logistic
  parameters, since it primarily governs asymmetry and tail behavior
  rather than a directly observable biological quantity.

## 6. Comparative Summary of all five models

**(a) Symmetry, scale, and shape**

- `logistic4` (4-parameter logistic in x)

  - Symmetric S-curve on the **linear x scale**

  - Appropriate when the empirical curve is roughly symmetric in x

- `loglogistic4` (4-parameter log-logistic)

  - Symmetric S-curve on a $log(concentration)$ scale.

  - Appropriate for standard dose–response data spanning orders of
    magnitude, where plotting vs. $log(x)$ gives a symmetric S.

  - `gompertz4` (Gompertz-type in $x$)

    - Intrinsically asymmetric on the $x$ scale.
    - Good when one tail is clearly more gradual or more abrupt than the
      other.

- `logistic5` (5-parameter generalized logistic in $x$):

  - Logistic-like but with an additional shape parameter $g$ to control
    asymmetry in $x$.
  - Bridges between symmetric logistic behavior $(g = 1)$ and a family
    of asymmetric logistic-like curves.

- `loglogistic5` (5-parameter generalized log-logistic):

  - Generalized asymmetric sigmoidal form with rich shape control via g.
  - Can mimic or approximate various asymmetric logistic/Gompertz-type
    shapes.

**(b) Interpretation of slope and shape parameters**

- `logistic4`

  - $b$ is a horizontal slope/scale in units of $x$; $c$ is the
    inflection point.

- `loglogistic4`

  - $b$ is a Hill slope on the log scale; $c$ is the inflection point;
    shape is symmetric in $ln(x)$.

- `gompertz4`

  - $b$ controls growth and decay in a nested exponential; $c$ controls
    location of the rapid-growth region; no simple “symmetric slope”
    interpretation.

- `logistic5`

- $b$ still controls global steepness in $x$, but $g$ adjusts asymmetry
  of the S-curve around $c$.

- When $g \approx 1$), it behaves like `logistic4`; deviations from 1
  tilt/stretch the curve’s shoulders.

- `loglogistic5`

  - b affects how quickly the transition occurs as $x$ passes $c$; $g$
    fine-tunes the asymmetry and tail behavior of the curve built from
    $exp\left( - b*(x - c) \right)$.

- Interpretation is more qualitative (shape control) than a simple
  single-slope concept.

**(c) Fit behavior across wide concentration ranges**

- `logistic4`

  - Best when data are reasonably symmetric on a linear $x$ axis, and
    $x$ does not span extremely wide ranges.

- `loglogistic4`

  - Typically best for classical dose–response (e.g., pharmacology,
    toxicology) where $log_{10}(x)$ spacing is used; symmetry in
    log-space.

- `gompertz4`

  - Useful when the approach to the upper plateau (or the rise from the
    lower plateau) is clearly more stretched than the other tail and
    simpler logistic forms cannot capture this.

- `logistic5`

  - Adds asymmetry to the logistic-in-x family, improving fit where
    `logisic4` is too rigid but the Gompertz-type `gompertz4` does not
    match the central region well.

- Good when the overall S-shape is logistic-like but visibly skewed

- `loglogistic5`

- The most flexible of the “logistic / log-logistic-like” families
  considered here.

- Can capture nuanced asymmetries and tail behaviors over wide $x$
  ranges, particularly useful when both tails and the center are
  well-sampled and clearly deviate from standard logistic or Hill forms.

**(d) Performance and numerical stability**

- `logistic4`
  - Simple and generally stable; minimal risk of
    over-parameterization.  
  - May misfit asymmetric curves.
- `loglogistic4`
  - Similarly stable and widely used; excellent for typical log-spaced
    dose–response data.
  - Interpretation of parameters is straightforward and biologically
    meaningful in many contexts.
  - Use when concentrations span orders of magnitude and the curve is
    symmetric in $log(x)$; this is often the default for classical
    dose–response analysis.
- `gompertz4`
  - More flexible for asymmetry but can be numerically sensitive due to
    nested exponentials; large values of $\left| b*(x - c) \right|$ can
    cause extreme function values.
  - Requires careful starting values and potentially parameter
    constraints.
- `logistic5`
  - More flexible than `logistic4` but with an extra parameter $g$,
    which may be weakly identified if data do not strongly exhibit
    asymmetry.
  - Potential for increased correlation between parameters $(b,c,g)$,
    leading to less stable fits if data are sparse or noisy in either
    tail.
- `loglogistic5`
  - Highest flexibility and highest parameter count among the five.
  - Can deliver superior fits on complex asymmetric data, but has the
    greatest risk of:
    - Overfitting,
    - Poor identifiability of $g$,
    - Numerical instability without good initial values and appropriate
      bounds.
    - Interpretation focuses on qualitative shape rather than simple
      inflection point and slope.
  - Interpretation focuses on qualitative shape rather than simple
    inflection point and slope.
