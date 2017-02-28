In the presence of sugar,
the organism needs to start producing something to digest sugar.

Input state variables are:

- $u_1(t)$ : concentration of sugar.

Input functions are:

- oscillating?
- step function, on/off?

# Option 1: internal states are only TFs.

First two internal states are:

- $x_1(t)$ concentration of the thing that binds to the promoter of sugar-digesting enzyme
- $x_2(t)$ concentration of toxin

Output is the same as these two, $y_1=x_1$ and $y_2 =x_2$.

Suppose that digesting sugar produces one unit of energy but also $\gamma$ units of toxin,
which degrades at rate $\lambda$: so these entries of $A$ are fixed.

With small rate $\delta$, the enzyme and the toxin combine to produce something fatal,
so it is killed with rate equal to $\delta x_1(t) x_2(t)$,
so the probability of survival is
$$
\exp\left( - \int_0^\infty \delta x_1(t) x_2(t) dt \right) .
$$

