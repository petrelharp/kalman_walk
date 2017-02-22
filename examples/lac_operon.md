Produces lactase, in the presence of lactose and the absence of glucose.

Input state variables are:

- lactose concentration $u_1(t)$
- glucose concentration $u_2(t)$

Inputs are concentration curves from E. coli eating up these things.

The first two internal state variables are concentrations of

- lactase $x_1(t)$, and 
- glucase $x_2(t).

These are also the output: $y_1=x_1$ and $y_2=x_2$.

Fitness function is total net amount of energy produced up until time $T$.
Suppose that one glucose is equal to $g$ lactoses with $g>1$.
The energy is
$$
E(T) = \int_0^T ( u_1(t) x_1(t) - \alpha_1 x_1(t)^2 ) dt
$$

BUT WHY THE QUADRATIC TERM?


