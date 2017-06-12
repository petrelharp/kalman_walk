Papers can be found [here](https://www.bibsonomy.org/user/peter.ralph/quantitative_genetics).

Generally, short-term response to selection is fairly robust;
the dependence of available genetic variation on mutation-selection balance (and associated G matrix) is less clear.

In our model, each *row* of $A$ is a single locus with relatively little recombination.
So, assuming that each entry of $A$ is an independent Gaussian, and the distribution is filled out,
is probably unreasonable.


- Environmental variance ($V_E$) can be incorporated into the fitness function,
    by replacing $V_S$ with $V_S+V_E$.

- The *infinitesimal model* (barton2016infinitesimal) has that an offspring's genetic value is
    (parent midpoint) plus independent Gaussian with (segregation) variance $V_0$.
    - from barton2016infinitesimal: 
        "[The infinitesimal model] should not be confused with two other models that have been used more extensively.
        Kimura (1965) investigated the distribution of effects of alleles at a single locus, and approximated
        this continuum-of-alleles model by a Gaussian; Lande (1976a) developed this model to investigate
        maintenance of variation by mutation, despite stabilizing selection. This is a quite different approach
        from the infinitesimal model, which requires no strong assumptions about the distribution
        of effects at each locus, and which does not assume a Gaussian distribution of trait values."

- lande1979quantitative does multidimensional things assuming always Gaussian:
    * change in population mean is $G \grad \log \bar W$, and $\bar W$ is the mean population fitness
    * in the absence of directional selection distribution of population mean after time $t$ has covariance matrix $tG/N_e$
    * and still without constraint the steady-state mutation-drift balance has $G = 2 U$, 
        where $U$ is the pleiotropic one-generation mutation covariance matrix

- "house of cards" means that mutants alleles have trait independent of parental trait

- It is not clear at all that the equilibrium level of diversity under mutation-selection balance
    is as in Lande.  turelli1984heritable ("Lerch's zeta meets the abdominal bristle")
    shows it depends on the modeling, and points out that in high dimensions there is too much trait space.

- barton2004effects looks at how genetic drift (during a bottleneck) affects genetic variance,
    including "converting" epistatic variance to additive variance

- turelli1994genetic (Turelli & Barton, "What, me normal?") show that the infinitesimal model does pretty good.
    "Our basic conclusion, which could not have been foreseen from our weak
    selection analysis, is that even intense truncation selection on an additive
    polygenic trait is not likely to produce significant departures from (1) caused
    by higher order disequilibria. " 

- Dominance makes things more complicated (Barton & Turelli 2004: barton2004effects).

- jones2003stability (Jones, Arnold, and Burger) look at stability of the G matrix with simulations.
    - and similarly, and jones2014epistasis look at how the mutational architecture G matrix) of traits evolves, with simulation
