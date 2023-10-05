# Snow crab instar analsyis

## Background:

Snow crab size distributions, especially for smaller sizes, typically exhibit characteristic peaks clustered around the mean sizes of succeessive instars (i.e. moult stages). The properties of
these instars can be derived from analysis of size-frequency data, leading to information about their growth characterisitcs and their abundances. Snow crab instars are numbered using roman
numericals, starting with the first benthic stage after larvae have settled on the bottom and moulted, referred to as **instar I**. Sexual maturation starts at **instar VIII**. Larger instars can 
undergo a final moult to maturity (called the terminal moult) whereupon they attain full sexual maturity.

A common approach is to apply **finite mixture model**, comprised of a fixed number of Gaussian components whose statistical properties, namely their mean sizes, standard errors and proportions 
are estimated from fitting the model to a set of observed size data.

However, there are issues with mixture models:
 - **identifiability** : a mixture fit to the data is invariant (i.e. does not change) when we exchange the component labels. So some ad hoc constraint is usually applied, such as
   forcing the component means to be increasing, needs to be applied to properly constrain the parameter space.
 - **multiple local solutions** : in that different solutions can give similar fits to the data. In addition, some of these are degenerate, e.g. shrinking the variance around a particular
   observations yields arbitrarily high likelihood. 
 - Analytical solutions are not available, iterative methods must be used.
 - Other than a size-frequency analysis, there is no way to independently determine the instar stage, as there is no internal structure that we can i.e. there is no structure that we can
   analyse and count the number of moults a crab has undergone.
 - The **size overlap** between successive instars increases as they grow, leading to difficulties in resolving larger instars (generally starting at instar **IX** or larger).
   Size overlap also increases with maturition, leading to sill more variation in growth. Because of this, unstructured mixture models tend converge to converge to solutions which are
   at odds with biological knowledge or assumptions for larger instar sizes. 

## Methods 

### Data 
Data used are from the southern Gulf of St. Lawrence (sGSL) **snow crab survey**. This is a scientific survey that specially targets snow crab in the sGSL. It uses a European-style Nephrops trawl 
that is designed to dig into soft sediment. 

### Modelling approach
- However, while larger, maturing instars are harder to separate, they are often more clearly resolved in specific areas in specific years, i.e. global samples are aggregated of growth variations
  which have occurred over space and time, but specific samples are subject to less growth variation. Although local sample sizes are generally not sufficient to perform an independent analysis.
-  Incorporate some structure  the mixture model which tailored for analyzing snow crab on a single sample, relatively fast.
•	Extend the structure to account for hierarchical samples, applicable for a generating space and time-specific (i.e. inferences). 

## References
•	Weldon, W.F.R. 1893. On certain correlated variations in Carcinus moenas. Proceedings of the Royal Society of London 54. Pp. 318-329.



