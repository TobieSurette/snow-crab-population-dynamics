# Snow crab instar analysis

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
 - **multiple local solutions** : in that different solutions can give similar fits to the data. Some of these solutions are also degenerate, e.g. shrinking the variance around a particular
   observations yields arbitrarily high likelihood. 
 - **instar stage is unknown** : Other than a size-frequency analysis, there is no way to independently determine the instar stage. There is no internal structure that we can analyse and count 
   the number of benthic moults a crab has undergone. 
 - The **size overlap** between successive instars increases as they grow, leading to difficulties in resolving larger instars (generally starting at instar **IX** or larger).
   Size overlap also increases with maturition, leading to sill more variation in growth. Because of this, unstructured mixture models tend converge to converge to solutions which are
   at odds with biological knowledge or assumptions for larger instar sizes. 

## Methods 

### Data 
Data used are from the southern Gulf of St. Lawrence (sGSL) **snow crab survey**. This is a scientific survey that specially targets snow crab in the sGSL. It uses a European-style Nephrops trawl 
that is designed to dig into soft sediment. 

### Modelling approach
- Incorporate some **biological growth structure** into the mixture model.
- Make use of the fact that sometimes, instars are more **easily resolved at the local level** than by analyzing an aggregate sample (say over a large area).
-	Build a **hierarchical model** which is able to model mixture models both at the local and global levels, leading to growth and instar abundance inferences at both the local and global level.

## References
Marianne Alunno-Bruscia and Bernard Sainte-Marie. 1998. Abdomen allometry, ovary development, and growth of female snow crab, *Chionoecetes opilio* (Brachyura, Majidae), in the northwestern Gulf of St. Lawrence. Canadian Journal of Fisheries and Aquatic Sciences. 55(2): 459-477. https://doi.org/10.1139/f97-241.

J. M. (Lobo) Orensanz, Billy Ernst, David A. Armstrong. 2007. Variation of Female Size and Stage at Maturity in Snow Crab (*Chionoecetes opilio*) (Brachyura: Majidae) from the Eastern Bering Sea, Journal of Crustacean Biology, Volume 27, Issue 4(1), pp. 576–591, https://doi.org/10.1651/S-2790.1.

Bernard Sainte-Marie, Sylvain Raymond, and Jean-Claude Brêthes. 1995. Growth and maturation of the benthic stages of male snow crab, *Chionoecetes opilio* (Brachyura: Majidae). Canadian Journal of Fisheries and Aquatic Sciences. 52(5): 903-924. https://doi.org/10.1139/f95-091.

Weldon, W.F.R. 1893. On certain correlated variations in *Carcinus moenas*. Proceedings of the Royal Society of London 54. pp. 318-329.


