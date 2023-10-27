# Size and instar-based population modelling

## Biological background

Female snow crab to 11 benthic instars, i.e. they moult up to 11 times before attaining sexual maturity at the terminal moult.
Benthic instars are numbered using roman numerals from **I** to **XI**. Each instar can be grouped into one of four stages, based on maturity and shell condition:

| Stage              | Description                                                                                            |
| ------------------ | :----------------------------------------------------------------------------------------------------- |
| immature           | Morphometrically **immature** with small **white**-coloured gonads.                                    |
| pubescent          | Morphometrically **immature** with larger **orange**-coloured gonads.                                  |
| new-shelled mature | Morphometrically **mature** with **new-shell** carapace condition (i.e. moulted in the current year).  |
| old-shelled mature | Morphometrically **mature** with **old-shell** carapace condition (i.e. moulted in past years).        |

The immature maturity stage lasts from instar **I** up to instar **IX**. Immature instars **VII**, **VIII** and **IX** can then transition to the pubescent maturity stages as instars 
**VIII**, **IX** and **X**, respectively. Pubsecent instars **VIII**, **IX** and **X** then transition to primiparous instars **IX**, **X** and **XI** during the following terminal moult, respectively. 
Mature females have attained full sexual maturity and longer moult. New-shelled mature females become old-shelled females the following year. 

Note : Primiparous females are new-shelled matures, but the converse is not necessarily true. Primiparous females are females coarrying their first clutch of eggs, but they do so for one or two years. Thus, a primiparous female in the second year of incubating its eggs would be identified as an old-shelled mature. All multiparous females are old-shelled matures.

The global maturity-stage can be represented as the following diagram:

<div align="center">
<img
  src="https://github.com/TobieSurette/snow-crab-population-dynamics/assets/14942142/b29991c4-1150-40c0-9166-f6cb24ac11f4"
  style="display: inline-block; margin: 0 auto"
  align="center" 
  width=150>
</div>

where arrows represent annual increments of time.

The global maturity-stage can be broken down into instars, as in the following diagram:

<div align="center">
<img
  src="https://github.com/TobieSurette/snow-crab-population-dynamics/assets/14942142/94519558-9c52-49b1-b799-6d971ab34438"
  style="display: inline-block; margin: 0 auto"
  align="center" 
  width=700>
</div>

## Modelling approach

Size-frequency distributions for female snow crab partitioned by year and maturity group were analyzed using a Gaussian mixture modelling framework. Females were chosen because they present a simpler target for analysis: 

1) females are smaller than males and so their instars are better resolved, i.e. larger instars overlap too much and become muddied by variations in growth and maturation.
2) pubescent females are easily identified through their orange gonads, pubescent males are not so easily identified.
3) skip-moulting among pubescents is negligible among females.
4) there is no direct impact of the fishery.

From the dynamics of females, we may glean some insight on population recruitment, survey catchability, maturation and natural mortality processes. 
These do have implcations on their male counterparts and will aid in interpreting percieved patterns for male snow crab categories.

For each maturity group, a Gaussian mixture model is defined :

$$\sum_{j=1}^k \pi_j \times \Phi \left( x | \mu, \sigma \right) $$


