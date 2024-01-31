
<!-- README.md is generated from README.Rmd. Please edit that file -->

# narwind: Offshore wind impacts on North Atlantic right whales

<!-- badges: start -->

[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg?style=flat-square)](https://www.tidyverse.org/lifecycle/#maturing)
[![status](https://img.shields.io/badge/repo%20status-active-green.svg?style=flat-square)](https://www.repostatus.org/#active)

<!-- badges: end -->

The `narwind` R package provides methods to forecast the abundance of
critically endangered North Atlantic right whales (*Eubalaena
glacialis*) in the context of user-defined offshore wind development
scenarios. `narwind` implements a spatially-explicit bioenergetic PCoMS
model (*sensu* Pirotta et al. (2018)), whereby the movements of
different NARW cohorts (juveniles, adult males, pregnant females,
resting females, lactating females + dependent calves) are simulated
throughout a full calendar year, and population size projections are
made over a time horizon relevant to management (25–50 years). The model
operates on a daily scale and accounts for the effects of multiple
anthropogenic stressors known/presumed to affect right whale health,
reproduction, and survival, including: (1) direct mortality from vessel
strikes, (2) behavioral responses to noise exposure leading to
short-term cessation of foraging activities, and (3) increased energy
expenditure resulting from entanglement in fishing gear.

## Getting started

If you are just getting started with `narwind`, we recommend reading the
[tutorial
vignette](https://offshore-wind.github.io/narwind/articles/narwind.html),
which provides a detailed introduction to the package.

## Installation

Install the GitHub development version to access the latest features and
patches.

``` r
# install.packages("remotes")
remotes::install_github("offshore-wind/narwind") # OR

# install.packages("devtools")
devtools::install_github("offshore-wind/narwind")
```

Note that the package relies on compiled code (C++) and functionality
provided by the
[Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) package.
The `Rtools` software may be needed on Windows machines.

Installation instructions can be found at
<https://cran.r-project.org/bin/windows/Rtools/rtools40.html>

## Background

Marine renewable energy sources are poised to make a vital contribution
in the fight against climate warming and global anthropogenic change.
While offshore wind power remains a nascent market in the United States
(DeCastro et al. 2019), strong and persistent oceanic wind regimes in
the North Atlantic offer vast potential for the expansion of the
industry along the U.S. eastern seaboard. However, concerns have been
raised that the deployment of thousands of turbines may adversely affect
iconic marine wildlife such as seabirds and cetaceans (Verfuss et al.
2016; Goodale and Milman 2019). Much of what is known about the impacts
of offshore wind developments on cetaceans comes from studies conducted
in Europe, where efforts to monitor animal behavior during both
construction and operation activities have been made for a limited
number of species common to that area, primarily harbor porpoises
(Bailey, Brookes, and Thompson 2014; Thompson et al. 2010). Little
information exists on the responses of baleen whales to offshore wind
activities, and in contrast to European case studies on resident
cetacean populations, many marine mammals occurring within U.S. waters
are migratory. This makes it difficult to directly apply lessons learnt
and determine with certainty what potential offshore wind effects may be
emergent in the North Atlantic in the future (Petruny, Wright, and Smith
2014). It is also critical to recognize that every new development
represents an incremental addition to the pre-existing and rapidly
expanding footprint of human activities on offshore marine ecosystems.
There is therefore a need to consider the cumulative impacts of offshore
wind in the context of additional concurrent human uses of ocean areas.
This is particularly important for sensitive wildlife species, such as
those in low abundance and/or in decline, and those occurring within
narrow habitat ranges.

One such species is the Critically Endangered North Atlantic right whale
(NARW, *Eubalaena glacialis*), which numbers only 340 individuals and
continues to experience a downward trajectory towards extinction (Linden
2023; Richard M. Pace, Corkeron, and Kraus 2017). Assessing the full
scope of potential offshore wind impacts on NARWs across the U.S.
Atlantic calls for a quantitative framework capable of integrating the
state of individuals (e.g., energy reserves, reproductive status) with
the state of the surrounding biophysical and anthropogenic environment
(e.g., resource density, wind farm characteristics) to project changes
in NARW vital rates and translate short-term patterns of disturbance
into long-term demographic outcomes. Bioenergetic models represent such
an approach and have been used as a major implementation of the
Population Consequences of Disturbance (PCoD) paradigm (Costa 2012).
PCoD is a conceptual model of the pathways through which the effects of
sub-lethal stressors may scale up to population-level impacts (Pirotta
et al. 2018); it has been implemented in various ways and has become a
major component of several regulatory frameworks (e.g., U.S. Marine
Mammal Protection Act,16 U.S. Code § 1361 et seq.; Endangered Species
Act,16 U.S. Code § 1531 et seq.). PCoD models have been developed for
several cetacean species such as gray whales, humpback whales,
bottlenose dolphins, beaked whales, as well as pinnipeds like the
northern elephant seal. Significant attention has also been given to
NARW in the last decade (Rolland et al. 2016; Schick et al. 2013), but
without explicit consideration of offshore wind energy areas, or of the
interactions that may arise between the myriad of pressures which NARW
are currently facing. An expanded framework known as the Population
Consequences of Multiple Stressors (PCoMS) (Pirotta et al. 2022) now
exists to accommodate scenarios where individuals within a population
are exposed to multiple stressors that may combine in unexpected ways to
generate cumulative impacts greater than the individual sum of their
parts (Crain, Kroeker, and Halpern 2008; Pirotta et al. 2019). PCoMS
models encompassing each key NARW habitat have been identified as a high
priority target for management in a recent review of NARW health and
monitoring needs (Moore et al. 2021).

## Funding

This R package was developed as part of a dedicated study funded by the
U.S. Bureau of Ocean Energy Management (BOEM, Contract 140M0121C0008).

## Acknowledgements

We are grateful for the technical support received from the following
individuals and organizations, all of whom provided data and/or
analytical expertise: NOAA’s National Marine Fisheries Service (Eric
Patterson, Lance Garrison, Jeff Adams, Hannah Blondin, Dan Linden, Laura
Solinger), the North Atlantic Right Whale Consortium, the New England
Aquarium (Heather Pettis, Amy Knowlton, Philip Hamilton), WHOI (Carolyn
Miller), IFAW (Sarah Sharp), the MERIDIAN Initiative at Dalhousie
University (Romina Gehrmann, Matthew Smith), Fisheries and Oceans Canada
(Caroline Lehoux, Stéphane Plourde), the North Atlantic Fisheries
Organization, and the ONR-SERDP funded PCOMS project. We also thank
Susanna Blackwell, Susan Parks, Elizabeth Henderson, Brandon Southall,
Paul Wensveen, and Aimee Darias O’Hara for their contributions to the
expert elicitation on right whale sensitivity to piling noise.

## Relevant literature

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Bailey2014" class="csl-entry">

Bailey, Helen, Kate L. Brookes, and Paul M. Thompson. 2014.
“<span class="nocase">Assessing environmental impacts of offshore wind
farms: lessons learned and recommendations for the future</span>.”
*Aquatic Biosystems* 10: 8. <https://doi.org/10.1186/2046-9063-10-8>.

</div>

<div id="ref-Costa2012" class="csl-entry">

Costa, Daniel P. 2012. “<span class="nocase">A Bioenergetics Approach to
Developing a Population Consequences of Acoustic Disturbance
Model</span>.” In *<span class="nocase">The Effects of Noise on Aquatic
Life</span>*, 423–26. New York, USA.
<https://doi.org/10.1007/978-1-4419-7311-5_96>.

</div>

<div id="ref-Crain2008" class="csl-entry">

Crain, Caitlin Mullan, Kristy Kroeker, and Benjamin S. Halpern. 2008.
“<span class="nocase">Interactive and cumulative effects of multiple
human stressors in marine systems</span>.” *Ecology Letters* 11 (12):
1304–15. <https://doi.org/10.1111/j.1461-0248.2008.01253.x>.

</div>

<div id="ref-DeCastro2019" class="csl-entry">

DeCastro, M., S. Salvador, M. Gómez-Gesteira, X. Costoya, D. Carvalho,
F. J. Sanz-Larruga, and L. Gimeno. 2019. “<span class="nocase">Europe,
China and the United States: Three different approaches to the
development of offshore wind energy</span>.” *Renewable and Sustainable
Energy Reviews* 109: 55–70.
<https://doi.org/10.1016/j.rser.2019.04.025>.

</div>

<div id="ref-Goodale2019" class="csl-entry">

Goodale, M. Wing, and Anita Milman. 2019.
“<span class="nocase">Assessing the cumulative exposure of wildlife to
offshore wind energy development</span>.” *Journal of Environmental
Management* 235: 77–83. <https://doi.org/10.1016/j.jenvman.2019.01.022>.

</div>

<div id="ref-Linden2023" class="csl-entry">

Linden, Daniel. 2023. “<span class="nocase">Population size estimation
of North Atlantic right whales from 1990-2022</span>.” 314. US Dept
Commer Northeast Fish Sci Cent Tech Memo.

</div>

<div id="ref-Moore2021" class="csl-entry">

Moore, Michael J., Teresa K. Rowles, Deborah A. Fauquier, Jason D.
Baker, Ingrid Biedron, John W. Durban, Philip K. Hamilton, et al. 2021.
“<span class="nocase">REVIEW: Assessing North Atlantic right whale
health: Threats, and development of tools critical for conservation of
the species</span>.” *Diseases of Aquatic Organisms* 143: 205–26.
<https://doi.org/10.3354/dao03578>.

</div>

<div id="ref-Petruny2014" class="csl-entry">

Petruny, Loren M., Andrew J. Wright, and Courtney E. Smith. 2014.
“<span class="nocase">Getting it right for the North Atlantic right
whale (Eubalaena glacialis): A last opportunity for effective marine
spatial planning?</span>” *Marine Pollution Bulletin* 85 (1): 24–32.
<https://doi.org/10.1016/j.marpolbul.2014.06.004>.

</div>

<div id="ref-Pirotta2018" class="csl-entry">

Pirotta, Enrico, Cormac G. Booth, Daniel P. Costa, Erica Fleishman,
Scott D. Kraus, David Lusseau, David Moretti, et al. 2018.
“<span class="nocase">Understanding the population consequences of
disturbance</span>.” *Ecology and Evolution* 8 (19): 9934–46.
<https://doi.org/10.1002/ece3.4458>.

</div>

<div id="ref-Pirotta2019" class="csl-entry">

Pirotta, Enrico, Marc Mangel, Daniel P. Costa, Jeremy Goldbogen, John
Harwood, Vincent Hin, Ladd M. Irvine, et al. 2019.
“<span class="nocase">Anthropogenic disturbance in a changing
environment: modelling lifetime reproductive success to predict the
consequences of multiple stressors on a migratory population</span>.”
*Oikos* 128 (9): 1340–57. <https://doi.org/10.1111/oik.06146>.

</div>

<div id="ref-Pirotta2022" class="csl-entry">

Pirotta, Enrico, Len Thomas, Daniel P. Costa, Ailsa J. Hall, Catriona M.
Harris, John Harwood, Scott D. Kraus, et al. 2022.
“<span class="nocase">Understanding the combined effects of multiple
stressors: A new perspective on a longstanding challenge</span>.”
*Science of The Total Environment* 821: 153322.
<https://doi.org/10.1016/j.scitotenv.2022.153322>.

</div>

<div id="ref-Pace2017" class="csl-entry">

Richard M. Pace, III, Peter J. Corkeron, and Scott D. Kraus. 2017.
“<span class="nocase">State–space mark–recapture estimates reveal a
recent decline in abundance of North Atlantic right whales</span>.”
*Ecology and Evolution* 7 (21): 8730.
<https://doi.org/10.1002/ece3.3406>.

</div>

<div id="ref-Rolland2016" class="csl-entry">

Rolland, Rosalind M., Robert S. Schick, Heather M. Pettis, Amy R.
Knowlton, Philip K. Hamilton, James S. Clark, and Scott D. Kraus. 2016.
“<span class="nocase">Health of North Atlantic right whales Eubalaena
glacialis over three decades: from individual health to demographic and
population health trends</span>.” *Marine Ecology Progress Series* 542:
265–82. <https://doi.org/10.3354/meps11547>.

</div>

<div id="ref-Schick2013" class="csl-entry">

Schick, Robert S., Scott D. Kraus, Rosalind M. Rolland, Amy R. Knowlton,
Philip K. Hamilton, Heather M. Pettis, Robert D. Kenney, and James S.
Clark. 2013. “<span class="nocase">Using Hierarchical Bayes to
Understand Movement, Health, and Survival in the Endangered North
Atlantic Right Whale</span>.” *PLOS ONE* 8 (6): e64166.
<https://doi.org/10.1371/journal.pone.0064166>.

</div>

<div id="ref-Thompson2010" class="csl-entry">

Thompson, Paul M., David Lusseau, Tim Barton, Dave Simmons, Jan Rusin,
and Helen Bailey. 2010. “<span class="nocase">Assessing the responses of
coastal cetaceans to the construction of offshore wind turbines</span>.”
*Marine Pollution Bulletin* 60 (8): 1200–1208.
<https://doi.org/10.1016/j.marpolbul.2010.03.030>.

</div>

<div id="ref-Verfuss2016" class="csl-entry">

Verfuss, Ursula K., Carol E. Sparling, Charlie Arnot, Adrian Judd, and
Michael Coyle. 2016. “<span class="nocase">Review of Offshore Wind Farm
Impact Monitoring and Mitigation with Regard to Marine Mammals</span>.”
*Advances in Experimental Medicine and Biology* 875:1175-82. (;): Adv.
<https://doi.org/10.1007/978-1-4939-2981-8_147>.

</div>

</div>
