.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "-----------------------------------------------------------------------------------\n",
    "Population Effects of Offshore Wind Development on North Atlantic Right Whales\n",
    "-----------------------------------------------------------------------------------\n",
    "narwind: version 1.0\n",
    "Package developed under funding from the U.S. Bureau of Ocean Energy Management\n",
    "(BOEM Contract No. 140M0121C0008)\n",
    "\n",
    "For more information, see the package vignette at: URL "
                        )
 options(tibble.width = Inf)
 options(pillar.sigfig = 5)
}
