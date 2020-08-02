# Connectedness and pushed/pulled range expansions


#### 

\[!\[DOI\](https://zenodo.org/badge/DOI/10.5281/zenodo.3969988.svg)\](https://doi.org/10.5281/zenodo.3969988)

 This repo contains all data and code needed to re-do the analyses and figures in our manuscript

"Shifts from pulled to pushed range expansions caused by reduction of landscape connectedness"
(by Maxime Dahirel, Aline Bertin, Marjorie Haond, Aurélie Blin, Eric Lombaert, Vincent Calcagno, Simon Fellous, Ludovic Mailleret, Thibaut Malausa, Elodie Vercken)

(link to bioRxiv preprint [here](https://doi.org/10.1101/2020.05.13.092775))

- raw experimental data in `csv` format are in the `data` folder
- the source code for the Netlogo model in the `Netlogo_model` folder
- the R scripts (including detailed information about the analyses) are in the `R` folder. There are four `Rmd` files: 
    - the first one (`Trichogramma-2020_1_experiment-main-text`) corresponds to the analysis of the experimental data
    - the second one (`Trichogramma-2020_2_ibm-creation.Rmd`) uses the `nlrx` package and the Netlogo model to run simulations and save some metrics of interest to `Netlogo_output`
    - the third one (`Trichogramma-2020_3_ibm-analysis-main-text`) corresponds to the main analysis of these simulated data
    - the fourth one corresponds to the "Supplementary Material" file. It can be run just as the other scripts, or knitted to produce the Supplementary Material `html` file. It is in its own subfolder along with bibliography files used when knitting.

This folder is a RStudio project folder, and the script uses the [`here`](https://here.r-lib.org/) package (see also [here](https://github.com/jennybc/here_here)). This means all file paths are relative, and the analysis should work on your computer no questions asked, whether you use the R project or not, no code line to change as long as you download the entire repository (you just need to install all the needed packages first, of course).

If you run the first, third or fourth`R`scripts for the first time, models and some other time-consuming outputs will be saved as `RData` files in the `R_output` folder to be reused as needed, saving you some time later (you can bypass this behaviour to re-run models manually if you want; see code for details).
