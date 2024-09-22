# Measuring receptivity to misinformation at scale on a social media platform

This is the code for the paper:
> Tokita CK, Aslett K, Godel WP, Sanderson Z, Tucker JA, Nagler J, Persily N, Bonneau RA. (2024). Measuring receptivity to misinformation at scale on a social media platform. PNAS Nexus.

Bibtex Citation:
```
@article{tokita2024measuring,
  title={Measuring receptivity to misinformation at scale on a social media platform},
  author={Tokita, Christopher K and Aslett, Kevin and Godel, William P and Sanderson, Zeve and Tucker, Joshua A and Nagler, Jonathan and Persily, Nathaniel and Bonneau, Richard},
  journal={PNAS Nexus},
  year={2024},
  publisher={Oxford University Press}
}
 ```
This is the living version of the code, although I don't expect many updates to take place going forward. The officially archived code for the paper as it was published can be found at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13826419.svg)](https://doi.org/10.5281/zenodo.13826419)

An archived copy of the data used in this project can be found at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13777171.svg)](https://doi.org/10.5281/zenodo.13777171). See the section "Data" below for more details.

## Project abstract
Using computational simulations, this project combines social media data and experimental survey data to create an estimate for exposure to and belief in top-trending news articles on Twitter.

The abstract from our paper:
> Measuring the impact of online misinformation is challenging. Traditional measures, such as user views or shares on social media, are incomplete because not everyone who is exposed to misinformation is equally likely to believe it. To address this issue, we developed a method that combines survey data with observational Twitter data to probabilistically estimate the number of users both exposed to and likely to believe a specific news story. As a proof of concept, we applied this method to 139 viral news articles and find that although false news reaches an audience with diverse political views, users who are both exposed and receptive to believing false news tend to have more extreme ideologies. These receptive users are also more likely to encounter misinformation earlier than those who are unlikely to believe it. This mismatch between overall user exposure and receptive user exposure underscores the limitation of relying solely on exposure or interaction data to measure the impact of misinformation, as well as the challenge of implementing effective interventions. To demonstrate how our approach can address this challenge, we then conducted data-driven simulations of common interventions used by social media platforms. We find that these interventions are only modestly effective at reducing exposure among users likely to believe misinformation, and their effectiveness quickly diminishes unless implemented soon after misinformation's initial spread. Our paper provides a more precise estimate of misinformation's impact by focusing on the exposure of users likely to believe it, offering insights for effective mitigation strategies on social media.

The significane statement from our paper:
> As social media platforms grapple with misinformation, our study offers a new approach to measure its spread and impact. By combining survey data with social media data, we estimate not only the number of users exposed to false (and true) news but also the number of users likely to believe these news stories. We find that the impact of misinformation is not evenly distributed, with ideologically extreme users being more likely to see and believe false content, often encountering it before others. Our simulations suggest that current interventions may have limited effectiveness in reducing the exposure of receptive users. These findings highlight the need to consider individual user receptiveness when measuring misinformation's impact and developing policies to combat its spread.

## About this repository and how to use the code
Most of this project is written in Python. Python scripts are used to construct the dataset, calculate exposure/belief, run simulations of platform interventions, and analyze the retweet network structure. R is used for generating figures/plots and summary statistics for the paper.

### Structure of the repository
* `scripts/`: contains all code for this project.
* `output/`: contains plots generated from the R scripts.
* `slurm_scripts/`: shell scripts needed to run python scripts on Princeton's high-performing computering clusters.
* `README.md`
* `requirements.txt`: required python packages (see subsection "Installing necessary packages" in this repository).
* `plotcolors.txt`: lists the exact hex codes for colors to denote user ideology in the plots.
* `news-belief-at-scale.Rproj`: R project for easier running of plotting scripts.

### Structure of the code

#### The main analysis and data visualization scripts
Python scripts are named with a leading number that indicates the order in which they should be run, e.g., `01_parse_tweets.py`. Thus, the python script are structured in the analysis pipeline for the paper. The general steps of the pipeline are as follows:

1. Parsing tweet data
2. Parsing twitter user data
3. Constructing the full twitter dataset by adding tweet/user metadata
4. Inferring twitter user ideologies
5. Estimating exposure to news articles
6. Estimating belief in news articles (with use of experimental survey data)
7. Constructing and analyzing retweet networks
8. Simulating platform interventions on misinformation

Some steps in the pipeline have multipel substeps, which is reflected in the script numbering, e.g., `07a_construct_retweet_networks.py` and `07b_analyze_retweet_networks.py`

All R scripts for generating plots and summary statistics are denoted by the leading part of the script name `plot_{description of script focus}.R`.

#### Other folders in `scripts/`
* `scripts/suppl_code`: scripts for supplemental/ad hoc analysis.
* `scripts/_plot_themes`: contains the custom ggplot2 and plotnine themes for quicker plot formatting.

### Installing necessary packages
* **Python**: All Python packages for this project can be found in `<requirements.txt>` and can easily be installed with `<pip install -r requirements.txt>` from command line.
* **R**: dplyr, ggplot2, tidyr, stringr, RColorBrewer, scales, brms

### Structure of folders in `output/`
The `output/` folder contains all of the results of our analysis and visualization scripts (almost exclusively our `.R` scripts). Generally, the results are plots and data visualizations.

* `belief/`: estimates for the scale of user belief in news articles (subfolders break this down by article veracity and news source type).
* `exposure/`: estimates for the scale of user exposure to news articles (subfolders break this down by article veracity and news source type).
* `ideology_basis/`: comparison of self-reported user ideology (e.g., "liberal", "very liberal", etc.) and estimated ideology from Twitter data. This serves as a check on our basis for aligning self-reported ideology categories with the continuous ideology scores estimated from Twitter data.
* `interventions/`: the impact of different platform-level interventions on user exposure/belief in misinformation artilces in our dataset.
* `networks/`: visualization retweetn networks and analyzing their structure (includes break down by article veracity and news source type).
* `tweeter_ideology/`: distribution of ideologies among users who tweeted articles, along with estimates for the ideology of their followers (subfolders break this down by article veracity and news source type).
* `tweeting`: how users shared articles on Twitter via tweets (subfolders break this down by article veracity and news source type).

## Data
Given the large size of the Twitter data, all the data in this study was held on an external harddrive. Therefore, the scripts expect to be pointed to the data storage path, `data_directory` that has a `data/` folder containing all raw data and `data_derived/` folder with all resulting data from analyzing and processing the raw data. A copy of this entire data directory has been uploaded to Zenodo at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13777171.svg)](https://doi.org/10.5281/zenodo.13777171)


The data directory is expected to have the following structure:

```
data_directory/
│
└───data/
│   └───article_belief/
│   └───articles/
│   └───followers/
│   └───friends/
│   └───ideology_check/
│   └───tweets/
│   
└───data_derived/
    └───articles/
    └───belief/
    └───exposure/
    └───friends/
    └───ideological_scores/
    └───interventions/
    └───networks/
    └───tweets/
```
