# Measuring and simulating belief in news articles at the scale of the social media ecosystem

This is the code for the forthecoming paper:
> Tokita CK, Aslett K, Godel WP, Sanderson Z, Tucker JA, Nagler J, Persily N, Bonneau RA. (In Prep). Measuring and mitigating belief in misinformation at the scale of the social media ecosystem.

Bibtex Citation:
```
 @article{
 tokita_aslett_godel_sanderson_tucker_nagler_persily_bonneau, 
 title={Measuring and mitigating belief in misinformation at the scale of the social media ecosystem}, 
 author={Tokita, Christopher K. and Aslett, Kevin and Godel, William P. and Sanderson, Zeve and Tucker, Joshua  A. and Nagler, Jonathan and Persily, Nathaniel and Bonneau, Richard A.}
  } 
 ```

## Project abstract
Using computational simulations, this project combines social media data and experimental survey data to create an estimate for exposure to and belief in top-trending news articles on Twitter.

The abstract from our paper:
> While online misinformation has become a growing topic of scholarly and public inquiry in recent years, the science of studying its impact and evaluating potential remedies at the scale of the social media ecosystem is still developing. On the one hand, research on the diffusion of and exposure to misinformation utilizes large-scale observational data, but these data provide little insight into belief. On the other, measurement of belief in misinformation is largely derived from surveys, which cannot capture the networked dynamics of belief's spread. These two methods have yet to be unified, leaving us without an ecosystem-level estimate of levels of belief in misinformation and the efficacy of interventions aimed at mitigating belief at scale. Here, we pair experimental survey data with observational Twitter data to estimate both exposure to and belief in trending true and false news. Using this new method, we find that the patterns of exposure and belief differ between true and false news: while true news is seen and believed by an ideologically balanced segment of the user base, false news is seen by an ideologically diverse set of users yet believed by an ideologically-skewed segment of the user base. Thus, efforts to infer the impact of misinformation by measuring exposure alone may not accurately capture the true impact of misinformation on social media. We extend this new method and conduct data-driven simulation to evaluate different interventions social media platforms deploy to curtail misinformation. We find that most interventions have a minimal impact on the number of people who see and believe misinformation and, importantly, that the effectiveness of interventions quickly dissipates with each hour delay in implementation. Our paper provides the first full-scale estimation of belief in and potential remedies for online misinformation at the scale of social media platforms.


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
* `scripts/analysis_functions`: contains the main `tweet_parser()` function that parses the tweet Json files.
* `scripts/_plot_themes`: contains the custom ggplot2 and plotnine themes for quicker plot formatting.

### Installing necessary packages
* **Python**: All Python packages for this project can be found in `<requirements.txt>` and can easily be installed with `<pip install -r requirements.txt>` from command line.
* **R**: dplyr, ggplot2, tidyr, stringr, RColorBrewer, scales, brms

### Data
Given the large size of the Twitter data, all the data in this study is currently held on an external harddrive, with a backup copy stored on Princeton's Della cluster.
