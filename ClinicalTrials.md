---
name: ClinicalTrials
topic: Clinical Trial Design, Monitoring, and Analysis
maintainer: Ed Zhang, W. G. Zhang, R. G. Zhang
email: ClinicalTrials.TaskView@yahoo.com
version: 2021-12-29
source: https://github.com/cran-task-views/ClinicalTrials/
---


This task view gathers information on specific R packages for design,
monitoring and analysis of data from clinical trials. It focuses on
including packages for clinical trial design and monitoring in general
plus data analysis packages for a specific type of design. Also, it
gives a brief introduction to important packages for analyzing clinical
trial data. Please refer to task views
`r view("ExperimentalDesign")`,
`r view("Survival")`,
`r view("Pharmacokinetics")`, `r view("Meta-analysis")` for more details on these
topics.

Contributions are always welcome and encouraged, either via e-mail to the
maintainer or by submitting an issue or pull request in the GitHub
repository linked above.


### Design and Monitoring

-   `r pkg("TrialSize", priority = "core")` This package has
    more than 80 functions from the book *Sample Size Calculations in
    Clinical Research* (Chow & Wang & Shao, 2007, 2nd ed., Chapman
    &Hall/CRC).
-   `r pkg("asd", priority = "core")` This Package runs
    simulations for adaptive seamless designs using early outcomes for
    treatment selection.
-   `r pkg("bcrm", priority = "core")` This package
    implements a wide variety of one and two-parameter Bayesian CRM
    designs. The program can run interactively, allowing the user to
    enter outcomes after each cohort has been recruited, or via
    simulation to assess operating characteristics.
-   `r pkg("blockrand", priority = "core")` creates
    randomizations for block random clinical trials. It can also produce
    a PDF file of randomization cards.
-   `r pkg("clusterPower")` Calculate power for cluster
    randomized trials (CRTs) that compare two means, two proportions, or
    two counts using closed-form solutions. In addition, calculate power
    for cluster randomized crossover trials using Monte Carlo methods.
    For more information, see Reich et al. (2012)
    [doi:10.1371/journal.pone.0035564](https://dx.doi.org/10.1371/journal.pone.0035564)
    .
-   `r pkg("conf.design")` This small package contains a
    series of simple tools for constructing and manipulating confounded
    and fractional factorial designs.
-   `r pkg("crmPack")` Implements a wide range of
    model-based dose escalation designs, ranging from classical and
    modern continual reassessment methods (CRMs) based on dose-limiting
    toxicity endpoints to dual-endpoint designs taking into account a
    biomarker/efficacy outcome. The focus is on Bayesian inference,
    making it very easy to setup a new design with its own JAGS code.
    However, it is also possible to implement 3+3 designs for comparison
    or models with non-Bayesian estimation. The whole package is written
    in a modular form in the S4 class system, making it very flexible
    for adaptation to new models, escalation or stopping rules.
-   `r pkg("cosa")` Implements bound constrained optimal
    sample allocation (BCOSA) framework described in Bulus & Dong (2019)
    for power analysis of multilevel regression discontinuity designs
    (MRDDs) and multilevel randomized trials (MRTs) with continuous
    outcomes. Separate tools for statistical power and minimum
    detectable effect size computations are provided.
-   `r pkg("dfcrm", priority = "core")` This package provides
    functions to run the CRM and TITE-CRM in phase I trials and
    calibration tools for trial planning purposes.
-   `r pkg("DTAT")` Dose Titration Algorithm Tuning (DTAT)
    is a methodologic framework allowing dose individualization to be
    conceived as a continuous learning process that begins in
    early-phase clinical trials and continues throughout drug
    development, on into clinical practice. This package includes code
    that researchers may use to reproduce or extend key results of the
    DTAT research programme, plus tools for trialists to design and
    simulate a '3+3/PC' dose-finding study.
-   `r pkg("ewoc")` An implementation of a variety of
    escalation with overdose control designs introduced by Babb, Rogatko
    and Zacks (1998)
    [doi:10.1002/(SICI)1097-0258(19980530)17:10%3C1103::AID-SIM793%3E3.0.CO;2-9](https://dx.doi.org/10.1002/(SICI)1097-0258(19980530)17:10%3C1103::AID-SIM793%3E3.0.CO;2-9)
    . It calculates the next dose as a clinical trial proceeds as well
    as performs simulations to obtain operating characteristics.
-   `r pkg("experiment", priority = "core")` contains tools
    for clinical experiments, e.g., a randomization tool, and it
    provides a few special analysis options for clinical trials.
-   `r pkg("FrF2")` This package creates regular and
    non-regular Fractional Factorial designs. Furthermore, analysis
    tools for Fractional Factorial designs with 2-level factors are
    offered (main effects and interaction plots for all factors
    simultaneously, cube plot for looking at the simultaneous effects of
    three factors, full or half normal plot, alias structure in a more
    readable format than with the built-in function alias). The package
    is currently subject to intensive development. While much of the
    intended functionality is already available, some changes and
    improvements are still to be expected.
-   `ldBand` from `r pkg("Hmisc", priority = "core")`
    computes and plots group sequential stopping boundaries from the
    Lan-DeMets method with a variety of a-spending functions using the
    ld98 program from the Department of Biostatistics, University of
    Wisconsin written by DM Reboussin, DL DeMets, KM Kim, and KKG Lan.
-   `r pkg("ldbounds", priority = "core")` uses Lan-DeMets
    Method for group sequential trial; its functions calculate bounds
    and probabilities of a group sequential trial.
-   `r pkg("longpower", priority = "core")`Compute power and sample size 
    for linear models of longitudinal data. The package is described in 
    Iddi and Donohue (2022) <doi:10.32614/RJ-2022-022>.  
-   `r pkg("Mediana")` Provides a general framework for
    clinical trial simulations based on the Clinical Scenario Evaluation
    (CSE) approach. The package supports a broad class of data models
    (including clinical trials with continuous, binary, survival-type
    and count-type endpoints as well as multivariate outcomes that are
    based on combinations of different endpoints), analysis strategies
    and commonly used evaluation criteria.
-   `r pkg("PowerTOST", priority = "core")` contains
    functions to calculate power and sample size for various study
    designs used for bioequivalence studies. See function
    known.designs() for study designs covered. Moreover the package
    contains functions for power and sample size based on 'expected'
    power in case of uncertain (estimated) variability. Added are
    functions for the power and sample size for the ratio of two means
    with normally distributed data on the original scale (based on
    Fieller's confidence ('fiducial') interval).
-   `r pkg("MinEDfind")` The nonparametric two-stage
    Bayesian adaptive design is a novel phase II clinical trial design
    for finding the minimum effective dose (MinED). This design is
    motivated by the top priority and concern of clinicians when testing
    a new drug, which is to effectively treat patients and minimize the
    chance of exposing them to subtherapeutic or overly toxic doses. It
    is used to design single-agent trials.
-   `r pkg("presize")` Bland (2009) recommended to base
    study sizes on the width of the confidence interval rather the power
    of a statistical test. The goal of 'presize' is to provide
    functions for such precision based sample size calculations. For a
    given sample size, the functions will return the precision (width of
    the confidence interval), and vice versa.
-   `r pkg("PowerUpR")` Includes tools to calculate
    statistical power, minimum detectable effect size (MDES), MDES
    difference (MDESD), and minimum required sample size for various
    multilevel randomized experiments with continuous outcomes. Some of
    the functions can assist with planning two- and three-level
    cluster-randomized trials (CRTs) sensitive to multilevel moderation
    and mediation (2-1-1, 2-2-1, and 3-2-1).
-   `r pkg("pwr", priority = "core")` has power analysis
    functions along the lines of Cohen (1988).
-   `r pkg("randomizeR")` This tool enables the user to
    choose a randomization procedure based on sound scientific criteria.
    It comprises the generation of randomization sequences as well the
    assessment of randomization procedures based on carefully selected
    criteria. Furthermore, 'randomizeR' provides a function for the
    comparison of randomization procedures.
-   `r pkg("replicateBE")` Performs comparative
    bioavailability calculations for Average Bioequivalence with
    Expanding Limits (ABEL). Implemented are 'Method A' and 'Method
    B' and the detection of outliers. If the design allows, assessment
    of the empiric Type I Error and iteratively adjusting alpha to
    control the consumer risk. Average Bioequivalence - optionally with
    a tighter (narrow therapeutic index drugs) or wider acceptance range
    (Gulf Cooperation Council, South Africa: Cmax) - is implemented as
    well.
-   `r pkg("rpact")` Design and analysis of confirmatory
    adaptive clinical trials with continuous, binary, and survival
    endpoints according to the methods described in the monograph by
    Wassmer and Brannath (2016). This includes classical group
    sequential as well as multi-stage adaptive hypotheses tests that are
    based on the combination testing principle.
-   `r pkg("samplesize")` computes sample size for
    Student's t-test with equal and nonequal variances and for the
    Wilcoxon-Mann-Whitney test for categorical data with and without
    ties.
-   `r pkg("simglm")` Simulates regression models, including
    both simple regression and generalized linear mixed models with up
    to three level of nesting. Power simulations that are flexible
    allowing the specification of missing data, unbalanced designs, and
    different random error distributions are built into the package.
-   `r pkg("UnifiedDoseFinding")` In many phase I trials,
    the design goal is to find the dose associated with a certain target
    toxicity rate. In some trials, the goal can be to find the dose with
    a certain weighted sum of rates of various toxicity grades. For
    others, the goal is to find the dose with a certain mean value of a
    continuous response. This package provides the setup and
    calculations needed to run a dose-finding trial with non-binary
    endpoints and performs simulations to assess design's operating
    characteristics under various scenarios.

### Design and Analysis

-   Package `r pkg("AGSDest")` This package provides tools
    and functions for parameter estimation in adaptive group sequential
    trials.
-   Package `r pkg("clinfun", priority = "core")` has
    functions for both design and analysis of clinical trials. For phase
    II trials, it has functions to calculate sample size, effect size,
    and power based on Fisher's exact test, the operating
    characteristics of a two-stage boundary, Optimal and Minimax 2-stage
    Phase II designs given by Richard Simon, the exact 1-stage Phase II
    design and can compute a stopping rule and its operating
    characteristics for toxicity monitoring based repeated significance
    testing. For phase III trials, it can calculate sample size for
    group sequential designs.
-   Package `r pkg("CRM")` Continual Reassessment Method
    (CRM) simulator for Phase I Clinical Trials.
-   Package `r pkg("dfpk")` Statistical methods involving PK
    measures are provided, in the dose allocation process during a Phase
    I clinical trials. These methods enter pharmacokinetics (PK) in the
    dose finding designs in different ways, including covariates models,
    dependent variable or hierarchical models. This package provides
    functions to generate data from several scenarios and functions to
    run simulations which their objective is to determine the maximum
    tolerated dose (MTD).
-   Package `r pkg("dfped")` A unified method for designing
    and analysing dose-finding trials in paediatrics, while bridging
    information from adults, is proposed in the dfped package. The dose
    range can be calculated under three extrapolation methods: linear,
    allometry and maturation adjustment, using pharmacokinetic (PK)
    data. To do this, it is assumed that target exposures are the same
    in both populations. The working model and prior distribution
    parameters of the dose-toxicity and dose-efficacy relationships can
    be obtained using early phase adult toxicity and efficacy data at
    several dose levels through dfped package. Priors are used into the
    dose finding process through a Bayesian model selection or adaptive
    priors, to facilitate adjusting the amount of prior information to
    differences between adults and children. This calibrates the model
    to adjust for misspecification if the adult and paediatric data are
    very different. User can use his/her own Bayesian model written in
    Stan code through the dfped package. A template of this model is
    proposed in the examples of the corresponding R functions in the
    package. Finally, in this package you can find a simulation function
    for one trial or for more than one trial.
-   Package `r pkg("DoseFinding")` provides functions for
    the design and analysis of dose-finding experiments (for example
    pharmaceutical Phase II clinical trials). It provides functions for:
    multiple contrast tests, fitting non-linear dose-response models,
    calculating optimal designs and an implementation of the
    `r pkg("MCPMod")` methodology. Currently only normally
    distributed homoscedastic endpoints are supported.
-   `r pkg("MCPMod")` This package implements a methodology
    for the design and analysis of dose-response studies that combines
    aspects of multiple comparison procedures and modeling approaches
    (Bretz, Pinheiro and Branson, 2005, Biometrics 61, 738-748). The
    package provides tools for the analysis of dose finding trials as
    well as a variety of tools necessary to plan a trial to be conducted
    with the MCPMod methodology.
-   Package `r pkg("TEQR", priority = "core")` The target
    equivalence range (TEQR) design is a frequentist implementation of
    the modified toxicity probability interval (mTPI) design and a
    competitor to the standard 3+3 design (3+3). The 3+3 is the work
    horse design in Phase I. It is good at determining if a safe dose
    exits, but provides poor accuracy and precision in estimating the
    level of toxicity at the maximum tolerated dose (MTD). The TEQR is
    better than the 3+3 when compared on: 1) the number of times the
    dose at or nearest the target toxicity level was selected as the
    MTD, 2) the number of subjects assigned to doses levels, at or
    nearest the MTD, and 3) the overall trial DLT rate. TEQR more
    accurately and more precisely estimates the rate of toxicity at the
    MTD because a larger number of subjects are studied at the MTD dose.
    The TEQR on average uses fewer subjects and provide reasonably
    comparable results to the continual reassessment method (CRM) in the
    number of times the dose at or nearest the target toxicity level was
    selected as the MTD and the number of subjects assigned doses, at,
    or nearest the target and in overall DLT rate.
-   Package `r pkg("ThreeArmedTrials")` Design and analyze
    three-arm non-inferiority or superiority trials which follow a
    gold-standard design, i.e. trials with an experimental treatment, an
    active, and a placebo control.

### Analysis for Specific Designs

-   `r pkg("adaptTest", priority = "core")` The functions
    defined in this program serve for implementing adaptive two-stage
    tests. Currently, four tests are included: Bauer and Koehne (1994),
    Lehmacher and Wassmer (1999), Vandemeulebroecke (2006), and the
    horizontal conditional error function. User-defined tests can also
    be implemented. Reference: Vandemeulebroecke, An investigation of
    two-stage tests, Statistica Sinica 2006.
-   `r pkg("adaptr")` simulates adaptive (multi-arm, multi-stage) 
    clinical trials using adaptive stopping, adaptive arm dropping, and/or adaptive randomisation. 
-   `r pkg("clinsig")` This function calculates both
    parametric and non-parametric versions of the Jacobson-Truax
    estimates of clinical significance.
-   `r pkg("clinicalsignificance")` The goal of this package is to provide all 
   necessary tools for analyses of clinical significance in clinical intervention studies. 
   In contrast to statistical significance, which assesses if it is probable that there 
   is a treatment effect, clinical significance can be used to determine if a treatment 
   effect is of practical use or meaningful for patients.
-   `r pkg("nppbib")` implements a nonparametric statistical
    test for rank or score data from partially-balanced incomplete
    block-design experiments.
-   `r pkg("speff2trial", priority = "core")`, the package
    performs estimation and testing of the treatment effect in a 2-group
    randomized clinical trial with a quantitative or dichotomous
    endpoint.
-   `r pkg("ThreeGroups")` This package implements the
    Maximum Likelihood estimator for three-group designs proposed by
    Gerber, Green, Kaplan, and Kern (2010).

### Analysis in General

-   Base R, especially the stats package, has a lot of functionality
    useful for design and analysis of clinical trials. For example,
    `chisq.test`, `prop.test`, `binom.test`, `t.test`, `wilcox.test`,
    `kruskal.test`, `mcnemar.test`, `cor.test`, `power.t.test`,
    `power.prop.test`, `power.anova.test`, `lm`, `glm`, `nls`, `anova`
    (and its `lm` and `glm` methods) among many others.
-   `r pkg("accrualPlot")` Tracking accrual in clinical trials is important for trial success.
    'accrualPlot' provides functions to aid the tracking of accrual and predict when a trial 
    will reach it's intended sample size.
-   `r pkg("binomSamSize")` is a suite of functions for
    computing confidence intervals and necessary sample sizes for the
    success probability parameter Bernoulli distribution under simple
    random sampling or under pooled sampling.
-   `r pkg("coin")` offers conditional inference procedures
    for the general independence problem including two-sample, K-sample
    (non-parametric ANOVA), correlation, censored, ordered and
    multivariate problems.
-   `r pkg("ctrdata")` is a system for querying, retrieving and analyzing 
    protocol- and results-related information on clinical trials from four public registers
-   `r pkg("epibasix")` has functions such as `diffdetect`,
    `n4means` for continuous outcome and `n4props` and functions for
    matched pairs analysis in randomized trials.
-   `ae.dotplot` from `r pkg("HH")` shows a two-panel
    display of the most frequently occurring adverse events in the
    active arm of a clinical study.
-   The `r pkg("Hmisc")` package contains around 200
    miscellaneous functions useful for such things as data analysis,
    high-level graphics, utility operations, functions for computing
    sample size and power, translating SAS datasets into S, imputing
    missing values, advanced table making, variable clustering,
    character string manipulation, conversion of S objects to LaTeX
    code, recoding variables, and bootstrap repeated measures analysis.
-   `r pkg("mmrm")` Implements mixed models for repeated measures (MMRM), 
    a popular choice for analyzing longitudinal continuous outcomes in 
    randomized clinical trials and beyond.
-   `r pkg("multcomp")` covers simultaneous tests and
    confidence intervals for general linear hypotheses in parametric
    models, including linear, generalized linear, linear mixed effects,
    and survival models.
-   `r pkg("rbmi")` Implements standard and reference based multiple 
    imputation allowing for the imputation of longitudinal datasets using 
    predefined strategies. The package is described in Gower-Page et al (2022) 
    <doi: 10.21105/joss.04251>.
-   `r pkg("survival", priority = "core")` contains
    descriptive statistics, two-sample tests, parametric accelerated
    failure models, Cox model. Delayed entry (truncation) allowed for
    all models; interval censoring for parametric models. Case-cohort
    designs.
-   `r pkg("ssanv")` is a set of functions to calculate
    sample size for two-sample difference in means tests. Does
    adjustments for either nonadherence or variability that comes from
    using data to estimate parameters.

### Meta-Analysis

-   `r pkg("metasens")` is a package for statistical methods
    to model and adjust for bias in meta-analysis
-   `r pkg("meta")` is for fixed and random effects
    meta-analysis. It has Functions for tests of bias, forest and funnel
    plot.
-   `r pkg("metafor")` consists of a collection of functions
    for conducting meta-analyses. Fixed- and random-effects models (with
    and without moderators) can be fitted via the general linear
    (mixed-effects) model. For 2x2 table data, the Mantel-Haenszel and
    Peto's method are also implemented.
-   `r pkg("metaLik")` Likelihood inference in meta-analysis
    and meta-regression models.
-   `r pkg("rmeta")` has functions for simple fixed and
    random effects meta-analysis for two-sample comparisons and
    cumulative meta-analyses. Draws standard summary plots, funnel
    plots, and computes summaries and tests for association and
    heterogeneity.



### Links
-   [Regulatory Compliance and Validation Issues (A Guidance Document for the Use of R in Regulated Clinical Trial Environments)](http://www.R-project.org/doc/R-FDA.pdf)
