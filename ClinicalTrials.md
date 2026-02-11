---
name: ClinicalTrials
topic: Clinical Trial Design, Monitoring, Analysis and Reporting
maintainer: Ya Wang, Thomas Jaki, Laura Pascasio Harris, Orla Doyle, Elias Laurin Meyer, Wilmar Igl
email: ya.wang10@gilead.com
version: 2026-02-11
source: https://github.com/cran-task-views/ClinicalTrials/
---


### Get Started

This task view provides an overview of R packages relevant to the design, monitoring, analysis and reporting of clinical trial data. 

Packages are grouped in the following categories:

- [**Design**](#design): tools to support a variety of clinical trial designs. The packages are further categorized into subgroups, such as [*adaptive designs*](#adaptive-designs), [*bioequivalence study designs*](#bioequivalence), [*dose-finding designs*](#dose-finding), [*factorial designs*](#factorial-designs), [*group sequential designs*](#group-sequential-designs), [*randomization*](#randomization), [*response adaptive randomization*](#response-adaptive-randomization), [*sample size and power calculations*](#sample-size-and-power-calculations), and [*simulation*](#simulation) for clinical trial designs.

- [**Monitoring**](#monitoring): tools dedicated to monitoring of clinical trials, such as computing the probability of crossing sequential boundaries along its interim analyses, sample size re-estimation, etc.

- [**Analysis**](#analysis): tools for implementing commonly used analysis method in clinical trials. The packages are further categorized into subgroups, such as [*general analysis*](#general-analysis), [*longitudinal data analysis*](#longitudinal-data-analysis), [*survival analysis*](#survival-analysis), [*meta-analysis*](#meta-analysis), [*missing data imputation*](#missing-data-imputation), as well as [*analysis for specific designs*](#other-analysis-for-specific-designs).

- [**Reporting**](#reporting): tools to facilitate the reporting of clinical trial results.


Here are several foundational books on clinical trial design and analysis that can help users gain a deeper understanding of the methods implemented in the relevant R packages: 

- *Clinical Trials: A Practical Approach* by Pocock (2013) `r doi("10.1002/9781118793916")`.

- *Fundamentals of Clinical Trials* by Friedman et al. (2015) `r doi("10.1007/978-3-319-18539-2")`.

- *Sample Sizes for Clinical Trials* by Julious (2023) `r doi("10.1201/9780429503658")`.

- *Group Sequential and Confirmatory Adaptive Designs in Clinical Trials* by Wassmer and Brannath (2016) `r doi("10.1007/978-3-319-32562-0")`.

- *Bayesian Adaptive Methods for Clinical Trials* by Berry et al. (2010) `r doi("10.1201/EBK1439825488")`.

In 2021, the R Foundation for Statistical Computing published the guidance document [*Regulatory Compliance and Validation Issues (A Guidance Document for the Use of R in Regulated Clinical Trial Environments)*](https://www.R-project.org/doc/R-FDA.pdf), which outlines its consensus on the use of R in GxP-regulated environments. This document provides a structured framework to help end users align with internal procedures, fulfill documentation requirements, and comply with regulatory standards.


### Inclusion Criteria

The packages included in this task view were carefully curated through a combination of expert recommendations from the pharmaceutical industry and academia, as well as an automated CRAN search using the `pkgsearch::pkg_search()` function. The search leveraged keywords aligned with our grouping categories, such as *clinical trials, study design, adaptive design, sample size calculation*, etc. Packages were considered within scope if they offered tools to facilitate the design, monitoring, or analysis of clinical trials. 

Some task views may include packages that are also relevant to clinical trials and will be listed within the grouping categories. Please refer to task views `r view("ExperimentalDesign")`, `r view("Meta-analysis")`, `r view("MissingData")`, `r view("Pharmacokinetics")`, `r view("Survival")` for a more comprehensive list of R packages related to these topics.  Note that while the `r view("Pharmacokinetics")` task view packages are closely related to clinical trials, they are not explicitly listed under the grouping categories in this task view to avoid duplication.

Contributions are always welcome and encouraged. You can contribute by emailing the maintainer directly or by submitting an issue or pull request in the GitHub repository linked above.  For further details see the [Contributing guide](https://github.com/cran-task-views/ctv/blob/main/Contributing.md).


### Design

#### *Adaptive Designs*

- `r pkg("adaptr")` facilitates simulation and comparison of adaptive clinical trial designs. It supports a flexible number of arms, use of a common control arm, pre-specified and user-defined outcome- and posterior probability distribution-generating functions, fixed- and response-adaptive randomisation, various adaptation rules for arm dropping and stopping, calculation of trial design performance metrics, and visualisation of results.
    
- `r pkg("adaptTest", priority = "core")` The functions defined in this program serve for implementing adaptive two-stage tests.  Currently, four tests are included, each based on Bauer and Koehne (1994) `r doi("10.2307/2533441")`, Lehmacher and Wassmer (1999) `r doi("10.1111/j.0006-341X.1999.01286.x")`, [Vandemeulebroecke (2006)](https://www.jstor.org/stable/24307582), and the horizontal conditional error function. User-defined tests can also be implemented.
    
- `r pkg("adestr")` provides methods to evaluate the performance characteristics of various point and interval estimators for optimal adaptive two-stage designs. Specifically, this package is written to work with trial designs created by the `adoptr` package based on Kunzmann et al. (2021) `r doi("10.18637/jss.v098.i09")` and Pilz et al. (2021) `r doi("10.1002/sim.8953")`.

- `r pkg("adpss")` provides the functions for planning and conducting a clinical trial with adaptive sample size determination. Maximal statistical efficiency will be exploited even when dramatic or multiple adaptations are made. 

- `r pkg("asd", priority = "core")` runs simulations for adaptive seamless designs with and without early outcomes for treatment selection and subpopulation type designs. It allows sample size modification in subpopulation selection.

- `r pkg("ASSISTant")` Clinical trial design for subgroup selection in three-stage group sequential trial as described in Lai et al. (2014) `r doi("10.1016/j.cct.2014.09.001")`. It includes facilities for design, exploration and analysis of such trials.

- `r pkg("BDP2")` Tools and workflow to choose design parameters in Bayesian adaptive single-arm phase II trial designs with binary endpoint (response, success) with possible stopping for efficacy and futility at interim analyses; and also contains routines to determine and visualize operating characteristics. See Kopp-Schneider et al. (2018) `r doi("10.1002/bimj.201700209")`. 

- `r pkg("cats")` simulates a cohort platform trial design whereby every cohort consists of two arms (control and experimental treatment). Endpoints are co-primary binary endpoints and decisions are made using either Bayesian or frequentist decision rules; and realistic trial trajectories are simulated with the operating characteristics of the designs calculated.

- `r pkg("CohortPlat")` is a collection of functions dedicated to simulating staggered entry platform trials whereby the treatment under investigation is a combination of two active compounds. A more detailed description of the design can be found in Meyer et al. `r doi("10.1002/pst.2194")` and a manual in Meyer et al. `r doi("10.48550/arXiv.2202.02182")`.

- `r pkg("esDesign")` is developed to implement the adaptive enrichment designs with sample size re-estimation presented in Lin et al. (2021) `r doi("10.1016/j.cct.2020.106216")`. In details, three-proposed trial designs are provided, including the AED1-SSR (or ES1-SSR), AED2-SSR (or ES2-SSR), AED3-SSR (or ES3-SSR); additionally, several widely used adaptive designs, such as the Marker Sequential Test (MaST) design proposed Freidlin et al. (2014) `r doi("10.1177/1740774513503739")`, the adaptive enrichment designs without early stopping (AED or ES), the sample size re-estimation procedure (SSR) based on the conditional power proposed by Proschan and Hunsberger (1995) `r doi("10.2307/2533262")`.

- `r pkg("eselect")` Endpoint selection and sample size reassessment for multiple binary endpoints based on blinded and/or unblinded data.  The implemented design is proposed in Roig et al. (2022) `r doi("10.48550/arXiv.2206.09639")`.

- `r pkg("gMCP")` provides functions and a graphical user interface for graphical described multiple test procedures. Examples of weighted tests that are available in gMCP are the weighted Bonferroni, parametric and Simes tests.

- `r pkg("graphicalMCP")` is a low-dependency implementation of graphical MCPs which allow mixed types of tests. It also includes power simulations and visualization of graphical MCPs.

- `r pkg("gsMAMS")` It provides functions to generate operating characteristics and to calculate Sequential Conditional Probability Ratio Tests(SCPRT) efficacy and futility boundary values along with sample/event size of Multi-Arm Multi-Stage(MAMS) trials for different outcomes. The package is based on Wu et al. (2023) `r doi("10.1002/sim.9682")`, Wu and Li (2023) `r doi(" 10.1002/sim.9682")`, and Wu et al. (2023) *Group Sequential Multi-Arm Multi-Stage Trial Design with Ordinal Endpoints* (In preparation). 

- `r pkg("MABOUST")` conducts and simulates the MABOUST design, including making interim decisions to stop a treatment for inferiority or stop the trial early for superiority or equivalency.

- `r pkg("MAMS")` designs multi-arm multi-stage studies with (asymptotically) normal endpoints and known variance. It could be used to determine the boundaries of a multi-arm multi-stage study for a given boundary shape and finds the required number of subjects, as well as simulates multi-arm multi-stage designs and estimates power and expected sample size. 

- `r pkg("MinEDfind")` The nonparametric two-stage Bayesian adaptive design is a novel phase II clinical trial design for finding the minimum effective dose (MinED) in single-agent trials. This design is motivated by the top priority and concern of clinicians when testing a new drug, which is to effectively treat patients and minimize the chance of exposing them to subtherapeutic or overly toxic doses. 

- `r pkg("NCC")` supports the design and analysis of flexible platform trials with non-concurrent controls. Functions for data generation, analysis, visualization and running simulation studies are provided. The implemented analysis methods are described in Bofill Roig et al. (2022) `r doi("10.1186/s12874-022-01683-w")`, Saville et al. (2022) `r doi("10.1177/17407745221112013")` and Schmidli et al. (2014) `r doi("10.1111/biom.12242")`.

- `r pkg("rpact", priority = "core")` Design and analysis of confirmatory adaptive clinical trials with continuous, binary, and survival endpoints according to the methods described in the monograph by Wassmer and Brannath (2016) `r doi("10.1007/978-3-319-32562-0")`. This includes classical group sequential as well as multi-stage adaptive hypotheses tests that are based on the combination testing principle.

- `r pkg("SAME")` allows design of a Bayesian seamless multi-arm biomarker-enriched phase II/III design with the survival endpoint with allowing sample size re-estimation, based on Wason et al. (2015) `r doi("10.1038/bjc.2015.278")`, Yin et al. (2018) `r doi("10.1007/s12561-017-9199-7")` and Yuan et al. (2016) `r doi("10.1002/sim.6971")`.


#### *Bioequivalence*

- `r pkg("adaptIVPT")` contains functions carrying out adaptive procedures using mixed scaling approach to establish bioequivalence for in-vitro permeation test (IVPT) data. Currently, the package provides procedures based on parallel replicate design and balanced data, according to the U.S. Food and Drug Administration's [Draft Guidance on Acyclovir](https://www.accessdata.fda.gov/drugsatfda_docs/psg/PSG_021478.pdf).

- `r pkg("PK")` contains methods to estimate PK parameters using non-compartmental theory and provides facilities to obtain confidence intervals and perform tests for single analysis as well as bioequivalence studies.

- `r pkg("PowerTOST", priority = "core")` is a specialized R package designed to support the planning and evaluation of bioequivalence studies. It offers a comprehensive set of functions for calculating statistical power, sample size, and confidence intervals across a wide range of study designs.

- `r pkg("replicateBE")` Performs comparative bioavailability calculations for Average Bioequivalence with Expanding Limits (ABEL). Implemented are 'Method A' and 'Method B' and the detection of outliers.
    

#### *Dose-Finding*

- `r pkg("BayesianMCPMod")` implements a Bayesian extension of MCPMod based on Fleischer et al. (2022) `r doi("10.1002/pst.2193")`, enabling the systematic inclusion of historical data. It supports simulation, analysis, and evaluation of trial designs across varying dose-response relationships and sample sizes. Users can specify robust mixture priors (e.g., via the Meta-Analytic-Predictive approach), apply weighted model averaging, and assess Minimally Efficacious and Target Doses. Bootstrapped estimates and visualizations of dose-response relationships are also provided.

- `r pkg("bcrm", priority = "core")` This package implements a wide variety of one and two-parameter Bayesian CRM designs. The program can run interactively, allowing the user to enter outcomes after each cohort has been recruited, or via simulation to assess operating characteristics.

- `r pkg("crmPack")` Implements a wide range of model-based dose escalation designs, ranging from classical and modern continual reassessment methods (CRMs) based on dose-limiting toxicity endpoints to dual-endpoint designs taking into account a biomarker/efficacy outcome. The focus is on Bayesian inference, making it very easy to setup a new design with its own JAGS code.
    
- `r pkg("dfcrm", priority = "core")` This package provides functions to run the CRM and TITE-CRM in phase I trials and calibration tools for trial planning purposes.    
    
- `r pkg("DTAT")` Dose Titration Algorithm Tuning (DTAT) is a methodologic framework allowing dose individualization to be conceived as a continuous learning process that begins in early-phase clinical trials and continues throughout drug development, on into clinical practice. This package includes code that researchers may use to reproduce or extend key results of the DTAT research programme, plus tools for trialists to design and simulate a '3+3/PC' dose-finding study.
    
- `r pkg("ewoc")` An implementation of a variety of escalation with overdose control designs introduced by Babb et al. (1998) `r doi("10.1002/(SICI)1097-0258(19980530)17:10%3C1103::AID-SIM793%3E3.0.CO;2-9")`. It calculates the next dose as a clinical trial proceeds as well as performs simulations to obtain operating characteristics.    

- `r pkg("OncoBayes2")` Bayesian Logistic Regression for Oncology Dose-Escalation Trials. It provides flexible functions for Bayesian meta-analytic modeling of the incidence of Dose Limiting Toxicities (DLTs) by dose level, under treatment regimes involving any number of combination partners.

- `r pkg("UnifiedDoseFinding")` includes three dose finding designs: unified phase I design based on Ivanova et al. (2009) `r doi("10.1111/j.1541-0420.2008.01045.x")`), Quasi-CRM/Robust-Quasi-CRM based on Yuan et al. (2007) `r doi("10.1111/j.1541-0420.2006.00666.x")` and Pan et al. (2014) `r doi("10.1371/journal.pone.0098147")`, and generalized BOIN design based on Mu et al. (2018) `r doi("10.1111/rssc.12263")`.  The design goal can include finding the dose associated with a certain target toxicity rate, or the goal can be to find the dose with a certain weighted sum of rates of various toxicity grades, or to find the dose with a certain mean value of a continuous response.

- `r pkg("DoseFinding")` provides functions for the design and analysis of dose-finding experiments (for example pharmaceutical Phase II clinical trials). It provides functions for: multiple contrast tests, fitting non-linear dose-response models, calculating optimal designs and an implementation of the `r pkg("MCPMod")` methodology, but currently only normally distributed homoscedastic endpoints are supported.

- `r pkg("pocrm")` implements functions to implement and simulate the partial order continual reassessment method (PO-CRM) for use in Phase I trials of combinations of agents.

- `r pkg("escalation")` Implements a range of different approaches for dose-finding clinical trials including the continual reassessment method (CRM), the modified TPI (mTPI) design, the Bayesian optimal interval design (BOIN), EffTox and the 3+3 design.

- `r pkg("MCPMod")` This package implements a methodology for the design and analysis of dose-response studies that combines aspects of multiple comparison procedures and modeling approaches based on Bretz et al. (2005) `r doi("10.1111/j.1541-0420.2005.00344.x")`. Please note: The `MCPMod` package will not be further developed, all future development of the MCP-Mod methodology will be done in `r pkg("DoseFinding")`.
    
- `r pkg("TEQR", priority = "core")` The TEQR package contains software to calculate the operating characteristics for the TEQR and the ACT designs.  The TEQR (toxicity equivalence range) design is a toxicity based cumulative cohort design with added safety rules; the ACT (Activity constrained for toxicity) design is also a cumulative cohort design with additional safety rules with the unique feature of this design is that dose is escalated based on lack of activity rather than on lack of toxicity and is de-escalated only if an unacceptable level of toxicity is experienced.
    
- `r pkg("BOIN")` The Bayesian optimal interval (BOIN) design is a novel phase I clinical trial design for finding the maximum tolerated dose (MTD). It can be used to design both single-agent and drug-combination trials. The BOIN design yields an average performance that is comparable to that of the continual reassessment method (CRM, one of the best model-based designs) in terms of selecting the MTD, but has a substantially lower risk of assigning patients to subtherapeutic or overly toxic doses. For tutorial, please check Yan et al. (2020) `r doi("10.18637/jss.v094.i13")`.

- `r pkg("SEARS")` A seamless design that combines phase I dose escalation based on toxicity with phase II dose expansion and dose comparison based on efficacy. A rich set of parameters can be used to explore various real scenarios. It can generate operating characteristics via simulation to examine the design's property.

- `r pkg("MinEDfind")` supports dose determination for upcoming patient cohort in single-agent trials designed to identify the minimum effective dose.
    
- `r pkg("dfmta")` Phase I/II adaptive dose-finding design for single-agent Molecularly Targeted Agent (MTA), according to the paper *Phase I/II Dose-Finding Design for Molecularly Targeted Agent: Plateau Determination using Adaptive Randomization* by Riviere Marie-Karelle et al. (2016) `r doi("10.1177/0962280216631763")`.

- `r pkg("iAdapt")` Simulate and implement early phase two-stage adaptive dose-finding design for binary and quasi-continuous toxicity endpoints. See Chiuzan et al. (2018) `r doi("10.1080/19466315.2018.1462727")` for further reading.


#### *Factorial Designs*

- `r pkg("conf.design")` This small package contains a series of simple tools for constructing and manipulating confounded and fractional factorial designs.

- `r pkg("FrF2")` This package creates regular and non-regular Fractional Factorial designs. Furthermore, analysis tools for Fractional Factorial designs with 2-level factors are offered (main effects and interaction plots for all factors simultaneously, cube plot for looking at the simultaneous effects of three factors, full or half normal plot, alias structure in a more readable format than with the built-in function alias).
    

#### *Group Sequential Designs*
    
- `r pkg("BinGSD")` supports the computation of boundaries and conditional power for single-arm group sequential test with binary endpoint, via either asymptotic or exact test. The package also provides functions to obtain boundary crossing probabilities given the design.   
    
- `r pkg("clinfun", priority = "core")` has functions for both design and analysis of clinical trials. For phase II trials, it has functions to calculate sample size, effect size, and power based on Fisher's exact test, the operating characteristics of a two-stage boundary, Optimal and Minimax 2-stage Phase II designs given by Richard Simon, the exact 1-stage Phase II design and can compute a stopping rule and its operating characteristics for toxicity monitoring based repeated significance testing; and can calculate sample size for Phase III group sequential designs.    
    
- `r pkg("GroupSeq")` computes probabilities related to group sequential designs for normally distributed test statistics. Enables to derive critical boundaries, power, drift, and confidence intervals of such designs. Supports the alpha spending approach by Lan-DeMets (1994) `r doi("10.1002/sim.4780131308")`.

- `r pkg("grpseq")` supports the design of group sequential trials, including non-binding futility analysis at multiple time points, based on Gallo et al. (2014) `r doi("10.1080/10543406.2014.932285")`.

- `r pkg("gscounts")` supports the design and analysis of group sequential designs for negative binomial outcomes, as described by Mütze et al. (2018) `r doi("10.1177/0962280218773115")`.

- `r pkg("gsDesign")` derives group sequential clinical trial designs and describes their properties. Particular focus on time-to-event, binary, and continuous outcomes. Largely based on methods described in the book *Group Sequential Methods with Applications to Clinical Trials* by Jennison et al. (Chapman & Hall/CRC, ISBN 0-8493-0316-8, 2000).

- `r pkg("GSED")` provides function to apply "Group sequential enrichment design incorporating subgroup selection" (GSED) method proposed by Magnusson and Turnbull (2013) `r doi("10.1002/sim.5738")`.

- `r pkg("gsrsb")` implements a gate-keeping procedure to test a primary and a secondary endpoint in a group sequential design with multiple interim looks, including refined secondary boundaries for a gate-keeping test. It supports computations for both standard boundaries and boundaries using error spending functions. See Tamhane et al. (2018) `r doi("10.1111/biom.12732")` for details.

- `r pkg("ldbounds", priority = "core")` uses Lan-DeMets Method for group sequential trial; its functions calculate bounds and probabilities of a group sequential trial.

- `r pkg("lrstat")` enables the design of adaptive group sequential trials, allowing flexibility in sample size adjustments, error spending functions, and the number and timing of interim analyses. It also supports a range of methods for adjusted p-values, including graphical approaches and gatekeeping procedures.

- `r pkg("rpact", priority = "core")` is a comprehensive validated R package for clinical research which enables the design and analysis of confirmatory adaptive group sequential designs with continuous, binary, and survival endpoints.

- `r pkg("SurrogateSeq")` provides functions to implement group sequential procedures that allow for early stopping to declare efficacy using a surrogate marker and the possibility of futility stopping. More details are available in Parast and Bartroff (2024) `r doi("10.1093/biomtc/ujae108")`.


#### *Randomization*

- `r pkg("blockrand", priority = "core")` creates randomizations for block random clinical trials. It can also produce a PDF file of randomization cards.

- `r pkg("experiment", priority = "core")` contains tools for clinical experiments, e.g., a randomization tool, and it provides a few special analysis options for clinical trials.
    
- `r pkg("randomizeR")` This tool enables the user to choose a randomization procedure based on sound scientific criteria. It comprises the generation of randomization sequences as well the assessment of randomization procedures based on carefully selected criteria.


#### *Response Adaptive Randomization*

- `r pkg("BAR")` Bayesian adaptive randomization is also called outcome adaptive randomization, which is increasingly used in clinical trials.

- `r pkg("brada")` provides access to a range of functions for analyzing, applying and visualizing Bayesian response-adaptive trial designs for a binary endpoint. Includes the predictive probability approach and the predictive evidence value designs for binary endpoints.

- `r pkg("carat")` provides functions and command-line user interface to generate allocation sequence by covariate-adaptive randomization for clinical trials. The package currently supports six covariate-adaptive randomization procedures. See Ma et al. (2023) `r doi("10.18637/jss.v107.i02")` for details. 

- `r pkg("CARM")` In randomized controlled trial (RCT), balancing covariate is often one of the most important concern. CARM package provides functions to balance the covariates and generate allocation sequence by covariate-adjusted Adaptive Randomization via Mahalanobis-distance (ARM) for RCT.  For details, please see Yang et al. (2023) `r doi("10.1016/j.csda.2022.107642")`.

- `r pkg("covadap")` Implements seven Covariate-Adaptive Randomization to assign patients to two treatments. Given a set of covariates, the user can generate a single sequence of allocations or replicate the design multiple times by simulating the patients' covariate profiles. See  Antognini et al. (2022) `r doi("10.1007/s00362-022-01381-1")` for details.

- `r pkg("grouprar")` implements group response-adaptive randomization procedures, which also integrates standard non-group response-adaptive randomization methods as specialized instances. It is also uniquely capable of managing complex scenarios, including those with delayed and missing responses, thereby expanding its utility in real-world applications. 

- `r pkg("RABR")` supports Response Adaptive Block Randomization (RABR) design through simulations to evaluate its type I error rate, power and operating characteristics for binary and continuous endpoints. For more details of the proposed method, please refer to Zhan et al. (2021) `r doi("10.1002/sim.9104")`.

- `r pkg("RARfreq")` provides functions and command-line user interface to generate allocation sequence by response-adaptive randomization for clinical trials. The package currently supports two families of frequentist response-adaptive randomization procedures, Doubly Adaptive Biased Coin Design ('DBCD') and Sequential Estimation-adjusted Urn Model ('SEU'), for binary and normal endpoints.


#### *Sample Size and Power Calculations*

- `r pkg("BayesCTDesign")` A set of functions to help clinical trial researchers calculate power and sample size for two-arm Bayesian randomized clinical trials that do or do not incorporate historical control data.  Outcomes considered are Gaussian, Poisson, Bernoulli, Lognormal, Weibull, and Piecewise Exponential. The methods are described in Eggleston et al. (2021) `r doi("10.18637/jss.v100.i21")`. 
    
- `r pkg("clinfun", priority = "core")` provides functions to determine sample sizes, effect sizes, and power based on Fisher’s exact tests, as well as functions to calculate the power of rank tests for animal studies. 
    
- `r pkg("cosa")` Implements bound constrained optimal sample allocation (BCOSA) framework described in Bulus & Dong (2019) for power analysis of multilevel regression discontinuity designs (MRDDs) and multilevel randomized trials (MRTs) with continuous outcomes. Separate tools for statistical power and minimum detectable effect size computations are provided.
    
- `r pkg("longpower", priority = "core")`Compute power and sample size for linear models of longitudinal data. Supported models include mixed-effects models and models fit by generalized least squares and generalized estimating equations. The package is described in Iddi and Donohue (2022) `r doi("10.32614/RJ-2022-022")`. 

- `r pkg("MKpower")` performs power and sample size calculations for various tests, including Welch and Hsu t-tests (via Monte Carlo simulation), Wilcoxon rank sum and signed rank tests, diagnostic test evaluation, single proportion tests, comparison of negative binomial rates, ANCOVA, reference ranges, multiple primary endpoints, and AUC.
    
- `r pkg("lrstat")` performs power and sample size calculation for non-proportional hazards model using the Fleming-Harrington family of weighted log-rank tests.

- `r pkg("pmvalsampsize")` computes the minimum sample size required for the external validation of an existing multivariable prediction model using the criteria proposed by Archer (2020) `r doi("10.1002/sim.8766")` and Riley (2021) `r doi("10.1002/sim.9025")`.
    
- `r pkg("PowerTOST", priority = "core")` contains functions to calculate power and sample size for various study designs used in bioequivalence studies. Power and sample size can be obtained based on different methods, amongst them prominently the TOST procedure (two one-sided t-tests).
    
- `r pkg("PowerUpR")` Includes tools to calculate statistical power, minimum detectable effect size (MDES), MDES difference (MDESD), and minimum required sample size for various multilevel randomized experiments with continuous outcomes. Some of the functions can assist with planning two- and three-level cluster-randomized trials (CRTs) sensitive to multilevel moderation and mediation (2-1-1, 2-2-1, and 3-2-1).
    
- `r pkg("presize")` Bland (2009) `r doi("10.1136/bmj.b3985")` recommended to base study sizes on the width of the confidence interval rather the power of a statistical test. The goal of 'presize' is to provide functions for such precision based sample size calculations. For a given sample size, the functions will return the precision (width of the confidence interval), and vice versa.

- `r pkg("pwr", priority = "core")` Power calculations along the lines of Cohen (1988) `r doi("10.4324/9780203771587")` using in particular the same notations for effect sizes. Examples from the book are given.

- `r pkg("rpact", priority = "core")` provides sample size and power calculations for a range of endpoints, including: means (continuous endpoint), rates (binary endpoint), survival trials with flexible recruitment and survival time options, and count data.
    
- `r pkg("samplesize")` computes sample size for Student's t-test with equal and nonequal variances and for the Wilcoxon-Mann-Whitney test for categorical data with and without ties.
    
- `r pkg("ssanv")` is a set of functions to calculate sample size for two-sample difference in means tests. Does adjustments for either nonadherence or variability that comes from using data to estimate parameters.

- `r pkg("TrialSize", priority = "core")` This package has more than 80 functions from the book *Sample Size Calculations in Clinical Research* by Chow et al. (2007) `r doi("10.1201/9781584889830")`.

    
#### *Simulation*

- `r pkg("adaptDiag")` simulate clinical trials for diagnostic test devices and evaluate the operating characteristics under an adaptive design with futility assessment determined via the posterior predictive probabilities.

- `r pkg("airship")`is an R package that contains an R Shiny App designed to plot simulation results of clinical trials. Its main feature is allowing users to simultaneously investigate the impact of several simulation input dimensions through dynamic filtering of the simulation results. A more detailed description of the core app can be found in Meyer et al. (2023) `r doi("10.1016/j.softx.2023.101347")`.

- `r pkg("asd", priority = "core")` provides functions to simulate adaptive seamless designs that either (i) compare multiple experimental treatments against a single control group, or (ii) evaluate a single experimental treatment versus a control with co-primary analyses in both a predefined subgroup and the full population.

- `r pkg("bcrm", priority = "core")` allows users to Simulate multiple trial scenarios to assess operating characteristics and evaluate different Bayesian CRM designs.
    
- `r pkg("cats")` Given trial-specific design parameters, this package performs multiple trial simulations and exports the results to an Excel file for further analysis.
    
- `r pkg("CohortPlat")` is designed to simulate cohort-based platform trials that evaluate combination therapies involving two active compounds. 
    
- `r pkg("esDesign")` is developed to implement adaptive enrichment designs.  esDesign functions to conduct simulation studies of adaptive enrichment designs include the following strategies:  without early stopping boundary; with sample size re-estimation procedure; with Sample Size Re-estimation Procedure based on Futility and Efficacy Stopping Boundaries for the continuous endpoint.  Simulation studies can also be conducted for these additional designs: Marker Sequential Test design; standard design; sample size re-estimation procedure..   

- `r pkg("ewoc")` performs simulations to obtain operating characteristics for dose escalation overdose control designs for phase I clinical trials such that the probability of overdose is controlled explicitly.   The design was first introduced by Babb et al. (1998) `r doi("10.1002/(SICI)1097-0258(19980530)17:10%3C1103::AID-SIM793%3E3.0.CO;2-9")`. 

- `r pkg("Mediana")` provides a general framework for clinical trial simulations based on the Clinical Scenario Evaluation (CSE) approach. The package supports a broad class of data models (including clinical trials with continuous, binary, survival-type and count-type endpoints as well as multivariate outcomes that are based on combinations of different endpoints), analysis strategies and commonly used evaluation criteria.

- `r pkg("NCC")` is an R package that allows users to simulate platform trials and perform treatment–control comparisons using non-concurrent control data based on Krotka et al. (2023) `r doi("10.1016/j.softx.2023.101437")`. The package supports simulation of complex platform trial designs with continuous or binary endpoints and a flexible number of treatment arms that enter the trial at different time points and accommodates different treatment effects among the arms and includes several patterns for time trends using frequentist approach (e.g., regression model adjusting for time as a fixed effect mixed model adjusting for time as a random factor, and regression splines), the Bayesian time machine a meta-analytic predictive prior separate analysis, and pooled analysis.

- `r pkg("RABR")` conducts simulations of the Response Adaptive Block Randomization (RABR) design to evaluate its type I error rate, power and operating characteristics for binary and continuous endpoints. For more details of the proposed method, please refer to Zhan et al. (2021) `r doi("10.1002/sim.9104")`.

- `r pkg("rpact", priority = "core")` provides simulation tools for means, rates, survival data, and count data, enabling the evaluation of adaptive sample size or event number recalculations based on conditional power, as well as the assessment of treatment selection strategies in multi-arm trials.

- `r pkg("simglm")` Simulates regression models, including both simple regression and generalized linear mixed models with up to three level of nesting. Power simulations that are flexible allowing the specification of missing data, unbalanced designs, and different random error distributions are built into the package.

- `r pkg("UnifiedDoseFinding")` provides the setup and calculations needed to run a dose-finding trial with non-binary endpoints and performs simulations to assess design's operating characteristics under various scenarios.  Three dose finding designs are included in this package: unified phase I design based on Ivanova et al. (2009) `r doi("10.1111/j.1541-0420.2008.01045.x")`, Quasi-CRM/Robust-Quasi-CRM based on Yuan et al. (2007) `r doi("10.1111/j.1541-0420.2006.00666.x")` and Pan et al. (2014) `r doi("10.1371/journal.pone.0098147")`, and generalized BOIN design by Mu et al. (2018) `r doi("10.1111/rssc.12263")`.


### Analysis

#### *General Analysis*

- `r pkg("clintrialX")` fetches clinical trial data from sources like [ClinicalTrials.gov](https://clinicaltrials.gov) and the [Clinical Trials Transformation Initiative - Access to Aggregate Content of ClinicalTrials.gov database](https://aact.ctti-clinicaltrials.org), with support for pagination, bulk downloads, and HTML report generation.

- `r pkg("coin")` offers conditional inference procedures for the general independence problem including two-sample, K-sample (non-parametric ANOVA), correlation, censored, ordered and multivariate problems.
    
- `r pkg("ctrdata")` is a system for querying, retrieving and analyzing protocol- and results-related information on clinical trials from four public registers.
    
- `r pkg("epibasix")` has functions such as `diffdetect`, `n4means` for continuous outcome and `n4props` and functions for matched pairs analysis in randomized trials.
    
- `r pkg("HH")` is support software for *Statistical Analysis and Data Display* (Second Edition, Springer, ISBN 978-1-4939-2121-8, 2015) and (First Edition, Springer, ISBN 0-387-40270-5, 2004) by Richard M. Heiberger and Burt Holland. `ae.dotplot` shows a two-panel display of the most frequently occurring adverse events in the active arm of a clinical study.

- `r pkg("logistf")` Firth's Bias-Reduced Logistic Regression. Firth's method was proposed as ideal solution to the problem of separation in logistic regression.

- `r pkg("multcomp")` covers simultaneous tests and confidence intervals for general linear hypotheses in parametric models, including linear, generalized linear, linear mixed effects, and survival models.

- Base R, especially the `stats` package, has a lot of functionality useful for design and analysis of clinical trials. For example, `chisq.test`, `prop.test`, `binom.test`, `t.test`, `wilcox.test`, `kruskal.test`, `mcnemar.test`, `cor.test`, `power.t.test`, `power.prop.test`, `power.anova.test`, `lm`, `glm`, `nls`, `anova` (and its `lm` and `glm` methods) among many others.
    
- `r pkg("TestDesign")` uses the optimal test design approach by Birnbaum (Addison-Wesley, ISBN 9781593119348, 1968) and van der Linden (2018) `r doi("10.1201/9781315117430")` to construct fixed, adaptive, and parallel tests. Supports the following mixed-integer programming (MIP) solver packages: `Rsymphony`, `gurobi`, `lpSolve`, and `Rglpk`. The `gurobi` package is not available from CRAN; see <https://www.gurobi.com/downloads/>.
    

#### *Longitudinal Data Analysis*

- `r pkg("brms.mmrm")` is a powerful and versatile package for fitting Bayesian regression models.  It leverages 'brms' to run MMRMs, and it supports a simplified interfaced to reduce difficulty and align with the best practices of the life sciences.

- `r pkg("glmmTMB")` fits linear and generalized linear mixed models with various extensions, including zero-inflation.

- `r pkg("lme4")` fits linear and generalized linear mixed-effects models.

- `r pkg("mmrm")` Implements mixed models for repeated measures (MMRM), a popular choice for analyzing longitudinal continuous outcomes in randomized clinical trials and beyond.

- `r pkg("multcomp")` covers simultaneous tests and confidence intervals for general linear hypotheses in parametric models, including linear, generalized linear, linear mixed effects, and survival models.
    
- `r pkg("nlme")` fits and compare Gaussian linear and nonlinear mixed-effects models.
    
    
#### *Meta-Analysis*

This task view focuses on packages relevant to clinical trials. For a more comprehensive list of packages on this topic, please refer to the `r view("Meta-analysis")` task view.

- `r pkg("meta")` is a user-friendly package offering standard meta-analysis methods as described in the book *Meta-Analysis with R* by Schwarzer et al. (2015) `r doi(" 10.1007/978-3-319-21416-0")`, featuring common and random effects models, various plots (e.g., forest, funnel), advanced models (e.g., three-level, GLMM), bias evaluation, meta-regression, cumulative and leave-one-out analysis, and subgroup forest plot summaries.
    
- `r pkg("metafor")` is a comprehensive package for meta-analyses, providing functions to calculate effect sizes, fit various (e.g., fixed-, and random-effects) models, perform moderator/meta-regression analyses, create meta-analytical plots, apply specialized methods (e.g., Mantel-Haenszel method, Peto's method), and fit meta-analytic multivariate/multilevel models accounting for non-independent sampling errors or clustering.
    
- `r pkg("metaLik")` Likelihood inference in meta-analysis and meta-regression models.

- `r pkg("metasens")` is a package for statistical methods to model and adjust for bias in meta-analysis.
    
- `r pkg("netmeta")` provides a comprehensive set of functions for frequentist methods in network meta-analysis, including additive models, analysis of binary data, ranking methods (SUCRA, P-scores), consistency checks, league tables, funnel plots, evidence flow measures, network graphs, treatment ranking diagrams, contribution matrices, meta-regression, and subgroup analysis.

- `r pkg("RBesT")` Tool-set to support Bayesian evidence synthesis. It facilitate the use of historical information in clinical trials. Once relevant historical information has been identified, it supports the derivation of informative priors via the Meta-Analytic-Predictive (MAP) approach and the evaluation of the trial’s operating characteristics.

- `r pkg("rmeta")` has functions for simple fixed and random effects meta-analysis for two-sample comparisons and cumulative meta-analyses. Draws standard summary plots, funnel plots, and computes summaries and tests for association and heterogeneity.
    

#### *Missing Data Imputation*

This task view focuses on packages relevant to clinical trials. For a more comprehensive list of packages on this topic, please refer to the `r view("MissingData")` task view.

- `r pkg("mice")` implements multiple imputation by chained equations using Fully Conditional Specification (FCS) implemented by the MICE algorithm as described in Van Buuren and Groothuis-Oudshoorn (2011) `r doi("10.18637/jss.v045.i03")`.

- `r pkg("rbmi")` implements standard and reference based multiple imputation allowing for the imputation of longitudinal datasets using predefined strategies. The package is described in Gower-Page et al. (2022) `r doi(" 10.21105/joss.04251")`.
    
- `r pkg("remiod")` implements Reference-based multiple imputation of ordinal and binary responses under Bayesian framework, as described in Wang and Liu (2022) `r doi("10.48550/arXiv.2203.02771")`. Methods for missing-not-at-random include Jump-to-Reference (J2R), Copy Reference (CR), and Delta Adjustment which can generate tipping point analysis.
    
    
#### *Survival Analysis*

This task view focuses on packages relevant to clinical trials. For a more comprehensive list of packages on this topic, please refer to the `r view("Survival")` task view.

-  `r pkg("maxcombo")` provides functions for comparing survival curves using the max-combo test at a single timepoint or repeatedly at successive respective timepoints while controlling type I error, as published by Prior (2020) `r doi("10.1177/0962280220931560")`.

-  `r pkg("multcomp")` allows simultaneous inference on general linear hypotheses within a parametric model, such as a Cox proportional hazards model or a parametric survival model.

- `r pkg("nphRCT")` performs a stratified weighted log-rank test in a randomized controlled trial. Tests can be visualized as a difference in average score on the two treatment arms.

- `r pkg("rpsftm")` provides functions to fit a rank preserving structural failure time model to a two-arm clinical trial with survival outcomes.
    
- `r pkg("survival", priority = "core")` contains descriptive statistics, two-sample tests, parametric accelerated failure models, Cox model. Delayed entry (truncation) allowed for all models; interval censoring for parametric models. Case-cohort designs.
  
    
#### *Other Analysis for Specific Designs*  

- `r pkg("clinicalsignificance")` The goal of this package is to provide all necessary tools for analyses of clinical significance in clinical intervention studies. In contrast to statistical significance, which assesses if it is probable that there is a treatment effect, clinical significance can be used to determine if a treatment effect is of practical use or meaningful for patients.
   
- `r pkg("clinsig")` This package contains functions to calculate both parametric and non-parametric versions of the Jacobson-Truax estimates of clinical significance.
  
- `r pkg("MatchIt")` is an R package that selects matched samples of the original treated and control groups with similar covariate distributions.  It can be used to match exactly on covariates, to match on propensity scores, or perform a variety of other matching procedures. The package also implements a series of recommendations offered in Ho et al. (2007) `r doi("10.1093/pan/mpl013")`.

- `r pkg("nppbib")` implements a nonparametric statistical test for rank or score data from partially-balanced incomplete block-design experiments.
    
- `r pkg("speff2trial", priority = "core")` performs estimation and testing of the treatment effect in a 2-group randomized clinical trial with a quantitative or dichotomous endpoint.
    
- `r pkg("ThreeGroups")` This package implements the Maximum Likelihood estimator for three-group designs proposed by Gerber et al. (2010) `r doi("10.1093/pan/mpq008")`.
    

### Monitoring

- `r pkg("accrualPlot")` Tracking accrual in clinical trials is important for trial success. `accrualPlot` provides functions to aid the tracking of accrual and predict when a trial will reach it's intended sample size.
    
- `r pkg("esDesign")` provides tools to monitor adaptive enrichment designs.  Users can calculate the futility and/or efficacy stopping boundaries, the sample size required, calibrate the value of the threshold of the difference between subgroup-specific test statistics.

- `r pkg("monitOS")` Monitoring Overall Survival in Pivotal Trials in Indolent Cancers.

- `r pkg("PwrGSD")` provides tools for the evaluation of interim analysis plans for sequentially monitored trials on a survival endpoint; tools to construct efficacy and futility boundaries, for deriving power of a sequential design at a specified alternative, template for evaluating the performance of candidate plans at a set of time varying alternatives.

- `r pkg("rpact", priority = "core")` enables automatic boundary recalculations during a trial using the alpha spending approach, accommodating both under-running and over-running scenarios.

- `r pkg("SAME")` allows monitoring of a Bayesian seamless multi-arm biomarker-enriched phase II/III design with the survival endpoint with allowing sample size re-estimation, based on Wason et al. (2015) `r doi("10.1038/bjc.2015.278")`, Yin et al. (2018) `r doi("10.1007/s12561-017-9199-7")` and Yuan et al. (2016) `r doi("10.1002/sim.6971")`.

- `r pkg("seqmon")` provides sequential monitoring of clinical trials. It calculates the efficacy and futility boundaries at each look. It allows modifying the design and tracking the design update history.

- `r pkg("tLagInterim")` supports interim monitoring in clinical trials with time-lagged outcomes, It implements inverse and augmented inverse probability weighted estimators for common treatment effect parameters at an interim analysis with time-lagged outcome that may not be available for all enrolled subjects. See Tsiatis and Davidian (2022) `r doi("10.1002/sim.9580")` for details.

### Reporting

- `r pkg("consort")` creates CONSORT diagrams for randomized clinical trials using standardized disposition data, with optional custom text labels for nodes.

- `r pkg("gridify")` is a simple, flexible tool for creating enriched figures and tables by adding surrounding text using predefined or custom layouts. Supports any input convertible to a `grob` (e.g., `ggplot`, `gt`, `flextable`) and is built on R's `grid` graphics system. For details, see Murrell (2018) `r doi("10.1201/9780429422768")`.

- `r pkg("junco")` provides additional tools to enhance `rtables`, `rlistings`, and `tern` for generating tables and listings. It offers general-purpose features like extended statistical analyses, a production-ready RTF exporter, support for spanning headers and risk difference columns, and font-aware auto column width adjustment.

- `r pkg("rlistings")` provides a framework for formatting large datasets, with features tailored for listings commonly included in clinical trial submissions for regulatory review.

- `r pkg("rtables")` provides a framework for defining and applying complex, multi-level tabulations using a hierarchical, tree-like structure. It supports flexible row and column splits, multi-value cells, contextual summaries, and a pipe-friendly interface for layout and computation.

- `r pkg("tern")` creates tables, listings, and graphs (TLG) library for common outputs used in clinical trials.

- `r pkg("tidytlg")` generates tables, listings, and graphs (TLG) using `tidyverse`, supporting both functional workflows and metadata-driven summaries. It can also integrate with the `envsetup` package for environment setup..

### Links
- [Regulatory Compliance and Validation Issues (A Guidance Document for the Use of R in Regulated Clinical Trial Environments)](https://www.R-project.org/doc/R-FDA.pdf)
