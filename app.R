#libraries
library(tidyverse)
library(shiny)
library(plotly)
library(shinythemes)
library(DT)
library(shinyWidgets)

#text values
Intro <- "
*By Lucas Zhang*, 2022 Spring


## Introduction

Gulf War illness (GWI) as an enduring multi-symptom syndrome affects more than one-third of the U.S. veterans 
who served during the Gulf War from 1990 to 1991 (White et al., 2016) . The illness is associated with a variety 
of acute and/or chronic symptoms including fatigue, widespread pain, cognitive difficulties, skin rashes and 
gastrointestinal and respiratory disorders (White et al., 2016). Considerable resources have been devoted to 
characterizing the illness, investigating its possible causes and exploring the underlying mechanisms. Potential 
causes for GWI have been discussed including psychological trauma and exposures to chemicals such as chemical 
warfare agents, pyridostigmine bromide (PB), pesticides, oil well fires and depleted uranium (White et al., 2016). 
However, how genetic factors and environmental exposures would jointly affect the risk of GWI has been insufficiently 
discussed, and the detailed molecular and cellular mechanisms remain uncharacterized. 

The interaction effect between a genetic factor and an environmental exposure can be described as the conditional 
dependence of these two factors on the risk for GWI. Specifically, the interaction effect implies that the disease 
risk given different genetic factors varies across different levels of environmental exposures, or equivalently, 
that the disease risk given different histories of environmental exposures varies in populations with different 
genotypes (Ottman, 1996). Therefore, exploring the gene-environmental interactions in the GWI would provide a more 
accurate estimate for the main effect of genetic and environmental risk factors. It may help reveal the biological 
pathways associated with GWI, and thus, possibly lead to preventive advice and novel therapies (Hunter, 2005).

To determine the molecular mechanisms of GWI and improve the well-being of veterans, the US Department of Veteran 
Affairs Cooperative Studies Program launched a study of the Gulf War Era Cohort and Biorepository (GWECB). The GWECB 
cohort comprises 1343 subjects. The veterans' medical records, blood samples and survey data including their 
demographic information, health behavior and environmental exposures have been collected (Gifford et al., 2021). 
Therefore, this dataset allows genetic epidemiological studies which aim to determine the penetrance function of 
GWI and estimate the gene-environmental interaction effect (Radhakrishnan et al., 2021). Specifically, GWI prevalence 
has been reported significantly associated with the interaction between the length of the deployment to the Gulf War 
and a candidate single nucleotide polymorphism (SNP), rs662, in PON1 gene in the GWECB cohort (Vahey et al., 2021).

In the present project, we aim to re-evaluate the gene-environmental interaction revealed in the study of Vahey et 
al. (2021) by comparing different parameters specified in the penetrance function. Specifically, we will develop an 
interactive web application which allows efficiently simulating mock results based on an alternative setting of 
parameters and comparing with the observation in the GWECB cohort. Additionally, we will adjust our model for the 
measurement error in GWI diagnosis to explore how the error would influence our interpretation of the results. 
Therefore, our project may ultimately contribute to our communal understanding of biological and environmental 
mechanisms, and consequently, the developments in novel therapies for the GWI.
"

Methods <- "
## Methods

A penetrance function relates the proportion of veterans identified with GWI to their different genotypes and environmental exposure histories. It can be parameterized as
$$logit\\left[P\\left(AFF=1\\middle|G,E\\right)\\right]=\\beta_0+\\beta_{1}G+\\beta_{2}E+\\beta_{3}G \\times {E}$$,
where $logit\\left(P\\right)=\\log{\\frac{P}{1-P}}$ , $AFF = 1$ indicates the subject having GWI, $G$ represents the three different possible genotypes (AA, Aa, aa) at a 
specific biallelic susceptibility locus, $E$ represents environmental exposure history, $G \\times E$ represents the interaction term between different genotypes and environmental 
exposure conditions and $\\beta_0 - \\beta_3$ are regression coefficients (Qin et al., 2007; Liu et al., 2012; Gauderman et al., 2017). 
Specifically, a penetrance function with genetic factors only at a biallelic susceptibility locus can be simplified as $logit(P(AFF=1|G))= \\beta_0 + \\beta_{1} G$. 
Assuming Hardy-Weinberg equilibrium, the penetrance function is associated with the prevalence in the overall population as   
$$K = p ^ 2 f_{AA} + 2p \\left( 1 - p \\right) f_{AG} + \\left( 1 - p \\right) ^ 2 f_{GG}$$,
where $K$ is the overall prevalence of GWI, $f_{AA}$, $f_{AG}$, and $f_{GG}$ are the penetrance function given genotype AA, AG and GG, respectively, and $p$ is the 
frequency of allele A.

We developed an application which allows interactively altering the parameters in the regression models, 
simulating new observation and updating figures using the `shiny`, `ggplot2` and `plotly` (Wickham, 2016; 
Sievert, 2020; Chang et al., 2021) packages in the R software (R Core Team, 2021) (version 4.1.2). The initial 
values for the coefficients were estimated from Vahey et al. (Vahey et al., 2021). While a detailed description 
can be accessed from the original publication, a brief summary of GWI diagnosis, genotyping and statistical 
selection of SNPs is summarized below. Subjects in the GWECB cohort were diagnosed for GWI according to 
different criteria, including CDC, CDC Severe, Kansas full case definition, and Kansas symptom (Gifford et 
al., 2021), and the CDC Severe Case Definition were used as the primary outcome. Candidate genes including 
*ACHE*, *BCHE*, *PON1*, and *SOD1* were selected according to previous studies in human or animal models, 
and SNPs in or near these genes were investigated. The associations between the SNPs and environmental 
exposure and their interactions were adjusted for confounders of age and sex and multiple testing using 
a Bonferroni's correction. 

We postulated two measurements errors associated with the diagnosis of GWI. We assumed a different 
participation rate between GWI patients and healthy veterans in the GWECB cohort; namely, the 
patients with GWI phenotype may become more likely than healthy population to be enrolled in 
the study. In addition, we assumed a misclassification rate of GWI diagnosis, including false 
negative and false positive cases. A false negative case occurs when a patient with GWI is 
misclassified as healthy person, while a false positive occurs when a healthy person is 
misclassified as a GWI patient. The application accepted these parameters as inputs from 
the users and subsequently, adjusted the crude probability of overserving GWI cases 
$P\\left(AFF=1|G,E\\right)$ for the different participation rate and misclassification rate and 
re-calculated the results.

Our application has been deployed on [https://lucaszhang.shinyapps.io/Masters_Project/](https://lucaszhang.shinyapps.io/Masters_Project/)
while the source code is available at [https://github.com/Lujun995/Masters_Project](https://github.com/Lujun995/Masters_Project). 

"

Results1 <- "
## Modules

#### 1. The relationship between prevalence and a penetrance function without environmental exposures 

The GWI have a very high prevanlance estimated at ~ 15%-80%. Here we want to explore the parameter 
space of different parameters in the genetic model, and therefore, to identify possible combinations 
of allele frequency and penetrance functions. As for the following genetic model for a specific biallelic 
susceptibility locus of 
$$K = p^2f_{AA}+2p(1-p)f_{AG}+(1-p)^2f_{GG}$$,
where $p$ indicates the allele frequency of *A*, $p^2$, $2p(1-p)$ and $(1-p)^2$ are frequencies of 
*AA*, *AG* and  *GG* genotypes, and $f_{AA}$, $f_{AG}$ and $f_{GG}$ are their penetrance functions, 
respectively.

Suppose the allele G would increase the risk for the GWI while allele A is the wild type. An additive 
genetic model assumes $f_{AG} = f_{AA} + f_{GG}$, while a recessive genetic model assumes 
$f_{AG} = f_{AA} > f_{GG}$. 

Because of a high GWI prevalence in the Gulf War veterans, it is required that the penetrance of AA 
and AG genotypes are not 0% and/or that the disease allele frequency of G is relatively high in all genetic models.
"
Results2.1 <- "
#### 2. The influence of the marginal and interaction effects of genetic and environmental 
factors on the GWI prevalence

Considering the following genetic model for Gene-Environment interactions:

$$logit(Y)=\\beta_0 + \\beta_1 \\cdot G + \\beta_2 \\cdot E + \\beta_3 \\cdot G \\cdot E$$

, where $G$ is the genetic factors, including the AA, AG and GG genotypes, $E$ is the environmental exposure, ranging from 0-3. 0: No exposure, 1: less than 7 days, 2: 7-30 days, 3: more than 30 days.

Observations from Vahey et al. (2021):
"

Results2.2 <-"
The maximum likelihood estimates are:
"

Results2.3 <- "
Remarkably, the gene x environmental interaction coefficients caused different tendencies of GWI frequency over 
the exposure history in different genotypes. Therefore, this parameter would be attributed as different disease 
risk given different genetic factors varies across different levels of environmental exposures.
"

Results3.1 <- "
#### 3. The bias in the GWI frequency caused by a differentiated participation rate 

We assumed a different participation rate between GWI patients and healthy veterans in the GWECB cohort; 
namely, the patients with GWI phenotype may become more likely than healthy population to be enrolled in the study.
The shiny app will adjust the crude probability of overserving GWI cases $P\\left(AFF=1|G,E\\right)$ for the 
different participation rate and re-calculate the results.
"

Results3.2 <- "
Our results showed that the observed GWI frequency would increase if veterans with GWI were more likely to 
participate; however, this type of bias would not change the overall tendency observed. Notably, this bias 
would become more significant when participation rate is relatively low. For example, when the participation 
rate is ~5% among healthy population, a 5% increase in the rate among GWI population would drastically increase 
the observed GWI frequency.
"

Results4.1 <- "
#### 4. The bias in the GWI frequency caused by a diagnostic error of GWI

We assumed a misclassification rate of GWI diagnosis, including false negative and false positive cases. 
A false negative case occurs when a patient with GWI is misclassified as healthy person, while a false 
positive occurs when a healthy person is misclassified as a GWI patient. This application accepts these 
parameters as inputs from the users and subsequently, adjusts the crude probability of overserving GWI 
cases $P\\left(AFF=1|G,E\\right)$ for the error and re-calculates the results.
"

Results4.2 <- "
Our application examined the bias in the observed GWI frequencies introduced by errors in 
GWI diagnosis (Figure 4). A decreased true negative rate, or namely, an increased false positive 
rate significantly increased the observed GWI prevalence. Worse, it altered the overall tendencies 
and decreased the observed marginal and interaction effects of genetic and environmental factors. A 
decrease in true positive rate would decrease the observed GWI prevalence, and it would still diminish 
the observed marginal and interaction effects. As the GWI is associated with a variety of acute and/or 
chronic symptoms (White et al., 2016) and its diagnosis relies excessively on individual answers to questionnaires 
(Gifford et al., 2021), a notable diagnostic error of GWI will likely be caused. This error can introduce a 
systematic bias to the results; therefore, an in-depth investigation will be critically required to quantify 
the misclassification rate.
"

Concl <- "
## Concluding Remarks

We found that a relative 
high frequency of disease allele G and/or a non-zero penetrance of AA genotype would be necessary to induce the 
high GWI prevalence, that gene x environmental interactions were derived from different GWI risk levels of 
environmental exposures given different genetic factors, and that measurement errors such as a differentiated 
participation rate and a high diagnostic error rate would contribute to the high GWI prevalence while change the 
observed patterns differently. Our application and results identified a non-zero penetrance of healthy genotype 
and measurement errors in GWI diagnosis would be potential causes for the high GWI prevalence, and therefore, 
improved our understanding of GWI genetic models.
"

Refs <- "
## References

Chang, W., Cheng, J. & Allaire, J.J., et al. (2021). shiny: Web Application Framework for R. 

Gauderman, W.J., Mukherjee, B. & Aschard, H., et al. (2017). Update on the State of the Science for Analytical Methods for Gene-Environment Interactions. Am. J. Epidemiol., 186, 762-770.

Gifford, E.J., Vahey, J. & Hauser, E.R., et al. (2021). Gulf War illness in the Gulf War Era Cohort and Biorepository: The Kansas and Centers for Disease Control definitions. Life Sci., 278, 119454.

Hunter, D.J. (2005). Gene-environment interactions in human diseases. Nat. Rev. Genet., 6, 287-298.

Liu, C., Maity, A. & Lin, X., et al. (2012). Design and analysis issues in gene and environment studies. Environ. Health-Glob., 11, 93-93.

Ottman, R. (1996). Gene-Environment Interaction: Definitions and Study Designs. Prev. Med., 25, 764-770.

Qin, X., Schmidt, S. & Martin, E., et al. (2007). Visualizing genotype x phenotype relationships in the GAW15 simulated data. BMC Proc, 1 Suppl 1, S132.

R Core Team (2021). R: A Language and Environment for Statistical Computing. Vienna, Austria.

Radhakrishnan, K., Hauser, E.R. & Polimanti, R., et al. (2021). Genomics of Gulf War Illness in U.S. Veterans Who Served during the 1990-1991 Persian Gulf War: Methods and Rationale for Veterans Affairs Cooperative Study #2006. Brain Sciences, 11, 845.

Sievert, C. (2020). Interactive Web-Based Data Visualization with R, plotly, and shiny. Chapman and Hall/CRC.

Vahey, J., Gifford, E.J. & Sims, K.J., et al. (2021). Gene-Toxicant Interactions in Gulf War Illness: Differential Effects of the PON1 Genotype. Brain Sciences, 11, 1558.

White, R.F., Steele, L. & O'Callaghan, J.P., et al. (2016). Recent research on Gulf War illness and other health problems in veterans of the 1991 Gulf War: Effects of toxicant exposures during deployment. Cortex, 74, 449-475.

Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York.
"

p = seq(0, 1, 0.01)

G = rep(c("AA","AG","GG"), each = 4)
E = rep(0:3, times= 3)
Freq = c(0.20, 0.23, 0.25, 0.23, 
         0.20, 0.26, 0.26, 0.40,
         0.20, 0.40, 0.35, 0.46)
OBS <- data.frame(G = G, E = E, Freq = Freq, 
                  logit_Freq = round(log(Freq/(1-Freq)), digits = 3))
OBS.model <- lm(logit_Freq ~ G + E + G * E, data = OBS)
EST <- as.data.frame(t(data.frame(Coefficients = round(OBS.model$coefficients,3))))
names(EST) <- c("Intercept", "G_AG", "G_GG", "E", "G_AG:E", "G_GG:E")
coef_i<-OBS.model$coefficients

# Define UI for application that draws a histogram
ui <- fixedPage(
    theme = "shinythemes/css/sandstone.min.css",
    title = "Maters' Project",
    withMathJax(),
    # section below allows in-line LaTeX via $ in mathjax.
    tags$div(HTML("<script type='text/x-mathjax-config' >
            MathJax.Hub.Config({
            tex2jax: {inlineMath: [['$','$']]}
            });
            </script >
            ")),
    titlePanel("An Interactive Application Evaluating Model Parameters for 
               the Gene-Environment Interaction Studies in the Gulf War Illness"),
    markdown(Intro),
    br(),
    
    markdown(Methods),
    br(),
    
    markdown(Results1),
    br(),
    sidebarLayout(
      sidebarPanel(
        shiny::sliderInput("f_AA", label = "Penetrance of AA",
                           min = 0, max = 1, value = 0,
                           step = 0.01),
        shiny::sliderInput("f_AG", label = "Penetrance of AG",
                           min = 0, max = 1, value = 0.5,
                           step = 0.01),
        shiny::sliderInput("f_GG", label = "Penetrance of GG",
                           min = 0, max = 1, value = 1,
                           step = 0.01)
      ),
      mainPanel(
        plotlyOutput("Prev_G")
      )
    ),
    br(),
    
    markdown(Results2.1),
    DTOutput('table_obs', width = "50%"),
    br(),
    markdown(Results2.2),
    DTOutput('table_est', width = "100%"),
    br(),
    sidebarLayout(
      sidebarPanel(
        shiny::sliderInput("beta0", label = "Intercept",
                           min = -3, max = 0, value = coef_i[1],
                           step = 0.01),
        shiny::sliderInput("beta1_1", label = "G_AG",
                           min = -0.3, max = 0.3, value = coef_i[2],
                           step = 0.002),
        shiny::sliderInput("beta1_2", label = "G_GG",
                           min = -0.3, max = 0.3, value = coef_i[3],
                           step = 0.002),
        shiny::sliderInput("beta2", label = "E",
                           min = 0, max = 0.3, value = coef_i[4],
                           step = 0.001),
        shiny::sliderInput("beta3_1", label = "G_AG:E",
                           min = 0, max = 0.5, value = coef_i[5],
                           step = 0.002),
        shiny::sliderInput("beta3_2", label = "G_GG:E",
                           min = 0, max = 0.5, value = coef_i[6],
                           step = 0.002)
      ),
      mainPanel(
        plotlyOutput("sim1")
      )
    ),
    br(),
    markdown(Results2.3),
    br(),
    
    markdown(Results3.1),
    br(),
    sidebarLayout(
      sidebarPanel(
        shiny::sliderInput("beta0_new", label = "New Intercept",
                           min = -3, max = 0, value = coef_i[1],
                           step = 0.01),
        shiny::sliderInput("alpha1", 
                           label = "Participation Rate of GWI Population",
                           min = 0, max = 1, value = 0.10,
                           step = 0.01),
        shiny::sliderInput("alpha2", 
                           label = "Participation Rate of Healthy Population",
                           min = 0, max = 1, value = 0.05,
                           step = 0.01)
      ),
      mainPanel(
        plotlyOutput("sim2")
      )
    ),
    br(),
    markdown(Results3.2),
    br(),
    
    markdown(Results4.1),
    br(),
    sidebarLayout(
      sidebarPanel(
        shiny::sliderInput("beta0_new2", label = "New Intercept",
                           min = -3, max = 0, value = coef_i[1],
                           step = 0.01),
        shiny::sliderInput("sensitivity", 
                           label = "True Positive Rate (% Diagonosed GWI among all GWI Patients)",
                           min = 0, max = 100, value = 100,
                           step = 1),
        shiny::sliderInput("specificity", 
                           label = "True Negative Rate (% Considered GWI-free among all healthy population)",
                           min = 0, max = 100, value = 100,
                           step = 1)
      ),
      mainPanel(
        plotlyOutput("sim3")
      )
    ),
    markdown(Results4.2),
    br(),
    
    markdown(Concl),
    br(),
    
    markdown(Refs),
    br(),
    hr()
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    output$Prev_G <- plotly::renderPlotly({
      K = p^2*(input$f_AA) + 2*p*(1-p)*(input$f_AG) +
        (1-p)^2*(input$f_GG)
      ggplotly(qplot(x= p, y=K) + theme_bw() + 
                 xlab("A Allele Frequency p") + 
                 ylab("GWI Prevalence in the Overall Population"))
    })
    output$table_obs <- renderDT(OBS)
    output$table_est <- renderDT(EST)
    output$sim1 <- renderPlotly({
      Geno_risk <- 
        data.frame(Genotype = factor(c("AA", "AG", "GG")), 
                   risk = c(0, input$beta1_1, input$beta1_2))
      p <- ggplot(Geno_risk, aes(Genotype, risk)) + geom_point()
      ggplotly(p)
      
      OBS.model.new <- OBS.model
      OBS.model.new$coefficients <- c(input$beta0, input$beta1_1, input$beta1_2,
                                      input$beta2, input$beta3_1, input$beta3_2)
      
      logit_Freq_sim = predict(obj = OBS.model.new, type = "response", newdata = OBS)
      Freq_sim = exp(logit_Freq_sim) / (1+exp(logit_Freq_sim))
      sim <- data.frame(Genotype = factor(G),
                        Exposure = factor(E),
                        Yes = Freq_sim, No = 1- Freq_sim)
      levels(sim$Genotype) = c("AA","GA","GG")
      levels(sim$Exposure) = c("None", "1-6 days", 
                               "7-30 days", "31+ days")
      sim2 <- reshape2::melt(sim, id.vars = c(1,2), value.name="Frequency")
      sim2$variable <- factor(sim2$variable, levels = c("No","Yes"))
      p <- ggplot(sim2, aes(x = Exposure, y = Frequency, color = variable, fill = variable)) + 
        geom_bar(stat="identity")+facet_grid(.~Genotype) + theme_bw() + 
        scale_color_brewer(palette = "Set1") + 
        scale_fill_brewer(palette = "Set1") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      ggplotly(p)
    })
    output$sim2 <- renderPlotly({
      Geno_risk <- 
        data.frame(Genotype = factor(c("AA", "AG", "GG")), 
                   risk = c(0, input$beta1_1, input$beta1_2))
      OBS.model.new <- OBS.model
      OBS.model.new$coefficients <- c(input$beta0_new, input$beta1_1, input$beta1_2,
                                      input$beta2, input$beta3_1, input$beta3_2)
      
      logit_Freq_sim = predict(obj = OBS.model.new, type = "response", newdata = OBS)
      Freq_sim = exp(logit_Freq_sim) / (1+exp(logit_Freq_sim))
      sim <- data.frame(Genotype = factor(G),
                        Exposure = factor(E),
                        Yes = Freq_sim, 
                        No = (1- Freq_sim) )
      sim$Yes <- sim$Yes * input$alpha1
      sim$No <- sim$No * input$alpha2
      sim$total <- sim$Yes + sim$No
      sim$Yes <- sim$Yes / sim$total
      sim$No <- sim$No / sim$total
      sim <- sim[, 1:4]
      levels(sim$Genotype) = c("AA","GA","GG")
      levels(sim$Exposure) = c("None", "1-6 days", 
                               "7-30 days", "31+ days")
      sim2 <- reshape2::melt(sim, id.vars = c(1,2), value.name="Frequency")
      sim2$variable <- factor(sim2$variable, levels = c("No","Yes"))
      p <- ggplot(sim2, aes(x = Exposure, y = Frequency, color = variable, fill = variable)) + 
        geom_bar(stat="identity")+facet_grid(.~Genotype) + theme_bw() +
        scale_color_brewer(palette = "Set1") + 
        scale_fill_brewer(palette = "Set1") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      ggplotly(p)
    })
    output$sim3 <- renderPlotly({
      Geno_risk <- 
        data.frame(Genotype = factor(c("AA", "AG", "GG")), 
                   risk = c(0, input$beta1_1, input$beta1_2))
      OBS.model.new <- OBS.model
      OBS.model.new$coefficients <- c(input$beta0_new2, input$beta1_1, input$beta1_2,
                                      input$beta2, input$beta3_1, input$beta3_2)
      
      logit_Freq_sim = predict(obj = OBS.model.new, type = "response", newdata = OBS)
      Freq_sim = exp(logit_Freq_sim) / (1+exp(logit_Freq_sim))
      sim <- data.frame(Genotype = factor(G),
                        Exposure = factor(E),
                        Yes = Freq_sim, 
                        No = (1- Freq_sim) )
      sim$Diagnosed_Yes <- sim$Yes * input$sensitivity/100 + sim$No * (1-input$specificity/100)
      sim$Diagnosed_No <- sim$Yes * (1-input$sensitivity/100) + sim$No * input$specificity/100
      sim$total <- sim$Diagnosed_Yes + sim$Diagnosed_No
      sim$Diagnosed_Yes <- sim$Diagnosed_Yes / sim$total
      sim$Diagnosed_No <- sim$Diagnosed_No / sim$total
      sim <- sim[, c("Genotype", "Exposure", "Diagnosed_Yes", "Diagnosed_No")]
      levels(sim$Genotype) = c("AA","GA","GG")
      levels(sim$Exposure) = c("None", "1-6 days", 
                               "7-30 days", "31+ days")
      sim2 <- reshape2::melt(sim, id.vars = c(1,2), value.name="Frequency")
      sim2$variable <- factor(sim2$variable, levels = c("Diagnosed_No","Diagnosed_Yes"))
      p <- ggplot(sim2, aes(x = Exposure, y = Frequency, color = variable, fill = variable)) + 
        geom_bar(stat="identity")+facet_grid(.~Genotype) + theme_bw() +
        scale_color_brewer(palette = "Set1") + 
        scale_fill_brewer(palette = "Set1") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      ggplotly(p)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
