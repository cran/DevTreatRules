## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "##"
)
def.chunk.hook  <- knitr::knit_hooks$get("chunk")
knitr::knit_hooks$set(chunk = function(x, options) {
  x <- def.chunk.hook(x, options)
  ifelse(options$size != "normalsize", paste0("\\", options$size,"\n\n", x, "\n\n \\normalsize"), x)
})
options(digits=3)

## ------------------------------------------------------------------------
library(DevTreatRules)
head(obsStudyGeneExpressions)

## ----split_data----------------------------------------------------------
set.seed(123)
example.split <- SplitData(data=obsStudyGeneExpressions, n.sets=3, split.proportions=c(0.5, 0.25, 0.25))
table(example.split$partition)

## ----make_subsets, message=FALSE-----------------------------------------
library(dplyr)
development.data <- example.split %>% filter(partition == "development")
validation.data <-  example.split %>% filter(partition == "validation")
evaluation.data <-  example.split %>% filter(partition == "evaluation")

## ------------------------------------------------------------------------
names.influencing.treatment=c("prognosis", "clinic", "age")

## ------------------------------------------------------------------------
names.influencing.rule=c("age", paste0("gene_", 1:10))

## ----build_one_rule------------------------------------------------------
one.rule <- BuildRule(development.data=development.data,
                      study.design="observational",
                      prediction.approach="split.regression",
                      name.outcome="no_relapse",
                      type.outcome="binary",
                      desirable.outcome=TRUE,
                      name.treatment="intervention",
                      names.influencing.treatment=c("prognosis", "clinic", "age"),
                      names.influencing.rule=c("age", paste0("gene_", 1:10)),
                      propensity.method="logistic.regression",
                      rule.method="glm.regression")

## ------------------------------------------------------------------------
vec.approaches=c("OWL", "split.regression", "OWL.framework", "direct.interactions")
vec.rule.methods=c("glm.regression", "lasso", "ridge")

## ------------------------------------------------------------------------
vec.propensity.methods="logistic.regression"

## ----model_selection_on_validation, warning=FALSE------------------------
set.seed(123)
model.selection <- CompareRulesOnValidation(development.data=development.data,
                validation.data=validation.data,
                study.design.development="observational",
                vec.approaches=c("split.regression", "OWL.framework", "direct.interactions"),
                vec.rule.methods=c("glm.regression", "lasso", "ridge"),
                vec.propensity.methods="logistic.regression",
                name.outcome.development="no_relapse",
                type.outcome.development="binary",
                name.treatment.development="intervention",
                names.influencing.treatment.development=c("prognosis", "clinic", "age"),
                names.influencing.rule.development=c("age", paste0("gene_", 1:10)),
                desirable.outcome.development=TRUE)

## ----look_at_validation_split, size="footnotesize"-----------------------
model.selection$list.summaries[["split.regression"]] 

## ----look_at_validation_OWL, size="footnotesize"-------------------------
model.selection$list.summaries[["OWL.framework"]] 

## ----look_at_validation_DI, size="footnotesize"--------------------------
model.selection$list.summaries[["direct.interactions"]] 

## ----selected_split_on_validation----------------------------------------
model.selection$list.summaries$split.regression["propensity_logistic.regression_rule_glm.regression", ]

## ----selected_OWL_on_validation------------------------------------------
model.selection$list.summaries$OWL.framework["propensity_logistic.regression_rule_glm.regression", ]

## ----selected_DI_on_validation-------------------------------------------
model.selection$list.summaries$direct.interactions["propensity_logistic.regression_rule_lasso", ]

## ----rule_split----------------------------------------------------------
rule.split <- model.selection$list.rules$split.regression[["propensity_logistic.regression_rule_glm.regression"]]
coef(rule.split$rule.object.control)
coef(rule.split$rule.object.treat)

## ----rule_OWL------------------------------------------------------------
rule.OWL <- model.selection$list.rules$OWL.framework[["propensity_logistic.regression_rule_glm.regression"]]
coef(rule.OWL$rule.object)

## ----rule_DI-------------------------------------------------------------
rule.DI <- model.selection$list.rules$direct.interactions[["propensity_logistic.regression_rule_lasso"]]
coef(rule.DI$rule.object)

## ----split_on_eval, warning=FALSE----------------------------------------
set.seed(123)
split.eval <- EvaluateRule(evaluation.data=evaluation.data,
                           BuildRule.object=rule.split,
                           study.design="observational",
                           name.outcome="no_relapse",
                           type.outcome="binary",
                           desirable.outcome=TRUE,
                           name.treatment="intervention",
                           names.influencing.treatment=c("prognosis", "clinic", "age"),
                           names.influencing.rule=c("age", paste0("gene_", 1:10)),
                           propensity.method="logistic.regression",
                           bootstrap.CI=FALSE)
split.eval$summaries

## ----OWL_on_eval, warning=FALSE------------------------------------------
set.seed(123)
OWL.eval <- EvaluateRule(evaluation.data=evaluation.data,
                              BuildRule.object=rule.OWL,
                              study.design="observational",
                              name.outcome="no_relapse",
                              type.outcome="binary",
                              desirable.outcome=TRUE,
                              name.treatment="intervention",
                              names.influencing.treatment=c("prognosis", "clinic", "age"),
                              names.influencing.rule=c("age", paste0("gene_", 1:10)),
                              propensity.method="logistic.regression",
                              bootstrap.CI=FALSE)
OWL.eval$summaries

## ----DI_on_eval, warning=FALSE-------------------------------------------
set.seed(123)
DI.eval <- EvaluateRule(evaluation.data=evaluation.data,
                              BuildRule.object=rule.DI,
                              study.design="observational",
                              name.outcome="no_relapse",
                              type.outcome="binary",
                              desirable.outcome=TRUE,
                              name.treatment="intervention",
                              names.influencing.treatment=c("prognosis", "clinic", "age"),
                              names.influencing.rule=c("age", paste0("gene_", 1:10)),
                              propensity.method="logistic.regression",
                              bootstrap.CI=FALSE)
DI.eval$summaries

## ------------------------------------------------------------------------
bootstrap.CI=TRUE
booststrap.CI.replications=1000 # or any other number of replications

