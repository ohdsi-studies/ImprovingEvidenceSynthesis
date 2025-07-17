# Code for setting the connection details and schemas for the various databases,
# as well as some folders in the local file system.
library(DatabaseConnector)
library(dplyr)

options(andromedaTempFolder = "e:/andromedaTemp")
rootFolder <- "e:/EvidenceSynthStudy"

connectionDetails <- createConnectionDetails(
  dbms = "spark",
  connectionString = keyring::key_get("databricksConnectionString"),
  user = "token",
  password = keyring::key_get("databricksToken")
)
options(sqlRenderTempEmulationSchema = "scratch.scratch_mschuemi")

# Note: These are the only databases with data for our problem:
databases <- tibble(
  name = c("CCAE",
           "MDCD",
           "MDCR",
           "OptumDoD",
           "OptumEhr"),
  cdmDatabaseSchema = c("merative_ccae.cdm_merative_ccae_v3046",
                        "merative_mdcd.cdm_merative_mdcd_v3038",
                        "merative_mdcr.cdm_merative_mdcr_v3045",
                        "optum_extended_dod.cdm_optum_extended_dod_v3039",
                        "optum_ehr.cdm_optum_ehr_v3037")
) |>
  mutate(cohortTable = paste("es_method_cohort", name, sep = "_"),
         cohortDatabaseSchema = "scratch.scratch_mschuemi",
         folder = file.path(rootFolder, name))


