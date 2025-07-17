library(Capr)

# Exposure
h1n1pdmVaccines <- cs(
  descendants(c(40166130, 40166144, 40166607, 40166608, 40166609, 40166611, 40213187, 40225009, 40240135, 42902936, 45776785)),
  name = "H1N1pdm vaccines"
)
h1n1pdmVaccinesCohort <- cohort(
  entry = entry(
    drugExposure(h1n1pdmVaccines, firstOccurrence())
  )
)

# Outcomes
gbsCohort <- PhenotypeLibrary::getPlCohortDefinitionSet(1080)

hypothermia <- cs(
  descendants(435371),
  name = "Hypothermia"
)
inpatient <- cs(
  descendants(262, 9201),
  name = "Inpatient or inpatient ER visit"
)
hypothermiaCohort <- cohort(
  entry = entry(
    conditionOccurrence(
      hypothermia, 
      nestedWithAny(
        atLeast(1, visit(inpatient), duringInterval(eventStarts(-Inf, 0), eventEnds(0, Inf)))
      )
    ),
    primaryCriteriaLimit = "All",
    qualifiedLimit = "All"
  ),
  attrition = attrition(
    "has no events in prior 'clean window'" =  withAll(
      exactly(0, 
              conditionOccurrence(hypothermia), 
              duringInterval(eventStarts(-180, -1))
      )
    ),
    expressionLimit = "All"
  ),
  exit = exit(
    endStrategy = fixedExit(offsetDays = 0L)
  )
)
# writeLines(Capr::makeCohortSet(hypothermiaCohort)$json[1])

# Create cohorts ---------------------------------------------------------------
allCohorts <- rbind(
  makeCohortSet(h1n1pdmVaccinesCohort, hypothermiaCohort),
  gbsCohort
)

source("RealWorldStudy/DatabaseDetails.R")
connection <- DatabaseConnector::connect(connectionDetails)
for (i in 3:nrow(databases)) {
  database <- databases[i, ]
  message(sprintf("Generating cohorts for %s", database$name))
  cohortTableNames <- CohortGenerator::getCohortTableNames(cohortTable = database$cohortTable)
  CohortGenerator::createCohortTables(
    connection = connection,
    cohortDatabaseSchema = database$cohortDatabaseSchema,
    cohortTableNames = cohortTableNames
  )
  CohortGenerator::generateCohortSet(
    connection = connection,
    cdmDatabaseSchema = database$cdmDatabaseSchema,
    cohortDatabaseSchema = database$cohortDatabaseSchema,
    cohortTableNames = cohortTableNames,
    cohortDefinitionSet = allCohorts
  )
}
DatabaseConnector::disconnect(connection)

# Count cohorts ----------------------------------------------------------------
source("RealWorldStudy/DatabaseDetails.R")
connection <- DatabaseConnector::connect(connectionDetails)
counts <- list()
for (i in 1:nrow(databases)) {
  database <- databases[i, ]
  message(sprintf("Counting cohorts for %s", database$name))
  
  sql <- "SELECT COUNT(*) AS cohort_count,
    cohort_definition_id
  FROM @schema.@table
  GROUP BY cohort_definition_id;"
  count <- DatabaseConnector::renderTranslateQuerySql(
    connection = connection, 
    sql = sql, 
    schema = database$cohortDatabaseSchema,
    table = database$cohortTable,
    snakeCaseToCamelCase = TRUE
  )
  count$databaseId <- rep(database$name, nrow(count))
  counts[[i]] <- count
}
DatabaseConnector::disconnect(connection)
counts <- do.call(rbind, counts)
