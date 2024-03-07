# A script to make the metadata.csv file located in inst/extdata of the package.
library(HubPub)

# kinase annotation
metadata <- HubPub::hub_metadata(
  Title = "JohnsonKinaseAnnotation",
  Description = "Annotation data for each kinase PWM",
  BiocVersion = "3.19",
  Genome = NA_character_,
  SourceType = "XLS/XLSX",
  SourceUrl = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05575-3/MediaObjects/41586_2022_5575_MOESM3_ESM.xlsx",
  SourceVersion = "31.01.2024",
  Species = "Homo sapiens",
  TaxonomyId = as.integer(9606),
  Coordinate_1_based = NA,
  DataProvider = "Nature publishing group",
  Maintainer = "Florian Geier <florian.geier@unibas.ch>",
  RDataClass = "character",
  DispatchClass = "FilePath",
  Location_Prefix = "https://zenodo.org/",
  RDataPath = "records/10628874/files/JohnsonKinaseAnnotation.csv",
  Tags = "Phosphoproteomics:Kinase activity:PWM"
)

HubPub::add_resource(here::here(), metadata)

# kinase PWMs
metadata <- HubPub::hub_metadata(
  Title = "JohnsonKinasePWM",
  Description = "List of kinase position weight matrices (PWMs)",
  BiocVersion = "3.19",
  Genome = NA_character_,
  SourceType = "XLS/XLSX",
  SourceUrl = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05575-3/MediaObjects/41586_2022_5575_MOESM4_ESM.xlsx",
  SourceVersion = "31.01.2024",
  Species = "Homo sapiens",
  TaxonomyId = as.integer(9606),
  Coordinate_1_based = NA,
  DataProvider = "Nature publishing group",
  Maintainer = "Florian Geier <florian.geier@unibas.ch>",
  RDataClass = "character",
  DispatchClass = "FilePath",
  Location_Prefix = "https://zenodo.org/",
  RDataPath = "records/10628874/files/JohnsonKinasePWM.csv.gz",
  Tags = "Phosphoproteomics:Kinase activity:PWM"
)

HubPub::add_resource(here::here(), metadata)

# kinase background scores
metadata <- HubPub::hub_metadata(
  Title = "JohnsonKinaseBackroundScores",
  Description = "Matrix of PWM scores on a set of background phosphosites",
  BiocVersion = "3.19",
  Genome = NA_character_,
  SourceType = "XLS/XLSX",
  SourceUrl = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05575-3/MediaObjects/41586_2022_5575_MOESM5_ESM.xlsx",
  SourceVersion = "31.01.2024",
  Species = "Homo sapiens",
  TaxonomyId = as.integer(9606),
  Coordinate_1_based = NA,
  DataProvider = "Nature publishing group",
  Maintainer = "Florian Geier <florian.geier@unibas.ch>",
  RDataClass = "character",
  DispatchClass = "FilePath",
  Location_Prefix = "https://zenodo.org/",
  RDataPath = "records/10628874/files/JohnsonKinaseBackgroundScores.csv.gz",
  Tags = "Phosphoproteomics:Kinase activity:PWM"
)

HubPub::add_resource(here::here(), metadata)

# kinase background score quantiles
metadata <- HubPub::hub_metadata(
  Title = "JohnsonKinaseBackroundQuantiles",
  Description = "Matrix of PWM score quantiles on a set of background phosphosites",
  BiocVersion = "3.19",
  Genome = NA_character_,
  SourceType = "XLS/XLSX",
  SourceUrl = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05575-3/MediaObjects/41586_2022_5575_MOESM5_ESM.xlsx",
  SourceVersion = "31.01.2024",
  Species = "Homo sapiens",
  TaxonomyId = as.integer(9606),
  Coordinate_1_based = NA,
  DataProvider = "Nature publishing group",
  Maintainer = "Florian Geier <florian.geier@unibas.ch>",
  RDataClass = "character",
  DispatchClass = "FilePath",
  Location_Prefix = "https://zenodo.org/",
  RDataPath = "records/10628874/files/JohnsonKinaseBackgroundQuantiles.csv.gz",
  Tags = "Phosphoproteomics:Kinase activity:PWM"
)

HubPub::add_resource(here::here(), metadata)
