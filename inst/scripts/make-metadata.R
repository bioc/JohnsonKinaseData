# A script to make the metadata.csv file located in inst/extdata of the package.
library(HubPub)
library(here)
   
here::i_am('inst/scripts/make-metadata.R')

metadata_path <- file.path(here::here('inst','extdata'), "metadata.csv")
writeLines(paste(names(hub_metadata()), collapse = ","), 
           metadata_path)

# Ser/Thr kinase data published in Johnson et al. 2023

# Ser/Thr kinase annotation
metadata <- HubPub::hub_metadata(
  Title = "JohnsonKinaseAnnotation",
  Description = "Annotation data for each ser/thr kinase PWM published in Johnson et al. 2023",
  BiocVersion = "3.20",
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
  Tags = "Proteomics:Phosohosite:Serine/Threonine:Kinase activity:Position weight matrix:PWM"
)

HubPub::add_resource(here::here(), metadata)

# Ser/Thr kinase PWMs
metadata <- HubPub::hub_metadata(
  Title = "JohnsonKinasePWM",
  Description = "List of Ser/Thr kinase position weight matrices (PWMs) published in Johnson et al. 2023",
  BiocVersion = "3.20",
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
  Tags = "Proteomics:Phosohosite:Serine/Threonine:Kinase activity:Position weight matrix:PWM"
)

HubPub::add_resource(here::here(), metadata)

# Ser/Thr kinase background scores
metadata <- HubPub::hub_metadata(
  Title = "JohnsonKinaseBackgroundScores",
  Description = "Matrix of PWM scores on a set of background phosphosites published in Johnson et al. 2023",
  BiocVersion = "3.20",
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
  Tags = "Proteomics:Phosohosite:Serine/Threonine:Kinase activity:Position weight matrix:PWM"
)

HubPub::add_resource(here::here(), metadata)

# Ser/Thr kinase background score quantiles
metadata <- HubPub::hub_metadata(
  Title = "JohnsonKinaseBackgroundQuantiles",
  Description = "Matrix of PWM score quantiles on a set of background phosphosites published in Johnson et al. 2023",
  BiocVersion = "3.20",
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
  Tags = "Proteomics:Phosohosite:Serine/Threonine:Kinase activity:Position weight matrix:PWM"
)

HubPub::add_resource(here::here(), metadata)

## Tyrosine kinase data published in Yaron-Barir et al. 2024

# Tyr kinase annotation
metadata <- HubPub::hub_metadata(
    Title = "TyrosineKinaseAnnotation",
    Description = "Annotation data for each tyrosine kinase PWM published in Yaron-Barir et al. 2024",
    BiocVersion = "3.20",
    Genome = NA_character_,
    SourceType = "XLS/XLSX",
    SourceUrl = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07407-y/MediaObjects/41586_2024_7407_MOESM3_ESM.xlsx",
    SourceVersion = "16.09.2024",
    Species = "Homo sapiens",
    TaxonomyId = as.integer(9606),
    Coordinate_1_based = NA,
    DataProvider = "Nature publishing group",
    Maintainer = "Florian Geier <florian.geier@unibas.ch>",
    RDataClass = "character",
    DispatchClass = "FilePath",
    Location_Prefix = "https://zenodo.org/",
    RDataPath = "records/13831966/files/TyrosineKinaseAnnotation.csv",
    Tags = "Proteomics:Phosohosite:Tyrosine:Kinase activity:Position weight matrix:PWM"
)

HubPub::add_resource(here::here(), metadata)

# Tyr kinase PWMs
metadata <- HubPub::hub_metadata(
    Title = "TyrosineKinasePWM",
    Description = "List of tyrosine kinase position weight matrices (PWMs) published in Yaron-Barir et al. 2024",
    BiocVersion = "3.20",
    Genome = NA_character_,
    SourceType = "XLS/XLSX",
    SourceUrl = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07407-y/MediaObjects/41586_2024_7407_MOESM4_ESM.xlsx",
    SourceVersion = "16.09.2024",
    Species = "Homo sapiens",
    TaxonomyId = as.integer(9606),
    Coordinate_1_based = NA,
    DataProvider = "Nature publishing group",
    Maintainer = "Florian Geier <florian.geier@unibas.ch>",
    RDataClass = "character",
    DispatchClass = "FilePath",
    Location_Prefix = "https://zenodo.org/",
    RDataPath = "records/13831966/files/TyrosineKinasePWM.csv.gz",
    Tags = "Proteomics:Phosohosite:Tyrosine:Kinase activity:Position weight matrix:PWM"
)

HubPub::add_resource(here::here(), metadata)

# Tyr kinase background scores
metadata <- HubPub::hub_metadata(
    Title = "TyrosineKinaseBackgroundScores",
    Description = "Matrix of PWM scores on a set of background phosphosites published in Yaron-Barir et al. 2024",
    BiocVersion = "3.20",
    Genome = NA_character_,
    SourceType = "XLS/XLSX",
    SourceUrl = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07407-y/MediaObjects/41586_2024_7407_MOESM5_ESM.xlsx",
    SourceVersion = "16.09.2024",
    Species = "Homo sapiens",
    TaxonomyId = as.integer(9606),
    Coordinate_1_based = NA,
    DataProvider = "Nature publishing group",
    Maintainer = "Florian Geier <florian.geier@unibas.ch>",
    RDataClass = "character",
    DispatchClass = "FilePath",
    Location_Prefix = "https://zenodo.org/",
    RDataPath = "records/13831966/files/TyrosineKinaseBackgroundScores.csv.gz",
    Tags = "Proteomics:Phosohosite:Tyrosine:Kinase activity:Position weight matrix:PWM"
)

HubPub::add_resource(here::here(), metadata)

# Tyr kinase background score quantiles
metadata <- HubPub::hub_metadata(
    Title = "TyrosineKinaseBackgroundQuantiles",
    Description = "Matrix of PWM score quantiles on a set of background phosphosites published in Yaron-Barir et al. 2024",
    BiocVersion = "3.20",
    Genome = NA_character_,
    SourceType = "XLS/XLSX",
    SourceUrl = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07407-y/MediaObjects/41586_2024_7407_MOESM5_ESM.xlsx",
    SourceVersion = "16.09.2024",
    Species = "Homo sapiens",
    TaxonomyId = as.integer(9606),
    Coordinate_1_based = NA,
    DataProvider = "Nature publishing group",
    Maintainer = "Florian Geier <florian.geier@unibas.ch>",
    RDataClass = "character",
    DispatchClass = "FilePath",
    Location_Prefix = "https://zenodo.org/",
    RDataPath = "records/13831966/files/TyrosineKinaseBackgroundQuantiles.csv.gz",
    Tags = "Proteomics:Phosohosite:Tyrosine:Kinase activity:Position weight matrix:PWM"
)

HubPub::add_resource(here::here(), metadata)

