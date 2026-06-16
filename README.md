# AMRgenConverter

A Shiny app for importing antibiogram data, and exporting it to formats suitable for submission to NCBI or EBI-ready.

## How to use
1. The easiest way to use the app is to access the online at [apps.amrverse.org/amrgenconverter](http://apps.amrverse.org/amrgenconverter/), which requires no installation and no knowledge of R. 

2. If you have Rstudio, you can download the app code from this repository and run it locally (see below). 

3. The app uses the import and export functions from the [AMRgen](https://amrgen.org/) R package. If you are an R user, you might prefer to use the AMRgen R package directly, as this is more flexible and can be integrated into work flows.

## Features

- Upload CSV, TSV, TXT, compressed AST input files, or Excel workbooks (`.xls`/`.xlsx`)
- Import supported source formats including NCBI, EBI, WHONET, Vitek, MicroScan, Phoenix, and Sensititre
- Optionally interpret AST values with EUCAST, CLSI, or ECOFF rules during import
- Review imported long-format phenotype data in the browser
- Generate summary tables from `summarise_pheno()`
- Preview and download `export_ncbi_ast()` output as TSV
- Preview and download `export_ebi_ast()` output as TSV
- Optionally generate EBI JSON submission files and download them as a zip archive

## Run the app

Install the required packages if needed:

```r
install.packages(c("shiny", "bslib", "DT", "remotes"))
remotes::install_github("AMRverse/AMRgen")
```

Then start the app from the project directory:

```r
shiny::runApp()
```

## Deploy to shinyapps.io

Install `rsconnect` if needed:

```r
install.packages("rsconnect")
```

In the shinyapps.io dashboard, copy your account token details and run:

```r
rsconnect::setAccountInfo(
  name = "<ACCOUNT>",
  token = "<TOKEN>",
  secret = "<SECRET>"
)
```

Then deploy this app from the project directory:

```r
source("deploy_shinyapps.R")
```

Or deploy directly:

```r
rsconnect::deployApp(appDir = ".", appName = "amrgenconverter")
```

Notes:

- This app depends on `AMRgen`, which is installed from GitHub rather than CRAN.
- If deployment fails on package restore, reinstall `AMRgen` locally with `remotes::install_github("AMRverse/AMRgen")` and redeploy.
