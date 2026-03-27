# shinyAMR

A Shiny app for importing antibiogram data with `AMRgen::import_pheno()`, summarising it with `AMRgen::summarise_pheno()`, and exporting NCBI or EBI-ready outputs.

## Features

- Upload CSV, TSV, TXT, compressed AST input files, or Excel workbooks (`.xls`/`.xlsx`)
- Import supported source formats including NCBI, EBI, WHONET, Vitek, MicroScan, and Sensititre
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
