if (!requireNamespace("rsconnect", quietly = TRUE)) {
  stop("Please install the 'rsconnect' package before deploying.")
}

app_name <- Sys.getenv("SHINYAPPS_APP_NAME", unset = "amrgenconverter")

rsconnect::deployApp(
  appDir = ".",
  appName = app_name,
  appTitle = "AMRgenConverter",
  launch.browser = TRUE
)
