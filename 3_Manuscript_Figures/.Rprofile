source("renv/activate.R")

if (interactive() && requireNamespace("renv", quietly = TRUE)) {
  renv::activate()
  
  st <- renv::status()
  if (!st$synchronized) {
    message("\n⚠️  renv: project is out-of-sync")
    message("    Run: renv::restore()  (to install missing packages)")
    message("    Run: renv::snapshot() (if you intentionally changed packages)\n")
  }
}
