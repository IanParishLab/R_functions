# source("/home/nsaw/.Rprofile")
# # Print current libPaths for reference
# message("Using .libPaths():")
# print(.libPaths())

# # Install renv if not already installed
# if (!requireNamespace("renv", quietly = TRUE)) {
#   install.packages("renv", repos = "https://cloud.r-project.org")
# }

# # Initialize renv in bare mode (without installing packages)
# renv::init(bare = TRUE)

# # Snapshot the currently installed packages from all active library paths
# message("Creating snapshot from current libPaths...")
# renv::snapshot(library = .libPaths(), prompt = FALSE)

# # Optional: persist current libPaths in .Rprofile for future sessions
# rprofile_path <- ".Rprofile"
# libpaths_code <- sprintf(
#   ".libPaths(c(%s))",
#   paste(sprintf('\"%s\"', normalizePath(.libPaths(), winslash = "/")), collapse = ", ")
# )

# writeLines(libpaths_code, rprofile_path)
# message(sprintf("Written .libPaths() to %s:\n%s", rprofile_path, libpaths_code))
