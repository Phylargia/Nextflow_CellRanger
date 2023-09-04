# Timer hook
knitr::knit_hooks$set(start_timer = function(before, options, envir) {
  if (before) {
    # Check if the chunk has a name supplied as an option
    if ("name" %in% names(options)) {
      name <- options$name
      tic(name) # Start timer
    } else {
      stop("The 'name' option is missing in the chunk options.")
    }
  } else {
    # Execute at the end of the chunk
    if ("name" %in% names(options)) {
      name <- options$name
      toc(log = TRUE) # End timer and store data in log variable' 
    } else {
      stop("The 'name' option is missing in the chunk options.")
    }
  }
})
