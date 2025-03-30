#' Update NEWS.md with a new entry
#'
#' This function adds a new entry to the NEWS.md file based on user input
#' or git commit history. When using auto_git=TRUE, it only includes commits
#' with specific prefixes (feat/feature/new, fix/bug).
#'
#' @param version Character string with the new version number (e.g., "0.1.1")
#' @param type Character string indicating entry type: "feature", "bugfix", or "other"
#' @param description Character string with the entry description
#' @param date Date of the entry. Defaults to current date
#' @param auto_git Logical. If TRUE, automatically extracts commit messages since last version
#'
#' @export
update_news_entry <- function(version = NULL, 
                              type = c("feature", "bugfix", "other"),
                              description = NULL,
                              date = Sys.Date(),
                              auto_git = FALSE) {
  
  type <- match.arg(type)
  
  # Path to NEWS.md
  news_path <- file.path(".", "NEWS.md")
  
  if (!file.exists(news_path)) {
    stop("NEWS.md file not found in current directory.")
  }
  
  # Read current NEWS.md content
  news_content <- readLines(news_path)
  
  # Check if we need to get commits from git
  if (auto_git && requireNamespace("git2r", quietly = TRUE)) {
    message("Extracting commits from git history...")
    repo <- git2r::repository(".")
    
    # Find the last version tag
    tags <- git2r::tags(repo)
    
    if (length(tags) > 0) {
      # Sort tags by date
      tag_dates <- sapply(tags, function(tag) {
        commit <- git2r::lookup(repo, git2r::target(tag))
        as.character(commit$author$when)
      })
      last_tag <- names(tag_dates)[which.max(tag_dates)]
      
      # Get commits since last tag
      last_tag_commit <- git2r::lookup(repo, git2r::target(tags[[last_tag]]))
      recent_commits <- git2r::commits(repo, n = 50)
      
      commit_idx <- which(sapply(recent_commits, function(c) identical(c$sha, last_tag_commit$sha)))
      if (length(commit_idx) > 0) {
        recent_commits <- recent_commits[1:(commit_idx-1)]
      }
      
      # Extract commit messages
      commit_msgs <- sapply(recent_commits, function(c) c$message)
      
      # Categorize commits based on prefixes - only include commits with specific prefixes
      features <- grep("^feat|^feature|^new", commit_msgs, value = TRUE, ignore.case = TRUE)
      bugfixes <- grep("^fix|^bug", commit_msgs, value = TRUE, ignore.case = TRUE)
      
      # Remove empty messages
      features <- features[nzchar(features)]
      bugfixes <- bugfixes[nzchar(bugfixes)]
      
      # Prepare entries
      feature_entries <- if (length(features) > 0) {
        paste("*", gsub("^[^:]*:\\s*", "", features))
      }
      
      bugfix_entries <- if (length(bugfixes) > 0) {
        paste("*", gsub("^[^:]*:\\s*", "", bugfixes))
      }
      
      # Check if we have any tagged commits
      if (length(features) == 0 && length(bugfixes) == 0) {
        message("No tagged commits found (feat:, fix:, etc.). Please provide changes manually.")
        auto_git <- FALSE
      } else {
        # Ask for version if not provided
        if (is.null(version)) {
          version <- readline(prompt = "Enter new version number: ")
        }
        
        # Create new entry content - only including sections with tagged commits
        new_entry <- c(
          paste("## R3D", version, paste0("(", date, ")")),
          ""
        )
        
        # Only add sections that have entries
        if (length(features) > 0) {
          new_entry <- c(new_entry, "### New Features", "", feature_entries, "")
        }
        
        if (length(bugfixes) > 0) {
          new_entry <- c(new_entry, "### Bug Fixes", "", bugfix_entries, "")
        }
      }
    } else {
      message("No tags found in repository. Please provide version and description manually.")
      auto_git <- FALSE
    }
  }
  
  # Manual entry if not using git or git extraction failed
  if (!auto_git) {
    if (is.null(version)) {
      version <- readline(prompt = "Enter new version number: ")
    }
    
    if (is.null(description)) {
      description <- readline(prompt = paste0("Enter ", type, " description: "))
    }
    
    # Create section headers based on type
    section_headers <- list(
      feature = "### New Features",
      bugfix = "### Bug Fixes",
      other = "### Other Changes"
    )
    
    # Format the entry
    new_entry <- c(
      paste("## R3D", version, paste0("(", date, ")")),
      "",
      section_headers[[type]],
      "",
      paste("*", description),
      ""
    )
  }
  
  # Find where to insert the new entry (after the title)
  title_line <- grep("^# R3D News$", news_content)
  
  if (length(title_line) == 0) {
    title_line <- 0  # If title not found, insert at the beginning
  }
  
  # Insert the new entry
  updated_content <- c(
    news_content[1:(title_line)],
    new_entry,
    news_content[(title_line+1):length(news_content)]
  )
  
  # Write the updated content
  writeLines(updated_content, news_path)
  
  message("NEWS.md updated successfully!")
}

# Example usage:
# In your console, run:
# source("update_news.R")
# update_news_entry(version = "0.1.1", type = "feature", description = "Added support for weighted data")
# 
# Or to automatically extract from git commits:
# update_news_entry(version = "0.1.1", auto_git = TRUE)