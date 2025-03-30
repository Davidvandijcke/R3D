# R3D Development Guide

This document outlines the development workflow for the R3D package.

## Development Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/R3D.git
   cd R3D
   ```

2. Install development dependencies:
   ```r
   install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown", "pkgdown", "git2r"))
   ```

3. Load the package in development mode:
   ```r
   devtools::load_all()
   ```

## Development Workflow

### Making Changes

1. Create a new branch for your feature or bugfix:
   ```bash
   git checkout -b feature/my-new-feature
   ```

2. Make your code changes.

3. Update documentation using roxygen2:
   ```r
   devtools::document()
   ```

4. Run tests:
   ```r
   devtools::test()
   ```

5. Check the package:
   ```r
   devtools::check()
   ```

### Commit Practices

Use prefixes in your commit messages to automatically categorize changes in NEWS.md:

- `feat:` or `feature:` for new features
- `fix:` or `bug:` for bug fixes
- For other changes, no specific prefix is needed

Example:
```bash
git commit -m "feat: Add support for weighted regression discontinuity"
```

### Version Numbering

We follow [semantic versioning](https://semver.org/):
- MAJOR version for incompatible API changes
- MINOR version for backward-compatible functionality additions
- PATCH version for backward-compatible bug fixes

### Updating NEWS.md

For manual updates:
```r
source("update_news.R")
update_news_entry(version = "0.1.1", type = "feature", description = "Added support for weighted data")
```

For automatic updates based on git commits:
```r
source("update_news.R")
update_news_entry(version = "0.1.1", auto_git = TRUE)
```

### Creating a New Release

1. Update the version number in the DESCRIPTION file.
2. Update NEWS.md using the helper function.
3. Run full checks:
   ```r
   devtools::check()
   ```
4. Commit the version bump:
   ```bash
   git commit -am "Bump version to x.y.z"
   ```
5. Tag the release:
   ```bash
   git tag -a vx.y.z -m "Version x.y.z"
   ```
6. Push the changes and tag:
   ```bash
   git push origin main --tags
   ```

## Website Development

### Local Website Preview

Build and preview the pkgdown site locally:
```r
pkgdown::build_site()
```

The site will be in the `docs/` directory. Open `docs/index.html` in your browser.

### Custom Articles

Create vignettes in the `vignettes/` directory:
```r
usethis::use_vignette("my-vignette-name")
```

Edit the resulting .Rmd file to create your article.

### Updating the Website Configuration

Edit `_pkgdown.yml` to customize the site structure, navigation, and appearance.

## Continuous Integration

The repository has GitHub Actions set up for:

1. `R-CMD-check`: Validates the package on multiple platforms
2. `pkgdown`: Builds and deploys the documentation website

These workflows run automatically when you push to the main branch or create a pull request.

## Need Help?

If you have questions about the development process, please open an issue on GitHub.