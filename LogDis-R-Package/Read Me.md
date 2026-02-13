# LogDis R Package

## Installation Guide

This package can be installed and used in several ways depending on your setup and needs. Choose the method that works best for you.

---

## Option 1 — Install directly from GitHub (Recommended)

This is the easiest and most reliable method.

```r
install.packages("remotes")
remotes::install_github("vidazamani/Metric-Skewness-for-Object-Data/LogDis-R-Package@V3")
library(LogDis)
```

This will:

* Download the package from the `V3` branch
* Build it automatically
* Install it on your system

---

## Option 2 — Install from a local package folder

You can install the package from a folder stored on your computer.

### Step 1: Get the folder

You can either:

* Download the package from GitHub (click **Code → Download ZIP**, then extract it), or
* Clone/download the repository and locate the `LogDis` folder inside it.

### Step 2: Install from the folder

**Method A — Using devtools**

```r
install.packages("devtools")
devtools::install("PATH/TO/LogDis")
library(LogDis)
```

Replace the path with the location of the `LogDis` folder on your own machine.

**Method B — Using RStudio Build tab**

1. Open the `LogDis` folder as a Project in RStudio (open the `.Rproj` file).
2. Go to the **Build** tab.

You have two main ways to install:

* Click **Install** (or **Install and Restart**)

  * RStudio will build the package from the source folder and install it directly.

* Or first click **Build Source Package**

  * This creates a `.tar.gz` (Mac/Linux) file.

* On Windows, you can also click **Build Binary Package**

  * This creates a `.zip` file (on Mac it creates a `.tgz` file).

You can then install using the built file:

```r
install.packages("path/to/LogDis_x.x.x.tar.gz", repos = NULL, type = "source")
```

or on Windows:

```r
install.packages("path/to/LogDis_x.x.x.zip", repos = NULL)
```

This option is useful for:

* Local development
* Testing changes
* Sharing a ready-to-install package file

---

## Option 3 — Compile required C++ files manually

If you prefer to work directly with the source C++ code:

```r
library(Rcpp)

sourceCpp("/PATH/TO/u2_statistic_rcpp.cpp")
sourceCpp("/PATH/TO/spd_distances.cpp")
```

Note: Update the file paths to match their location on your computer.

This option is intended for advanced users who want direct control over the compiled components.

---

## Notes

* Option 1 is recommended for most users.
* Option 2 is useful for local installation and development.
* Option 3 is for advanced use when working directly with the C++ source files.
