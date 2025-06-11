
# heatmap_R

This repository contains R scripts for generating heatmaps with labeled mutations. The main script, `labelled_mutation_in_heatmap.R`, offers a way to visualize mutation data as a heatmap, with the ability to label specific mutations directly on the plot.

## Features

- Generate customizable heatmaps from mutation or numeric matrix data
- Label specific mutations or data points for clarity
- Designed for biological data visualization (e.g., genomics, transcriptomics)

## Getting Started

### Prerequisites

- R (version 4.0 or newer recommended)
- The following R packages:
  - `ggplot2`
  - `pheatmap`
  - `ComplexHeatmap`
  - `tidyverse`
  - (Other dependencies as used in your script)

Install required packages in R with:

```r
install.packages(c("ggplot2", "pheatmap", "tidyverse"))
# For ComplexHeatmap, use:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
```

### Usage

1. Clone this repository or download the script:
    ```bash
    git clone https://github.com/BeyzaCanakci/heatmap_R.git
    ```

2. Prepare your input data (see example format below).

3. Run the script in R:

    ```r
    source("labelled_mutation_in_heatmap.R")
    ```

   Or run from the command line:
   ```bash
   Rscript labelled_mutation_in_heatmap.R
   ```

### Example Input

Your input data should be a matrix or data frame with samples as columns and features (such as mutations or genes) as rows. Example:

|         | Sample1 | Sample2 | Sample3 |
|---------|---------|---------|---------|
| GeneA   |   1     |   0     |   1     |
| GeneB   |   0     |   1     |   1     |
| GeneC   |   0     |   0     |   0     |

### Output

The script will generate a heatmap plot (e.g., as a PNG or PDF file), labeling mutations or features of interest.

## Customization

- You can modify which mutations or features are labeled by editing the script or providing a list of features.
- Change color schemes and clustering options as needed.


## Contact

For questions, please open an issue or contact [@BeyzaCanakci](https://github.com/BeyzaCanakci).

---
