# Custom scoring based on ecological topology of gut microbiota associated with cancer immunotherapy outcome. 
The gut microbiota influences clinical responses of cancer patients to immune checkpoint inhibitors. However, there is no consensus definition of detrimental dysbiosis. Based on metagenomics (MG) sequencing of 245 non-small cell lung cancer (NSCLC) patient feces, we constructed species-level co-abundance networks that were clustered into species interacting groups (SIG) correlating with overall survival. Thirty-seven and 45 MG species were associated with resistance (SIG1) and response (SIG2) to ICI, respectively. If combined with the quantification of Akkermansia species, this procedure allowed a person-based calculation of a topological score (TOPOSCORE) that was validated in additional 254 NSCLC patients and in 216 genitourinary cancers. Finally, this TOPOSCORE was translated into a 21 bacterial probe set-based qPCR- scoring that was validated in a prospective cohort of NSCLC patients, as well as in colorectal and melanoma patients. This approach could represent a dynamic diagnosis tool of intestinal dysbiosis to guide personalized microbiota-centered interventions.

[LINK](https://www.cell.com/cell/fulltext/S0092-8674(24)00538-5?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867424005385%3Fshowall%3Dtrue) to the orginial article

If you find it useful, please cite:
Derosa, Lisa, Valerio Iebba, Carolina Alves Costa Silva, Gianmarco Piccinno, Guojun Wu, Leonardo Lordello, Bertrand Routy, et al. “Custom Scoring Based on Ecological Topology of Gut Microbiota Associated with Cancer Immunotherapy Outcome.” Cell (Elsevier BV, June 2024). https://doi.org/10.1016/j.cell.2024.05.029.
 

# Toposcore Calculator

A Shiny application for calculating Toposcores from microbiome data to predict immunotherapy response.

Live demo: [Toposcore Calculator](https://toposcore.shinyapps.io/toposcore-calculator/)

## USING THE WEB INTERFACE

### Quick Start
1. Visit [https://toposcore.shinyapps.io/toposcore-calculator/](https://toposcore.shinyapps.io/toposcore-calculator/)
2. Upload your microbiome data CSV file
3. Click "Calculate Toposcore"
4. View and analyze results

### Detailed Usage Guide

#### 1. Data Upload
- Click "Browse..." to select your microbiome data file
- File must be in CSV format
- Required columns:
  - `Sample_id`: Unique identifier for each sample
  - Bacterial species columns with relative abundance values
  - `Akkermansia_muciniphila`: Required for classification
- Use "Download Template File" button to get a sample format

#### 2. Calculation
- After uploading data, click "Calculate Toposcore"
- Results will appear in the table below

#### 3. Viewing Results
- Results table shows:
  - Sample IDs
  - S-scores
  - SIG classifications
  - Toposcores
- Summary statistics display:
  - Distribution of SIG1/Gray/SIG2 classifications
  - Count of SIG1+ and SIG2+ Toposcores

#### 4. Visualization
- Select samples using the dropdown menu or "Select All"
- Score Distribution Plot shows:
  - Reference KDE plot
  - Selected samples as red vertical lines
  - Sample IDs labeled on the plot
- Species Heatmap displays:
  - Relative abundance of signature species
  - Akkermansia muciniphila levels
  - Color-coded by abundance levels

#### 5. Downloading Results
- "Download Results Table": Save complete results as CSV
- "Download Plot": Save the score distribution plot as PDF

## BUILDING AND MODIFYING THE APP

### Prerequisites
- R (>= 4.0.0)
- RStudio (recommended)

### Dependencies

Required R packages
```
packages <- c(
"shiny",
"DT",
"ggplot2",
"dplyr",
"readr",
"ggnewscale",
"ggrepel"
)
```

### App Project Structure
toposcore-calculator/
├── app/
│ ├── app.R # Main Shiny application
│ ├── R/
│ │ ├── toposcore_calc.R # Calculation functions
│ │ └── kde_plot.R # Plotting functions
│ └── data/
│ ├── sig1.txt # SIG1 species list
│ ├── sig2.txt # SIG2 species list
│ ├── scores_disc.csv # Reference scores
│ ├── clin_disc.csv # Clinical data
│ └── met4_valid_10rows.csv # Template file
├── dependencies.R # Package dependencies
└── README.md

### Installation

1. Clone the repository:
```
git clone https://github.com/valerioiebba/TOPOSCORE.git
cd TOPOSCORE
```
2. Install required packages in R:
   
```source("app/dependencies.R")```

### Running Locally
1. Open RStudio
2. Open `app/app.R`
3. Click "Run App" or run: ```shiny::runApp("app")```

### Deployment
To deploy to shinyapps.io:
1. Create an account at [shinyapps.io](https://www.shinyapps.io)
2. Install rsconnect package: ```rsconnect::deployApp("app")```

### Customization
- Modify `R/toposcore_calc.R` to adjust calculation methods
- Update `R/kde_plot.R` to customize visualizations
- Edit `app.R` to modify the UI/UX
- Update species lists in `data/sig1.txt` and `data/sig2.txt`

### Contributing
1. Fork the repository
2. Create a feature branch
3. Commit your changes
4. Push to the branch
5. Create a Pull Request

## License

## Contact
Valerio Iebba: valerio.iebba@gmail.com

Abhilash Dhal: adhalbiophysics@gmail.com
