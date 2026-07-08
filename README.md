This repository contains the code and supporting materials for the paper: 'Establishing an automated Cross-Platform Evaluation of Patients Prostate Cancer Risk: Leveraging Local German Hospital Data'. Includes web-scraping scripts, data processing pipelines, and the interactive R Shiny application.

# Shiny app

The code for the [Shiny app](https://91mq7l-chernovaok.shinyapps.io/TUM-RC/) that accompanies the manuscript is in app.R

# Web interface

The redesigned web front-end of the TUM-RC calculator (Fig 4) is in the [`web/`](web/) folder and is live at https://oksana-pcrc.melihy.net. It computes the TUM-RC, SPCC, and PCRC-MRI risks client-side; see [`web/README.md`](web/README.md) for details.

## Web-Scraping: Execution Guide

The code for the web-scraping is in the folder `Python`.

For BNC2-RC

```         
python bcn2_rc.py sample_inputs.csv bcn2_rc_output.csv
```

For ERSPC34-RC

```         
python erspc34_rc.py sample_inputs.csv erspc34_rc_output.csv
```

For MSP-RC

```         
python msp_rc_calculator_parallel.py sample_inputs.csv msp_rc_output.csv 3
```

For PCRC-MRI

```         
python pcrc_mri.py sample_inputs.csv pcrc_mri_output.csv
```

For SPCC

```         
python spcc.py sample_inputs.csv spcc_output.csv
```

# Software versions

We used R version 4.5.2 and the following R package versions:

```         
pmcalibration_0.2.0 patchwork_1.3.2     ggExtra_0.11.0     
ggplot2_4.0.1       dplyr_1.1.4         openxlsx_4.2.8.1 
```

Additionally, to perform web-scraping we used Python version 3.8.
