# GDM-PAF CI Explorer (PAF 95% CI Explorer: Greenland, Delta, Monte Carlo methods)

The **GDM-PAF CI Explorer** is a Shiny web application designed for exploring the confidence intervals of Population Attributable Fraction (PAF). This app supports three key functionalities:

1. **Single Input Calculator**: Calculate PAF and its 95% confidence intervals using Greenland, Delta, and Monte Carlo methods.
2. **Polytomous Exposure Calculator**: Compute PAF for multiple exposure categories using the Polytomous Exposure Equation.
3. **PAF Summation Calculator**: Aggregate PAF values across multiple exposures.

## How to Run
To use the application, visit the Shiny app at the following link:  
[https://sjunlee.shinyapps.io/GDMPAFCIExplorer/](https://sjunlee.shinyapps.io/GDMPAFCIExplorer/)

---

## Features

### Single Input Calculator
- Enter your **Name**, **Affiliation**, and **Email** (not stored).
- Input required values:
  - **Tp**: Total population
  - **Pe**: Prevalence of the risk factor
  - **RR**: Relative risk
  - **Lower** and **Upper** bounds of the 95% confidence interval of RR
- Alternatively, upload an Excel file with the required data.
- Outputs include:
  - PAF and confidence intervals (Greenland, Delta, Monte Carlo methods)
  - Sensitivity analysis results
- Download results in Excel format.

### Polytomous Exposure Calculator
- Upload an Excel file with data for multiple exposure categories or enter values directly.
- Compute PAF using the **Polytomous Exposure Equation**:

  $$PAF = \frac{\sum_{i=1}^{k} p_i (RR_i - 1)}{\sum_{i=1}^{k} p_i (RR_i - 1) + 1}$$

  where:
  - $p_i$: Prevalence of the $i$-th exposure category  
  - $RR_i$: Relative risk of the $i$-th exposure category

- Perform Monte Carlo simulations to estimate:
  - Median \( PAF \)
  - 95% confidence intervals
- Download results in Excel format.

### PAF Summation Calculator
- Input multiple PAF values manually or via Excel file upload.
- Summarize PAF values using the formula:

  $$Summed \, PAF = 1 - \prod_{i=1}^{n} (1 - PAF_i)$$

  This formula is derived from the study:  
  **"A Method for Summing Attributable Risks to Estimate the Total Risk Attributable to Multiple Risk Factors"**  
  [10.1093/oxfordjournals.aje.a121617](https://doi.org/10.1093/oxfordjournals.aje.a121617)

- View results and download in Excel format.


---

## How to Use

1. **Single Input Calculator**:
   - Fill in your details and input values directly or upload an Excel file.
   - View the calculated PAF, confidence intervals, and sensitivity analysis results.
   - Download the results.

2. **Polytomous Exposure Calculator**:
   - Upload an Excel file with the required format or enter data manually.
   - Compute PAF and Monte Carlo confidence intervals.
   - Download the results.

3. **PAF Summation Calculator**:
   - Provide PAF values manually or via Excel file upload.
   - Compute the summarized PAF value.
   - Download the results.

---

## Input File Format

### Single Input Calculator
| **Tp**  | **Pe** | **RR** | **Lower** | **Upper** |
|---------|--------|--------|-----------|-----------|
| 1000.00 | 0.01   | 1.20   | 1.05      | 1.37      |
| 1000.00 | 0.90   | 1.20   | 1.05      | 1.37      |

### Polytomous Exposure Calculator
| **Category** | **Tp**  | **Pe** | **RR** | **Lower** | **Upper** |
|--------------|---------|--------|--------|-----------|-----------|
| <18.5        | 1000.00 | 0.20   | 1.50   | 1.20      | 1.80      |

### PAF Summation Calculator
| **PAF** |
|---------|
| 0.10    |
| 0.20    |

---

## Output Format

### Single Input Calculator
Outputs include:
- PAF values and 95% confidence intervals (Delta, Greenland, Monte Carlo methods)
- Sensitivity analysis results for inputs \( Tp \), \( Pe \), \( RR \), and variance of \( RR \).

### Polytomous Exposure Calculator
Outputs include:
- Summed \( p_i (RR_i - 1) \)
- PAF value
- Monte Carlo median and 95% confidence intervals.

### PAF Summation Calculator
Outputs include:
- Individual PAF values
- Summed \( PAF \).

---

## Citation
Please cite this application as: [https://doi.org/10.3961/jpmph.24.272](https://doi.org/10.3961/jpmph.24.272)
