# GDM-PAF CI Explorer (PAF 95% CI Explorer: Greenland, Delta, Monte Carlo Edition)

The **GDM-PAF CI Explorer** is a Shiny web application for exploring the confidence intervals of Population Attributable Fraction (PAF). This app supports two key functionalities:
1. Calculating PAF and its 95% confidence intervals using Greenland, Delta, and Monte Carlo methods.
2. Computing PAF for multiple exposure categories using the Polytomous Exposure Equation.

## How to Run
To run the app, you can visit the Shiny app at the following link:
[https://sjunlee.shinyapps.io/GDMPAFCIExplorer/](https://sjunlee.shinyapps.io/GDMPAFCIExplorer/)

---

## Features

### **Single Input Calculator**
This tab allows users to:
1. Enter their **Name**, **Affiliation**, and **Email** in the respective fields.
2. Provide inputs directly for:
   - **Tp**: Total population
   - **RR**: Relative risk
   - **CIR**: Ratio of upper-to-lower 95% confidence interval of RR
   - **Pe**: Prevalence of the risk factor in the population
3. Alternatively, upload an Excel file containing the input data.
4. Compute the PAF and its confidence intervals using Greenland, Delta, and Monte Carlo methods.
5. Download the results in an Excel format.

---

### **Polytomous Exposure Calculator**
This tab allows users to:
1. Upload an Excel file containing data for multiple exposure categories.
2. Compute the PAF using the **Polytomous Exposure Equation**:

![PAF Equation](https://latex.codecogs.com/png.latex?\dpi{150}\bg_white%20PAF%20%3D%20%5Cfrac%7B%5Csum_%7Bi%3D1%7D%5Ek%20p_i%20%28RR_i%20-%201%29%7D%7B%5Csum_%7Bi%3D1%7D%5Ek%20p_i%20%28RR_i%20-%201%29%20%2B%201%7D)

**where**:
   
   ![k](https://latex.codecogs.com/png.latex?\dpi{200}\bg_white%20k%3A%20%5Ctext%7BNumber%20of%20exposure%20categories%7D)
   
   ![pi](https://latex.codecogs.com/png.latex?\dpi{200}\bg_white%20p_i%3A%20%5Ctext%7BPrevalence%20of%20the%20%7Di%5Ctext%7B-th%20exposure%20category%7D)
   
   ![RRi](https://latex.codecogs.com/png.latex?\dpi{200}\bg_white%20RR_i%3A%20%5Ctext%7BRelative%20Risk%20of%20the%20%7Di%5Ctext%7B-th%20exposure%20category%7D)

3. Perform Monte Carlo simulations to estimate:
   - Median \( PAF \)
   - 95% confidence intervals for \( PAF \)
4. Download the results in an Excel format.

---

## How to Use

### **Single Input Calculator**
1. Go to the "Single Input Calculator" tab.
2. Fill in your **Name**, **Affiliation**, and **Email**.
3. Enter the required inputs or upload an Excel file.
4. View the calculated results on the page.
5. Download the results using the **Download Results** button.

### **Polytomous Exposure Calculator**
1. Go to the "Polytomous Exposure Calculator" tab.
2. Upload an Excel file with the required format (see below).
3. The app will automatically compute the PAF using the Polytomous Exposure Equation and display the results, including Monte Carlo confidence intervals.
4. Download the results using the **Download Results** button.

---

## Input File Format

### Single Input Calculator
The format for the Excel file in the **Single Input Calculator** tab should match the following structure:

| Tp        | Pe  | RR  | Lower | Upper |
| --------- | --- | --- | ----- | ----- |
| 1000.00   | 0.01| 1.20| 1.05  | 1.37  |
| 1000.00   | 0.90| 1.20| 1.05  | 1.37  |
| 10000.00  | 0.50| 2.00| 1.75  | 2.28  |

### Polytomous Exposure Calculator
The format for the Excel file in the **Polytomous Exposure Calculator** tab should match the following structure:

| Category     | Tp        | Pe  | RR  | Lower | Upper |
| ------------ | --------- | --- | --- | ----- | ----- |
| <18.5        | 1000.00   | 0.20| 1.50| 1.20  | 1.80  |
| 25.0-26.9    | 2000.00   | 0.30| 2.00| 1.80  | 2.20  |
| 27.0-29.9    | 1500.00   | 0.25| 2.50| 2.20  | 2.80  |
| ≥30          | 500.00    | 0.10| 3.00| 2.70  | 3.30  |

---

## Output Format

### Single Input Calculator
| **Term**             | **Description**                                                                                 |
|-----------------------|-----------------------------------------------------------------------------------------------|
| **TP**               | Total population                                                                              |
| **Pe**               | Prevalence of the risk factor in the population                                                |
| **RR**               | Relative risk                                                                                 |
| **CIR**              | Ratio of upper-to-lower 95% confidence interval of relative risk                               |
| **Var.Pe**           | Variance of Pe                                                                                |
| **O**                | Odds of Pe                                                                                    |
| **logse**            | The standard error of log(RR)                                                                 |
| **Z**                | Absolute value of beta divided by logse                                                       |
| **Pval**             | P-value                                                                                       |
| **Var.beta**         | Square of logse                                                                               |
| **AF**               | Attributable fraction                                                                         |
| **Delta.Var.AF**     | Variance of AF using Delta method                                                             |
| **Delta.low**        | Lower 95% CI of AF using Delta method                                                         |
| **Delta.up**         | Upper 95% CI of AF using Delta method                                                         |
| **Green.Var.AF**     | Variance of AF using Greenland method                                                         |
| **Green.low**        | Lower 95% CI of AF using Greenland method                                                     |
| **Green.up**         | Upper 95% CI of AF using Greenland method                                                     |
| **Monte.RR**         | Median and 95% CI for the RR from the Monte Carlo method                                      |
| **Monte.Pe**         | Median and 95% CI for the Pe from the Monte Carlo method                                      |
| **Monte.AF**         | Median and 95% CI for the AF from the Monte Carlo method                                      |
| **Monte.low**        | Lower 95% CI from the Monte Carlo method                                                     |
| **Monte.up**         | Upper 95% CI from the Monte Carlo method                                                     |

### Polytomous Exposure Calculator
| **Term**             | **Description**                                                                                 |
|-----------------------|-----------------------------------------------------------------------------------------------|
| **Sum_Pi_RR_minus1**  | The summation of \( p_i \cdot (RR_i - 1) \) for all categories                                  |
| **PAF**               | Population Attributable Fraction using the Polytomous Exposure Equation                       |
| **Monte.PAF.median**  | Median \( PAF \) value from Monte Carlo simulations                                            |
| **Monte.PAF.low**     | Lower 95% confidence limit of \( PAF \) from Monte Carlo simulations                           |
| **Monte.PAF.up**      | Upper 95% confidence limit of \( PAF \) from Monte Carlo simulations                           |

---

## Citation
Please cite this application and method as follows: https://doi.org/10.3961/jpmph.24.272
