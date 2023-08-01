## GDM-PAF CI Explorer (PAF 95% CI Explorer: Greenland, Delta, Monte Carlo Edition)
The GDM-PAF CI Explorer is a Shiny web application for exploring the confidence intervals of GDM-PAF. This app allows users to either provide their inputs directly or upload an Excel file to see the analysis results. 

## How to run
To run the app, you can visit the Shiny app as below:
https://sjunlee.shinyapps.io/GDMPAFCIExplorer/

## How to use
1. Enter your information in the "Name", "Affiliation", "Email" fields and click the "Confirm" button.

2. Provide your "Tp", "RR", "CIR", "Pe" values in the "Please enter your inputs:" section, or

3. Upload your Excel file with the data in the "Or upload an Excel file with your inputs:" section.

4. Click the "Download Data" button to download the results of the analysis.

## Excel File format
You need to upload an Excel format input data with column names matching the provided table, taking into account case sensitivity.

| Tp        | Pe  | RR  | Lower | Upper |
| --------- | --- | --- | ----- | ----- |
| 1000.00   | 0.01| 1.20| 1.05  | 1.37  |
| 1000.00   | 0.90| 1.20| 1.05  | 1.37  |
| 10000.00  | 0.50| 2.00| 1.75  | 2.28  |
| 100000.00 | 0.10| 5.00| 4.39  | 5.70  |
| 1000.00   | 0.10| 10.00| 8.77 | 11.40 |
| 10000.00  | 0.01| 1.20| 0.95  | 1.52  |
| 100000.00 | 0.90| 1.20| 0.95  | 1.52  |
| 1000.00   | 0.10| 1.20| 0.85  | 1.70  |
| 100000.00 | 0.90| 2.00| 1.41  | 2.83  |
| 100000.00 | 0.10| 1.20| 0.69  | 2.08  |
| 10000.00  | 0.01| 5.00| 2.89  | 8.66  |
| 10000.00  | 0.50| 1.20| 0.38  | 3.79  |

## Output Format
| **Term** | Description |
| --- | --- |
| **TP** | Total population |
| **Pe** | Prevalence of the risk factor in the population |
| **RR** | Relative risk |
| **CIR** | Ratio of upper-to-lower 95% confidence interval of relative risk |
| **Var.Pe** | Variance of Pe |
| **O** | Odds of Pe |
| **logse** | The standard error of log(RR) |
| **Z** | Absolute value of beta divided by logse |
| **Pval** | P-value |
| **Var.beta** | Square of logse |
| **AF** | Attributable fraction |
| **Delta.Var.AF** | Variance of AF using Delta method |
| **Delta.low** | Lower 95% CI of AF using Delta method |
| **Delta.up** | Upper 95% CI of AF using Delta method |
| **Green.Var.AF** | Variance of AF using Greenland method |
| **Green.low** | Lower 95% CI of AF using Greenland method |
| **Green.up** | Upper 95% CI of AF using Greenland method |
| **Monte.RR** | Median and 95% CI for the RR from the Monte Carlo method |
| **Monte.Pe** | Median and 95% CI for the Pe from the Monte Carlo method |
| **Monte.AF** | Median and 95% CI for the AF from the Monte Carlo method |
| **Monte.low** | Lower 95% CI from the Monte Carlo method |
| **Monte.up** | Upper 95% CI from the Monte Carlo method |
| **Sensitivity results** | Derivatives of AF, Delta.Var.AF, and Green.Var.AF with respect to each of the inputs (beta, Var.beta, TP, Pe) are calculated |

## Citation
Please cite as:
