## GDM-PAF CI Explorer
The GDM-PAF CI Explorer is a Shiny web application for exploring the confidence intervals of GDM-PAF. This app allows users to either provide their inputs directly or upload an Excel file to see the analysis results. 

## How to run
To run the app, you can visit the Shiny app as below:
https://sjunlee.shinyapps.io/GDMPAFCIExplorer/

## How to use
1. Enter your information in the "Name", "Affiliation", "Email" fields and click the "Confirm" button.

2. Provide your "Tp", "RR", "CIR", "Pe" values in the "Please enter your inputs:" section, or
- "CIR" is defined as the ratio of upper-to-lower 95% confidnce inverfal of ralative risk
- (e.g. if RR (95% CI) = 10.0 (2.43-41.22), then CIR = 41.20/2.43=17)  

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

## Citation
Please cite as:
