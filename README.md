# Surv-tmle-sim
Code for replicating the simulation study and real-data-analysis. Includes all data and scripts for model implementation, performance evaluation, and visualization.

---
## 📁 Folder Structure

```
real_data_analysis/
├── pbc.R                         # Real-Data-Analysis Code for Data PBC
├── rhc.R                         # Real-Data-Analysis Code for Data RHC
simulation/
├── data1/                        # Data Scenario 1 (Baseline Scenario)
├── data2/                        # Data Scenario 2 (Binary Time-Varying Covariate)
├── data3/                        # Data Scenario 3 (Continuous Time-Varying Covariate)
├── data4/                        # Data Scenario 4 (Complex Covariate and Misspecification)
source/                           # Data Generating Function and Survival TMLE Function Implementation
```

### 📂 Mapping of Folder Names

Each data scenario folder of name prefixed with `data` contains subfolders for methods. The (sub-)folder names correspond to the following methods used in the study. Each (sub-)folder contains code for a specific method.
| **Folder Name**        | **Method**                                                |
|-------------------------|------------------------------------------------------------|
| `para`                 | TMLE with Parametric Estimator                             |
| `sl2`                  | TMLE with Super Learner Estimator (2-learner library)      |
| `sl4`                  | TMLE with Super Learner Estimator (4-learner library)      |

---

## ⚙️ How to Run Simulations

### 1. Setup

- **R required** (version ≥ 4.1 recommended)
- Install dependencies:
  ```r
  install.packages(c(
  "mlogit",
  "polspline",   
  "SuperLearner",
  "survival",
  "dplyr",
  "tidyr",
  "data.table",
  "ggplot2",
  "rstudioapi",
  ))
  ```

---

### 2. Choose a Data Simulation Scenario

Go to the relevant folder under `simulation/`:

- `data1`: Baseline Scenario
- `data2`: Binary Time-Varying Covariate Scenario
- `data3`: Continuous Time-Varying Covariate Scenario
- `data4`: Complex Covariate and Misspecification Scenario

---

### 3. Run Simulations for Each Method

**Step 1: Navigate to the Method Folder**

- Inside each scenario folder, go to the relevant method folder.

**Step 2: Run the Script**

- Open and run the script to generate results.
- To control the number of simulation iterations, set:
  - `ite`
- Total iterations = `ite`, i.e., if `ite` is set to 10, the number of simulation iterations is 10.

**Step 3: View the Results**

- Individual results are saved to `result/DATA_METHOD_X.Rds` files in the method folder. `DATA` = the data scenario folder name; `METHOD` = the method folder name; `X` = the simulation iteration index.
- These are aggregated into the `results` object in memory.
- To save the combined output, manually uncomment the final line of the code.

---

## 📝 Notes

- The output `results` is created in memory and not saved unless modified.
- Keep the file structure unchanged unless necessary.


---

## ⚙️ How to Run Real-World Data Analysis

**Step 1: Complete the setup as for the simulation.**

**Step 2: Navigate to `real_data_analysis` folder.**

**Step 3: Open and run the file `pbc.R` and / or `rhc.R`.**

- Results are aggregated into the `ATE_DATA` object in memory, representing the Average Treatment Effect Estimation of the corresponding `DATA`. `DATA` = `pbc` or `rhc`.
- The resulting condifence interval plot of ATE over timepoint by method (for the specific data) is saved as `DATA.ATE_over_Timepoint.png`. `DATA` = `PBC` or `RHC`.
- The resulting elapsed time plot by method (for the specific data) is saved as `PBC.Elapsed_Time_by_Method.png`. `DATA` = `PBC` or `RHC`.

---

## 📄 License

This project is licensed under the GPL-3.0 license.
