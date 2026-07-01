# Surv-tmle-sim
Code for replicating the simulation study and real-data-analysis. Includes all data and scripts for model implementation, performance evaluation, and visualization.

---
## 📁 Folder Structure

```
real_data_analysis/
├── pbc.R                         # Real-Data-Analysis Code for Data PBC
├── rhc.R                         # Real-Data-Analysis Code for Data RHC
simulation/
├── data1/                        # Scenario 1 (Baseline)
├── data2/                        # Scenario 2 (Binary Time-Varying Covariate)
├── data3/                        # Scenario 3 (Continuous Time-Varying Covariate)
├── data4/                        # Scenario 4 (Complex Covariates and Misspecification)
├── data5/                        # Scenario 5 (Time-Varying Treatment Effect)
├── continuous/                   # Continuous-time generation, discretized at varying resolutions (K = 4, 8, 16, 30)
source/                           # Data-generating functions (generateData1-5.R, generateDataCont.R)
                                  #   and the discrete-time TMLE implementation (SurvTMLE.R / SurvTMLE2.R)
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
- `data5`: Time-Varying Treatment Effect Scenario (extends Scenario 2; the conditional treatment effect declines across follow-up)
- `continuous`: Continuous-time generation discretized at varying resolutions (event/censoring times are drawn in continuous time and binned at K = 4, 8, 16, 30 intervals)

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

- **Reproducibility / seeding.** Every driver seeds each Monte Carlo iteration with `set.seed(20250616 + i)`, and the same seed is used across the `para`, `sl2`, and `sl4` drivers, so the three estimators are evaluated on identical data sets per iteration (a paired design). Re-running reproduces the `data5` and `continuous` results exactly. For Scenarios 1–4, whose figures in the manuscript were generated from the original (unseeded) Monte Carlo runs, a seeded re-run reproduces equivalent results within Monte Carlo error.
- The **`data5`** drivers save the full per-interval ATE trajectory (one row per timepoint `t = 1..4`) for each iteration, so the recovery of the time-varying effect can be visualised.
- The **`continuous`** drivers generate survival and censoring times in continuous time, discretize each data set at `K = 4, 8, 16, 30`, run the estimator at every resolution, and save one `.Rds` per iteration with the estimated ATE at each resolution (one row per `K`). The true continuous-time ATE is a fixed constant obtained from `true_ate_cont()` in `source/generateDataCont.R`. These use `n = 2000` and 250 iterations by default.
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
