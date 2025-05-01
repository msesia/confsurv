# confsurv

**confsurv** is an R package that provides a method for constructing conformal survival bands under right censoring.  
It includes a consistent interface for a variety of survival models (parametric, Cox, random forest-based),  
and implements the core methods described in the accompanying paper.

### 📄 Accompanying Paper

> **Conformal Survival Bands for Risk Screening under Right-Censoring**  
> Matteo Sesia and Vladimir Svetnik (2025)  
> [link to arXiv preprint will be added soon]

The repository containing examples and code to reproduce the results in the paper using this package is available here:  
👉 [https://github.com/msesia/conformal_survival_screening](https://github.com/msesia/conformal_survival_screening)

---

## 📦 Installation

To install the latest version from GitHub:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install confsurv from GitHub
devtools::install_github("msesia/confsurv")