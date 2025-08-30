# freq_iden_low_order

This repository contains the MATLAB code accompanying the paper [**“Frequency Response Identification of Low-Order Systems: Finite-Sample Analysis”**](https://arxiv.org/abs/2508.17142).

The code reproduces all experiments and numerical validations presented in the paper. It is organized around a few main scripts and helper functions.

---

## Experiments

### 1. `Low_Order_Identification.m`
This experiment reproduces the finite-sample identification results:
- **Identification error comparison** between averaging and the Low-Rank Nuclear Norm Minimization (LNNM) estimator.  
- **Bode plots** comparing the true system, averaging baseline, and LNNM.  
- **Singular values of the Loewner matrix** constructed from true, noisy, and estimated responses.

---

### 2. `Numerical_Verification_of_Sample_Complexity.m`
This experiment verifies the theoretical sample-complexity bound from the paper.  
It shows how the $\mathcal{H}_\infty$ identification error depends on different parameters:
- **Number of samples per frequency ($N$)**  
- **Number of frequency points ($M$)**  
- **Noise level ($\bar\eta$)**  

The scaling laws predicted in the theory are validated numerically.

---

## Operator Implementations

- **`loewner.m`** – Constructs the Loewner matrix from frequency points and frequency responses.  
- **`adjoint_loewner.m`** – Implements the adjoint of the Loewner operator, as defined in the paper.  
- **`inv_adj_loewner.m`** – Implements an inverse of the Loewner adjoint operator.
- **`inverse_loewner.m`** – Computes an inverse mapping from a Loewner matrix back to a frequency-response vector.

---

## Auxiliary Checks

- **`inv_adj_loewner_play.m`** and **`TEST_inv_adj_loewner_play.m`**  
  Used to examine the smallest nonzero eigenvalues of the matrices $E_{\text{in}}$ and $F_{\text{in}}$ (constructed from the Loewner operator bases).  
  These scripts demonstrate how the algebraic connectivity of these matrices are lower bounded with $M \over 2$.

- **`sum_1_z_squared.m`**  
  Validates **Lemma 7** in the paper, which states
  
  $${1 \over M^2} \sum_{r=1}^M \sum_{s=1}^M {1 \over \lvert \bar z_r - z_s \rvert^2} \le {-2 \ln\sin(\delta_m) \over (\pi - 2\delta_m)^2}.$$

  The script confirms numerically that this inequality holds and is tight.
  
---

## License

This code is released MIT License. See the [LICENSE](LICENSE) file for details.

---

## Citation
If you use this work, please cite:

```bibtex
@article{honarpisheh2025frequency,
  title={Frequency Response Identification of Low-Order Systems: Finite-Sample Analysis},
  author={Honarpisheh, Arya and Sznaier, Mario},
  journal={arXiv preprint arXiv:2508.17142},
  year={2025}
}


