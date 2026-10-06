# Independent recomputation

Overall: **PASS**

- F1_C1_M2_vs_M0: n=72, b/c=39/1, RD=+0.5278, exact p=7.45786e-11, check=PASS
- F1_C2_M3_vs_M1: n=72, b/c=53/0, RD=+0.7361, exact p=2.22045e-16, check=PASS
- F2B_ReadyBind_vs_RCFree: n=54, b/c=3/4, RD=-0.0185, exact p=1, check=PASS
- F2C_ReadyBind_vs_DGFree: n=54, b/c=16/3, RD=+0.2407, exact p=0.00442505, check=PASS
- F3_Full_vs_ReadinessOnly: n=36, b/c=18/0, RD=+0.5000, exact p=7.62939e-06, check=PASS
- F3_Full_vs_BinderOnly: n=36, b/c=3/0, RD=+0.0833, exact p=0.25, check=PASS
- F4_ReadyBind_vs_TypedGuard: n=36, b/c=0/1, RD=-0.0278, exact p=1, check=PASS

## Task-template sensitivity

- F1_C1_M2_vs_M0: G=12, RD=+0.5278, 95% cluster bootstrap CI=[+0.3194, +0.7222], exact sign-flip p=0.00195312, check=PASS
- F1_C2_M3_vs_M1: G=12, RD=+0.7361, 95% cluster bootstrap CI=[+0.6111, +0.8611], exact sign-flip p=0.000488281, check=PASS
- F2B_ReadyBind_vs_RCFree: G=9, RD=-0.0185, 95% cluster bootstrap CI=[-0.1481, +0.1296], exact sign-flip p=1, check=PASS
- F2C_ReadyBind_vs_DGFree: G=9, RD=+0.2407, 95% cluster bootstrap CI=[-0.0556, +0.5556], exact sign-flip p=0.25, check=PASS
- F3_Full_vs_ReadinessOnly: G=6, RD=+0.5000, 95% cluster bootstrap CI=[+0.2500, +0.7500], exact sign-flip p=0.0625, check=PASS
- F3_Full_vs_BinderOnly: G=6, RD=+0.0833, 95% cluster bootstrap CI=[+0.0000, +0.2500], exact sign-flip p=1, check=PASS
- F4_ReadyBind_vs_TypedGuard: G=12, RD=-0.0278, 95% cluster bootstrap CI=[-0.0833, +0.0000], exact sign-flip p=1, check=PASS
