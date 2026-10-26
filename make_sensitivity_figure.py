#!/usr/bin/env python3
"""Regenerate Supplementary Figure 15 (parameter sensitivity) at the PHYSICAL
baseline. Panel (a): local OAT S_WT (uniform +/-10%, results_local_OAT.csv).
Panel (b): Sobol S1/ST (results_sobol.csv). Matches the original image26.png style
and 1003x408 aspect. Writes sensitivity_figure.png."""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE="#2E6E8E"; RED="#BE463A"; GREEN="#2AA88B"; ORANGE="#E69637"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11,
                     "axes.linewidth":0.9,"xtick.direction":"out"})

oat=pd.read_csv("results_local_OAT.csv")
oat=oat.reindex(oat.S_WT.abs().sort_values(ascending=False).index)  # largest |S| first
sob=pd.read_csv("results_sobol.csv").sort_values("ST",ascending=False)

fig,(axA,axB)=plt.subplots(1,2,figsize=(10.03,4.08),dpi=150)

# ---- panel (a): local OAT ----
y=np.arange(len(oat))
cols=[BLUE if v>0 else RED for v in oat.S_WT]
axA.barh(y,oat.S_WT,color=cols,height=0.72)
axA.axvline(0,color="0.35",lw=0.8)
axA.set_yticks(y); axA.set_yticklabels(oat.param,fontsize=9.5)
axA.invert_yaxis()
axA.set_xlabel(r"local sensitivity  $S=(\Delta Y/Y)/(\Delta\theta/\theta)$",fontsize=10.5)
axA.spines[["top","right"]].set_visible(False)
axA.tick_params(axis="x",labelsize=9.5)
axA.text(0.97,0.06,"$k_{off0}$ and $n$ dominate;\n$F_0$, $k_0$, $K_0$ ≈ 0",
         transform=axA.transAxes,ha="right",va="bottom",fontsize=9,color="0.25")
axA.text(-0.17,1.02,"a",transform=axA.transAxes,fontsize=15,fontweight="bold",va="bottom")

# ---- panel (b): Sobol ----
yb=np.arange(len(sob)); hh=0.38
axB.barh(yb-hh/2,sob.S1,height=hh,color=GREEN,label=r"$S_1$ (first-order)")
axB.barh(yb+hh/2,sob.ST,height=hh,color=ORANGE,label=r"$S_T$ (total-effect)")
axB.axvline(0,color="0.35",lw=0.8)
axB.set_yticks(yb); axB.set_yticklabels(sob.param,fontsize=9.5)
axB.invert_yaxis()
axB.set_xlabel("Sobol index",fontsize=10.5)
axB.spines[["top","right"]].set_visible(False)
axB.tick_params(axis="x",labelsize=9.5)
axB.legend(frameon=False,fontsize=9,loc="lower right")
axB.text(-0.17,1.02,"b",transform=axB.transAxes,fontsize=15,fontweight="bold",va="bottom")

fig.tight_layout(w_pad=2.0)
fig.savefig("sensitivity_figure.png",dpi=150,facecolor="white")
print("wrote sensitivity_figure.png")
print("\npanel a order:", list(oat.param))
print("panel b order:", list(sob.param))
