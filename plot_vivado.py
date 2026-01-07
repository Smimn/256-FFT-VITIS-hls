import csv
import math
import matplotlib.pyplot as plt
from matplotlib import gridspec

# =========================
# LOAD TIME DOMAIN
# =========================
n = []
in_re = []
in_im = []

with open("time_domain_sv.csv", "r", newline="") as f:
    r = csv.DictReader(f)
    for row in r:
        n.append(int(row["n"]))
        in_re.append(float(row["re"]))
        in_im.append(float(row["im"]))

# =========================
# LOAD SPECTRUM
# =========================
k = []
pwr = []  # power = |X[k]|^2

with open("spectrum_sv.csv", "r", newline="") as f:
    r = csv.DictReader(f)
    for row in r:
        k.append(int(row["k"]))
        pwr.append(float(row["power"]))

N = len(k)

# =========================
# AMPLITUDE SPECTRUM (0 .. N-1)
# =========================
amp = [math.sqrt(max(p, 0.0)) / N for p in pwr]

k_axis = list(range(N))

# =========================
# FIGURE LAYOUT
# =========================
fig = plt.figure(figsize=(12, 8))
gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1.2])

# ---- Time domain real ----
ax1 = fig.add_subplot(gs[0, 0])
ax1.plot(n, in_re, label="Re{x[n]}")
ax1.set_title("Time domain – real part")
ax1.set_xlabel("n")
ax1.set_ylabel("amplitude")
ax1.grid(True)
ax1.legend()

# ---- Time domain imag ----
ax2 = fig.add_subplot(gs[0, 1])
ax2.plot(n, in_im, label="Im{x[n]}")
ax2.set_title("Time domain – imaginary part")
ax2.set_xlabel("n")
ax2.set_ylabel("amplitude")
ax2.grid(True)
ax2.legend()

# ---- Frequency domain (0 .. N-1, LINEAR scale) ----
ax3 = fig.add_subplot(gs[1, :])
ax3.stem(k_axis, amp, basefmt=" ")
ax3.set_title("Frequency domain – full range (0 … N−1)")
ax3.set_xlabel("k (FFT bin)")
ax3.set_ylabel("amplitude |X[k]| / N")
ax3.grid(True)
ax3.legend(["|X[k]| / N"])

# Optional: zoom y-axis a bit
ax3.set_ylim(bottom=0)

plt.tight_layout()
plt.show()
