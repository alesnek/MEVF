# This script is a forward model of 36Cl concentrations for four rock samples collected at the 
# Mount Edgecumbe Volcanic Field (MEVF). The script calculates 36Cl concentrations over time
# given exposure and burial constraints and creates a 4-panel plot that compares the measured
# and modeled 36Cl concentrations for each sample. In the plots, black lines with gray shading are
# the modeled 36Cl concentrations for all simulations (n=1000), the yellow dots are measured
# 36Cl concentrations, red shading indicates intervals of burial by tephra, and blue shading
# indicates intervals of partial shielding by seawater. Variables in the "GLOBAL PARAMETERS" section
# can be modified to test different combinations of the marine limit, emplacement (eruption) age,
# tephra deposit properties.
#
# To run this script, you must have the MEVF_data.xlsx file in the same folder as the script.
#
#
#
# IMPORT LIBRARIES
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import math

mpl.rcParams['font.family'] = 'Arial'

# GLOBAL PARAMETERS
ML_elev = 10.3        # marine limit elevation; m asl
ML_age = 8000         # marine limit age; yrs ago
eruption_age = 14800  # basalt emplacement age; yrs ago
burial_age   = 13160  # tephra emplacement age; yrs ago
depth_p = 1420        # tephra thickness; cm
erate_pc = 0.115      # tephra erosion rate for MEVF21-05 ONLY; cm/yr
rho_p = 1.7           # tephra density; g/cm3
dec36 = 2.3e-6        # decay constant for 36Cl; yr^-1
Lsp = 170             # spallation e-folding length; g cm^2

# SIMULATION SETUP
filename = 'MEVF_data.xlsx'
sheet = 'model_data'
timestep = 10         # simulate concentration every 10 years
n_sim = 1000          # number of simulations per sample

# Convert ages to ka
eruption_ka = eruption_age / 1000
burial_ka   = burial_age / 1000

# LOAD DATA
df = pd.read_excel(filename, sheet_name=sheet, usecols="A:E", skiprows=0)
df = df.dropna()

# pick the ID column name
if 'Sample ID' in df.columns:
    id_col = 'Sample ID'
else:
    id_col = df.columns[0]

# DEFINE FUNCTIONS
def water_depth(t):
    return ML_elev - (ML_elev / ML_age) * t

# simulate 36Cl concentration for samples below the marine limit, which were
# partially shielded by seawater during the Holocene
def sample_below_ML(row):
    Sample_ID   = str(row.iloc[0])
    Sample_elev = float(row.iloc[1])
    cronus_P36  = float(row.iloc[2])
    meas_conc   = float(row.iloc[3])
    conc_unc    = float(row.iloc[4])

    P36_high = cronus_P36 * 1.07        # high end of production rates; at/g/yr (±7% PR uncertainty)
    P36_low  = cronus_P36 * 0.85        # low end of production rates; at/g/yr (±7% PR uncertainty and intermittent shielding by 60 g cm^-2 mass depth)              
    rho_w = 1.03                        # seawater density; g cm^-3

    exposure1_time   = eruption_age - burial_age
    full_burial_time = burial_age - ML_age

    total_time = eruption_age
    num_steps = round(total_time / timestep)
    time_vector = np.arange(0, num_steps) * timestep
    time_vector = total_time - time_vector
    time_vector_ka = time_vector / 1000

    all_Cl36 = np.zeros((n_sim, num_steps))

    # time when sample first emerges above water
    sample_above_water = np.nan
    for t in range(0, ML_age + 1, timestep):
        if water_depth(t) <= Sample_elev:
            sample_above_water = ML_age - t
            break

    for sim in range(n_sim):
        P36 = np.round(np.random.uniform(P36_low, P36_high), 2)
        Cl36 = np.zeros(num_steps)
        z_hist = np.zeros(num_steps)

        for i in range(1, num_steps):
            current_time = i * timestep

            if current_time <= exposure1_time:
                z_hist[i] = 0
                prod_rate = P36

            elif current_time <= exposure1_time + full_burial_time:
                z_hist[i] = depth_p
                prod_rate = P36 * np.exp(-depth_p * rho_p / Lsp)

            else:
                time_since_ML = current_time - (exposure1_time + full_burial_time)
                sample_depth = max(0, water_depth(time_since_ML))  # m
                if Sample_elev < sample_depth:
                    shielding_depth = (sample_depth - Sample_elev) * 100  # cm
                    prod_rate = P36 * np.exp(-shielding_depth * rho_w / Lsp)
                else:
                    shielding_depth = 0
                    prod_rate = P36
                z_hist[i] = shielding_depth

            if i == 1:
                Cl36[i-1] = 0

            Cl36[i] = Cl36[i-1] * np.exp(-dec36 * timestep) + \
                      (prod_rate / dec36) * (1 - np.exp(-dec36 * timestep))

        all_Cl36[sim, :] = Cl36

    return {
        "Sample_ID": Sample_ID,
        "Sample_elev": Sample_elev,
        "time_vector_ka": time_vector_ka,
        "all_Cl36": all_Cl36,
        "meas_conc": meas_conc,
        "conc_unc": conc_unc,
        "sample_above_water": sample_above_water
    }

# simulate 36Cl concentration for MEVF21-05, which is above the marine limit and was
# exposed due to erosion of the overlying pyroclastic deposits
def sample_above_ML(row):
    Sample_ID   = str(row.iloc[0])
    Sample_elev = float(row.iloc[1])
    cronus_P36  = float(row.iloc[2])
    meas_conc   = float(row.iloc[3])
    conc_unc    = float(row.iloc[4])

    P36_high = cronus_P36 * 1.07        # high end of production rates; at/g/yr (±7% PR uncertainty)
    P36_low  = cronus_P36 * 0.93        # low end of production rates; at/g/yr (±7% PR uncertainty)

    pyro_eroded_t = depth_p / erate_pc  # yrs until pyroclastics fully removed

    total_time = eruption_age
    num_steps = round(total_time / timestep)
    time_vector = np.arange(0, num_steps) * timestep
    time_vector = total_time - time_vector
    time_vector_ka = time_vector / 1000

    all_Cl36 = np.zeros((n_sim, num_steps))

    for sim in range(n_sim):
        P36 = np.round(np.random.uniform(P36_low, P36_high), 2)
        Cl36 = np.zeros(num_steps)
        z_hist = np.zeros(num_steps)

        for i in range(1, num_steps):
            current_time = time_vector[i]

            if current_time >= burial_age:
                depth = 0
                prod_rate = P36

            elif burial_age > current_time >= burial_age - (depth_p / erate_pc):
                burial_duration = burial_age - current_time
                depth = max(depth_p - erate_pc * burial_duration, 0)
                prod_rate = P36 * np.exp(-depth * rho_p / Lsp)

            else:
                depth = 0
                prod_rate = P36

            z_hist[i] = depth

            if i == 1:
                Cl36[i-1] = 0

            Cl36[i] = Cl36[i-1] * np.exp(-dec36 * timestep) + \
                      (prod_rate / dec36) * (1 - np.exp(-dec36 * timestep))

        all_Cl36[sim, :] = Cl36

    return {
        "Sample_ID": Sample_ID,
        "Sample_elev": Sample_elev,
        "time_vector_ka": time_vector_ka,
        "all_Cl36": all_Cl36,
        "meas_conc": meas_conc,
        "conc_unc": conc_unc,
        "pyro_eroded_t": pyro_eroded_t
    }

# subplots for samples below the marine limit
def plot_below_ML(ax, out):
    avg_Cl36 = np.mean(out["all_Cl36"], axis=0)
    for spine in ax.spines.values():
        spine.set_color('dimgray')

    for sim in range(out["all_Cl36"].shape[0]):
        ax.plot(out["time_vector_ka"], out["all_Cl36"][sim, :], color='silver', linewidth=0.5)
    ax.plot(out["time_vector_ka"], avg_Cl36, 'k', linewidth=1.2, label='Simulated $^{36}$Cl')
    ax.errorbar(0, out["meas_conc"], yerr=out["conc_unc"], fmt='o', color='black',
                markersize=7, markerfacecolor='#ffed85', markeredgewidth=1.2,
                label='Measured $^{36}$Cl', clip_on=False, zorder=5)

    ax.set_xlabel('Age (ka)', fontsize=12, fontweight='bold')
    ax.set_ylabel('$^{36}$Cl concentration (atoms g$^{-1}$)', fontsize=12, fontweight='bold')
    ax.invert_xaxis()
    ax.grid(True, linewidth=0.3)
    ax.tick_params(labelsize=10)
    ax.set_xlim(right=0)
    ax.set_ylim(bottom=0)
    ax.set_title(out["Sample_ID"], fontsize=13, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    # shading
    ax.axvspan(ML_age/1000.0, burial_age/1000.0, color='#9e2a2b', alpha=0.2)
    sab = out["sample_above_water"]
    if not np.isnan(sab):
        shield_start = sab / 1000.0
        shield_width = (ML_age - sab) / 1000.0
        ax.axvspan(shield_start, shield_start + shield_width, color='#0072BD', alpha=0.2)

    # annotations
    sab_str = f"{sab/1000.0:.1f} ka" if not np.isnan(sab) else "N/A"
    annotation_text = (
        #f"{out['Sample_ID']}\n"
        f"Eruption age: {eruption_ka:.1f} ka\n"
        #f"Burial by pyroclastics: {burial_ka:.1f} ka\n"
        f"Sample elevation: {out['Sample_elev']:.1f} m asl\n"
        f"Marine limit: {ML_elev:.1f} m asl at {ML_age/1000.0:.1f} ka\n"
        f"Sample fully above water: {sab_str}"
    )
    ax.text(0.04, 0.75, annotation_text, transform=ax.transAxes, fontsize=11,
            bbox=dict(facecolor='white', edgecolor='black'))

# subplot for sample above the marine limit
def plot_above_ML(ax, out):
    avg_Cl36 = np.mean(out["all_Cl36"], axis=0)
    for spine in ax.spines.values():
        spine.set_color('dimgray')

    for sim in range(out["all_Cl36"].shape[0]):
        ax.plot(out["time_vector_ka"], out["all_Cl36"][sim, :], color='silver', linewidth=0.5)
    ax.plot(out["time_vector_ka"], avg_Cl36, 'k', linewidth=1.2, label='Simulated $^{36}$Cl')
    ax.errorbar(0, out["meas_conc"], yerr=out["conc_unc"], fmt='o', color='black',
                markersize=7, markerfacecolor='#ffed85', markeredgewidth=1.2,
                label='Measured $^{36}$Cl', clip_on=False, zorder=5)

    ax.set_xlabel('Age (ka)', fontsize=12, fontweight='bold')
    ax.set_ylabel('$^{36}$Cl concentration (atoms g$^{-1}$)', fontsize=12, fontweight='bold')
    ax.invert_xaxis()
    ax.grid(True, linewidth=0.3)
    ax.tick_params(labelsize=10)
    ax.set_xlim(right=0)
    ax.set_ylim(bottom=0)
    ax.set_title(out["Sample_ID"], fontsize=13, fontweight='bold')
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    # burial while pyroclastics present
    pyro_eroded_t = out["pyro_eroded_t"]
    ax.axvspan(burial_age/1000.0, (burial_age - pyro_eroded_t)/1000.0, color='#9e2a2b', alpha=0.2)

    # annotations
    fully_exposed_ka = (burial_age - pyro_eroded_t) / 1000.0
    annotation_text = (
        #f"{out['Sample_ID']}\n"
        f"Eruption age: {eruption_ka:.1f} ka\n"
        #f"Burial by pyroclastics: {burial_ka:.1f} ka\n"
        f"Sample elevation: {out['Sample_elev']:.1f} m asl\n"
        f"Pyroclastics erosion rate: {erate_pc*10.0:.2f} mm/yr\n"
        f"Sample fully exposed after burial: {fully_exposed_ka:.1f} ka"
    )
    ax.text(0.04, 0.75, annotation_text, transform=ax.transAxes, fontsize=11,
            bbox=dict(facecolor='white', edgecolor='black'))


# CREATE PLOTS
n = len(df)
ncols = math.ceil(math.sqrt(n))
nrows = math.ceil(n / ncols)
fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows))
axes = np.atleast_1d(axes).ravel()

# Reserve bottom-right for MEVF21-05 if grid has at least 4 axes
reserve_special = ('MEVF21-05' in df[id_col].values) and (len(axes) >= 4)
special_ax = axes[-1] if reserve_special else None

# Plot submerged samples
ax_i = 0
for _, row in df.iterrows():
    sid = str(row[id_col])
    if reserve_special and sid == 'MEVF21-05':
        continue
    if ax_i >= len(axes) - (1 if reserve_special else 0):
        break
    out = sample_below_ML(row)
    plot_below_ML(axes[ax_i], out)
    ax_i += 1

# Plot MEVF21-05 in bottom-right if reserved; else next available
if 'MEVF21-05' in df[id_col].values:
    row_special = df[df[id_col] == 'MEVF21-05'].iloc[0]
    out_sp = sample_above_ML(row_special)
    ax_target = special_ax if reserve_special else axes[ax_i]
    plot_above_ML(ax_target, out_sp)
    if not reserve_special:
        ax_i += 1

# Remove unused axes
used = ax_i + (1 if 'MEVF21-05' in df[id_col].values else 0 if reserve_special else 0)
for j in range(used, len(axes)):
    try:
        fig.delaxes(axes[j])
    except Exception:
        pass

fig.tight_layout(pad=2.5)
#fig.savefig("MEVF_15ka.pdf", bbox_inches='tight')
plt.show()