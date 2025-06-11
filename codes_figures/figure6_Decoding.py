import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy import stats
from scipy.stats import ttest_ind, wilcoxon
from matplotlib.gridspec import GridSpecFromSubplotSpec


# load data from npz file
loaded_data = np.load('data_v2.1.npz', allow_pickle=True)

# extract data
y_win = loaded_data['y_win']
predictions = loaded_data['predictions']
relevant_r2 = loaded_data['relevant_r2']
fsi_data = loaded_data['fsi_data']
msn_data = loaded_data['msn_data']
neuron_ttest = loaded_data['neuron_ttest']
other_data = loaded_data['other_data']
speed_data = loaded_data['speed_data']
loco_ttest = loaded_data['loco_ttest']
r2_data = loaded_data['r2_data']
# speed_r2 = loaded_data['speed_r2']
surrogate_data = loaded_data['surrogate_data']

# Calculate p-value for FSI vs MSN comparison
neuron_ttest = stats.ttest_ind(fsi_data, msn_data)

# Calculate p-value for Speed vs Other cells comparison 
loco_ttest = stats.ttest_ind(other_data, speed_data)

# Calculate p-value for surrogate data
# chance_comparison_stat, chance_comparison_pvalue = wilcoxon(speed_r2, surrogate_data)

# Create figure with the new layout
fig = plt.figure(figsize=(13, 7))

fig.text(0.02, 0.95, 'A', fontsize=20, fontweight='bold', fontname='Arial')
fig.text(0.02, 0.45, 'B', fontsize=20, fontweight='bold', fontname='Arial')
fig.text(0.345, 0.45, 'C', fontsize=20, fontweight='bold', fontname='Arial')
fig.text(0.673, 0.45, 'D', fontsize=20, fontweight='bold', fontname='Arial')

outer_gs = GridSpec(2, 3, figure=fig, height_ratios=[1, 1])
outer_gs.update(top=0.95, bottom=0.08, left=0.07, right=0.99, 
	hspace=.4, wspace=.3)
# Outer GridSpec has 2 rows: top (row 1) and bottom (row 2)
# outer_gs = GridSpec(2, 1, height_ratios=[1, 1], figure=fig)

###############################################################################
# First row: 2 columns (75% left, 25% right)
###############################################################################
# gs_top = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer_gs[0], width_ratios=[3, 1])

# --- Speed Decoding (Row 1, Column 1) ---
ax_speed = fig.add_subplot(outer_gs[0, 0:3])
ax_speed.plot(range(400, 500), y_win[400:500], label='Actual', alpha=0.7, color='black',linewidth=3)
ax_speed.plot(range(400, 500), predictions[400:500], label='Decoded', alpha=0.7, color='red',linewidth=3)
ax_speed.set_xlabel('Time (s)', fontname='Arial', fontweight='bold', fontsize=14)
ax_speed.set_ylabel('Speed (Z-score)', fontname='Arial', fontweight='bold', fontsize=14)
# ax_speed.legend(frameon=False, prop={'family': 'Arial', 'size': 12})
ax_speed.legend(frameon=False, prop={'family': 'Arial', 'size': 13, 'weight':'bold'}, loc='upper right',bbox_to_anchor=(1.0, 1.1))
ax_speed.spines['top'].set_visible(False)
ax_speed.spines['right'].set_visible(False)
ax_speed.set_ylim(-1.2, 2)
xticks = np.linspace(400, 500, 6)
xticklabels = [f'{x:g}' for x in np.linspace(0, 100, 6)]
ax_speed.set_xticks(xticks)
ax_speed.set_xticklabels(xticklabels, fontsize=12)
ax_speed.set_xlim(400, 500)
# ax_speed.text(- 0.1, 1.05, 'A', fontsize=20, fontweight='bold', fontname='Arial', transform=ax_speed.transAxes)

for label in ax_speed.get_yticklabels():
    label.set_fontweight('bold')
    label.set_fontname('Arial')
    label.set_fontsize(13)

for label in ax_speed.get_xticklabels():
    label.set_fontweight('bold')
    label.set_fontname('Arial')
    label.set_fontsize(13)

# --- Cell Type Analysis (Row 1, Column 2) ---
ax_cell = fig.add_subplot(outer_gs[1, 0])
bp_cell = ax_cell.boxplot([fsi_data, msn_data], showfliers=False, notch=False, patch_artist=True)
bp_cell['boxes'][0].set_facecolor('none')
bp_cell['boxes'][0].set_edgecolor('#FC4366')
bp_cell['boxes'][1].set_facecolor('none')
bp_cell['boxes'][1].set_edgecolor('#AEB2FF')
ax_cell.scatter([1.2] * len(fsi_data), fsi_data, alpha=0.5, color='#FC4366')
ax_cell.scatter([1.8] * len(msn_data), msn_data, alpha=0.5, color='#AEB2FF')
ax_cell.set_ylabel('Decoding accuracy (R²)', fontname='Arial', fontweight='bold', fontsize=14)
ax_cell.set_xticks([1, 2])
ax_cell.set_xticklabels(['FSI', 'MSN'], fontname='Arial', fontweight='bold', fontsize=14)
ax_cell.spines['top'].set_visible(False)
ax_cell.spines['right'].set_visible(False)
for median in bp_cell['medians']:
    median.set_color('black')
for whisker in bp_cell['whiskers']:
    whisker.set_color(whisker.get_color())
for cap in bp_cell['caps']:
    cap.set_color(cap.get_color())
y_max_cell = max(np.max(fsi_data), np.max(msn_data))
bar_height_cell = y_max_cell + 0.1
ax_cell.plot([1, 2], [bar_height_cell, bar_height_cell], '-k', linewidth=1)
ax_cell.plot([1, 1], [bar_height_cell - 0.01, bar_height_cell], '-k', linewidth=1)
ax_cell.plot([2, 2], [bar_height_cell - 0.01, bar_height_cell], '-k', linewidth=1)
if neuron_ttest.pvalue < 0.05:
    ax_cell.text(1.5, bar_height_cell + 0.01, '*', ha='center', va='bottom', fontsize=20)
else:
    ax_cell.text(1.5, bar_height_cell + 0.01, 'n.s.', ha='center', va='bottom', fontsize=18)

# ax_cell.text(1.05, 1.05, 'B', fontsize=20, fontweight='bold', fontname='Arial', transform=ax_speed.transAxes)
ax_cell.set_ylim(-0.05, 0.7)

for label in ax_cell.get_yticklabels():
    label.set_fontweight('bold')
    label.set_fontname('Arial')
    label.set_fontsize(13)

for label in ax_cell.get_xticklabels():
    label.set_fontweight('bold')
    label.set_fontname('Arial')
    label.set_fontsize(13)

###############################################################################
# Second row: 3 columns (25%, 25%, 50%)
###############################################################################
# gs_bottom = GridSpecFromSubplotSpec(1, 3, subplot_spec=outer_gs[1], width_ratios=[1, 1, 1])

# --- Locomotion Cell Analysis (Row 2, Column 1) ---
ax_loco = fig.add_subplot(outer_gs[1, 1])
bp_loco = ax_loco.boxplot([speed_data, other_data, surrogate_data], showfliers=False, notch=False, patch_artist=True)
bp_loco['boxes'][0].set_facecolor('none')
bp_loco['boxes'][0].set_edgecolor('#37a259')
bp_loco['boxes'][1].set_facecolor('none')
bp_loco['boxes'][1].set_edgecolor('gray')
bp_loco['boxes'][2].set_facecolor('none')
bp_loco['boxes'][2].set_edgecolor('#CDCDCD')
ax_loco.scatter([1.3] * len(speed_data), speed_data, alpha=0.5, color='#37a259')
ax_loco.scatter([2.3] * len(other_data), other_data, alpha=0.5, color='gray')
ax_loco.scatter([3.3] * len(surrogate_data), surrogate_data, alpha=0.5, color='#CDCDCD')
ax_loco.set_ylabel('Decoding accuracy (R²)', fontname='Arial', fontweight='bold', fontsize=14)
ax_loco.set_xticks([1, 2, 3])
ax_loco.set_xticklabels(['Speed\ncells', 'Non-speed\ncells', 'Surrogate'], fontname='Arial', fontweight='bold', fontsize=14)
ax_loco.spines['top'].set_visible(False)
ax_loco.spines['right'].set_visible(False)
ax_loco.set_ylim(-0.05, 0.7)
for median in bp_loco['medians']:
    median.set_color('black')
for whisker in bp_loco['whiskers']:
    whisker.set_color(whisker.get_color())
for cap in bp_loco['caps']:
    cap.set_color(cap.get_color())
# if loco_ttest.pvalue < 0.05:
    # y_max_loco = max(np.max(other_data), np.max(speed_data))
bar_height_loco = .52
ax_loco.plot([1, 3], [bar_height_loco + .06, bar_height_loco + .06], '-k', linewidth=1)
ax_loco.plot([1, 1], [bar_height_loco + .06, bar_height_loco + .05], '-k', linewidth=1)
ax_loco.plot([3, 3] ,[bar_height_loco + .06, bar_height_loco + .05], '-k', linewidth=1)

ax_loco.plot([1, 2], [bar_height_loco - .01, bar_height_loco - .01], '-k', linewidth=1)
ax_loco.plot([1, 1], [bar_height_loco - .01, bar_height_loco - .02], '-k', linewidth=1)
ax_loco.plot([2, 2] ,[bar_height_loco - .01, bar_height_loco - .02], '-k', linewidth=1)

# ax_loco.plot([2, 3], [bar_height_loco-.3, bar_height_loco-.3], '-k', linewidth=1)
ax_loco.text(1.5, bar_height_loco - 0.03, '*', ha='center', va='bottom', fontsize=20)
ax_loco.text(2, bar_height_loco + 0.04, '*', ha='center', va='bottom', fontsize=20)
     
# ax_loco.text(- 0.1, -0.35, 'C', fontsize=20, fontweight='bold', fontname='Arial', transform=ax_speed.transAxes)


for label in ax_loco.get_yticklabels():
    label.set_fontweight('bold')
    label.set_fontname('Arial')
    label.set_fontsize(13)

for label in ax_loco.get_xticklabels():
    label.set_fontweight('bold')
    label.set_fontname('Arial')
    label.set_fontsize(13)

# --- Decoding with Multiple Cells (Row 2, Column 3) ---
ax_multi = fig.add_subplot(outer_gs[1, 2])
cell_numbers = np.arange(1, 6)
positions = np.arange(1, 11, 2)
if isinstance(r2_data, np.ndarray) and r2_data.ndim == 0:
    r2_data = r2_data.item()
box_data = [r2_data[i] for i in cell_numbers]
bp_multi = ax_multi.boxplot(box_data, positions=positions, widths=0.6,
                            showfliers=False, patch_artist=True)
for box in bp_multi['boxes']:
    box.set_facecolor('none')
    box.set_edgecolor('#37a259')
for median in bp_multi['medians']:
    median.set_color('black')
ax_multi.set_xlabel('Number of cells', fontname='Arial', fontweight='bold', fontsize=14)
ax_multi.set_ylabel('Decoding accuracy (R²)', fontname='Arial', fontweight='bold', fontsize=14)
ax_multi.set(xlim=(0.5, 5.5))
ax_multi.set_xticks(positions)
ax_multi.set_xticklabels(cell_numbers, fontsize=12, fontweight='bold')
ax_multi.spines['top'].set_visible(False)
ax_multi.spines['right'].set_visible(False)
medians = [median.get_ydata()[0] for median in bp_multi['medians']]
sems = [np.std(data, ddof=1) / np.sqrt(len(data)) for data in box_data]
offset = -0.6
adjusted_positions = [pos + offset for pos in positions]
ax_multi.plot(adjusted_positions, medians, color='#37a259', alpha=0.4, linewidth=2)
ax_multi.scatter(adjusted_positions, medians, color='#37a259', zorder=5)
ax_multi.errorbar(adjusted_positions, medians, yerr=sems, fmt='none', ecolor='#37a259',
                  alpha=0.4, capsize=5)
p_values = []
for i in range(len(box_data) - 1):
    t_stat, p_val = ttest_ind(box_data[i], box_data[i + 1])
    p_values.append(p_val)
for i, p_val in enumerate(p_values):
    x1, x2 = positions[i], positions[i + 1]
    if i == 0:
        y1 = bp_multi['caps'][1].get_ydata()[0]
        y_max_multi = y1 + 0.13
    elif i == 3:
        y1 = bp_multi['caps'][2 * i + 1].get_ydata()[0]
        y_max_multi = y1 + 0.092
    else:
        y1 = bp_multi['caps'][2 * i + 1].get_ydata()[0]
        y_max_multi = y1 + 0.06
    bar_height_multi = y_max_multi + 0.06
    ax_multi.plot([x1, x2], [bar_height_multi, bar_height_multi], '-k', linewidth=1)
    ax_multi.plot([x1, x1], [bar_height_multi - 0.01, bar_height_multi], '-k', linewidth=1)
    ax_multi.plot([x2, x2], [bar_height_multi - 0.01, bar_height_multi], '-k', linewidth=1)
    if p_val < 0.05:
        ax_multi.text((x1 + x2) / 2, bar_height_multi + 0.01, '*', ha='center', va='bottom', fontsize=20)
    else:
        ax_multi.text((x1 + x2) / 2, bar_height_multi + 0.01, 'n.s.', ha='center', va='bottom', fontsize=18)
ax_multi.set_xlim(positions[0] - 1, positions[-1] + 0.5)

# ax_multi.text(0.72, -0.35, 'D', fontsize=20, fontweight='bold', fontname='Arial', transform=ax_speed.transAxes)
ax_multi.set_ylim(-0.05, 0.7)

for label in ax_multi.get_yticklabels():
    label.set_fontweight('bold')
    label.set_fontname('Arial')
    label.set_fontsize(13)

for label in ax_multi.get_xticklabels():
    label.set_fontweight('bold')
    label.set_fontname('Arial')
    label.set_fontsize(13)
# Save and show the figure
# plt.savefig('E:/OneDrive/DocData/codes/paperFigures_v1.0/figures/figure6_decoding.tiff', dpi=300)
plt.savefig('E:/OneDrive/DocData/codes/paperFigures_v1.0/figures/figure6_decoding.tiff',
            dpi=300, compression='tiff_lzw')
plt.tight_layout()
plt.subplots_adjust(hspace=0.4)
# plt.show()
