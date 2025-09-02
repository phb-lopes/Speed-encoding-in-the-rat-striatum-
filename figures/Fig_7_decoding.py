import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.patches import Patch
import pickle

# Load the data
with open('decoding_file_figure.pkl', 'rb') as f:
    data = pickle.load(f)

# Extract all variables
y_win = np.array(data['y_win'])
predictions = np.array(data['predictions'])
fsi_data = np.array(data['fsi_data'])
msn_data = np.array(data['msn_data'])
neuron_ttest_pvalue = data['neuron_ttest_pvalue']
other_data = np.array(data['other_data'])
speed_data = np.array(data['speed_data'])
loco_ttest_pvalue = data['loco_ttest_pvalue']
speed_r2 = np.array(data['speed_r2'])
surrogate_data = np.array(data['surrogate_data'])
chance_comparison_pvalue = data['chance_comparison_pvalue']
pos_r2 = np.array(data['pos_r2'])
neg_r2 = np.array(data['neg_r2'])
pos_neg_p = data['pos_neg_p']
r2_data = data['r2_data']
multi_cell_chance_r2 = data['multi_cell_chance_r2']
actual_p_values = data['actual_p_values']
multi_cell_comparisons_p = data['multi_cell_comparisons_p']
cell_numbers = np.array(data['cell_numbers'])
positions_actual = np.array(data['positions_actual'])
positions_surrogate = np.array(data['positions_surrogate'])

print("Data loaded successfully!")
print(f"FSI vs MSN p-value: {neuron_ttest_pvalue}")
print(f"Speed vs Other p-value: {loco_ttest_pvalue}")
print(f"Speed vs Surrogate p-value: {chance_comparison_pvalue}")


# Create figure with the new layout
fig = plt.figure(figsize=(13, 7))

outer_gs = GridSpec(2, 3, figure=fig, height_ratios=[1, 1])
outer_gs.update(top=0.95, bottom=0.08, left=0.06, right=0.99, 
	hspace=.4, wspace=.3)

# --- Speed Decoding (Row 1, Column 1) ---
ax_speed = fig.add_subplot(outer_gs[0, 0:2])

ax_speed.plot(y_win, label='Actual', color='black', linewidth=2.5)
ax_speed.plot(predictions, label='Decoded', color='red', linewidth=2.5)
ax_speed.set_xlabel('Time (s)', fontname='Arial', fontsize=14,weight='bold')
ax_speed.set_ylabel('Locomotion speed (cm/s)', fontname='Arial', fontsize=14, weight='bold')
ax_speed.legend(frameon=False, prop={'family': 'Arial', 'size': 13, 'weight':'bold'}, loc='upper right',bbox_to_anchor=(1.0, 1.1))
ax_speed.spines['top'].set_visible(False)
ax_speed.spines['right'].set_visible(False)
ax_speed.set_ylim(0, 50)
xticks = np.linspace(0, 1000, 6)
xticklabels = [f'{x:g}' for x in np.linspace(0, 100, 6)]
ax_speed.set_xticks(xticks)
ax_speed.set_xticklabels(xticklabels, fontsize=13)
ax_speed.set_xlim(0, 1000)
ax_speed.text(-0.1, 1.05, 'A', fontsize=18, fontweight='bold', fontname='Arial', transform=ax_speed.transAxes)

for label in ax_speed.get_yticklabels():
    label.set_fontweight('bold')
    label.set_fontname('Arial')
    label.set_fontsize(13)

for label in ax_speed.get_xticklabels():
    label.set_fontweight('bold')
    label.set_fontname('Arial')
    label.set_fontsize(13)

# --- Cell Type Analysis (Row 1, Column 2) ---
ax_cell = fig.add_subplot(outer_gs[0, 2])
bp_cell = ax_cell.boxplot([fsi_data, msn_data], showfliers=False, notch=False, patch_artist=True)
bp_cell['boxes'][0].set_facecolor('none')
bp_cell['boxes'][0].set_edgecolor('#FC4366')
bp_cell['boxes'][1].set_facecolor('none')
bp_cell['boxes'][1].set_edgecolor('#AEB2FF')

for median in bp_cell['medians']:
    median.set_color('black')

# Whiskers
bp_cell['whiskers'][0].set_color('#FC4366')
bp_cell['whiskers'][1].set_color('#FC4366')
bp_cell['whiskers'][2].set_color('#AEB2FF')
bp_cell['whiskers'][3].set_color('#AEB2FF')

# Caps
bp_cell['caps'][0].set_color('#FC4366')
bp_cell['caps'][1].set_color('#FC4366')
bp_cell['caps'][2].set_color('#AEB2FF')
bp_cell['caps'][3].set_color('#AEB2FF')


ax_cell.scatter([1.2] * len(fsi_data), fsi_data, alpha=0.5, color='#FC4366')
ax_cell.scatter([1.8] * len(msn_data), msn_data, alpha=0.5, color='#AEB2FF')
ax_cell.set_ylabel('Decoding accuracy (R²)', fontname='Arial', fontsize=14, fontweight='bold')
ax_cell.set_xticks([1, 2])
ax_cell.set_xticklabels(['FSI', 'MSN'], fontname='Arial', fontsize=14, fontweight='bold')
ax_cell.spines['top'].set_visible(False)
ax_cell.spines['right'].set_visible(False)
ax_cell.set_ylim(-.01, .4)
y_max_cell = max(np.max(fsi_data), np.max(msn_data))
bar_height_cell = y_max_cell + 0.02
ax_cell.plot([1, 2], [bar_height_cell, bar_height_cell], '-k', linewidth=1)
ax_cell.plot([1, 1], [bar_height_cell - 0.01, bar_height_cell], '-k', linewidth=1)
ax_cell.plot([2, 2], [bar_height_cell - 0.01, bar_height_cell], '-k', linewidth=1)
if neuron_ttest_pvalue < 0.05:
    ax_cell.text(1.5, bar_height_cell + 0.01, '*', ha='center', va='bottom', fontsize=14)
else:
    ax_cell.text(1.5, bar_height_cell + 0.01, 'n.s.', ha='center', va='bottom', fontsize=14)

for label in ax_cell.get_yticklabels():
    label.set_fontweight('bold')
    label.set_fontname('Arial')
    label.set_fontsize(13)

for label in ax_cell.get_xticklabels():
    label.set_fontweight('bold')
    label.set_fontname('Arial')
    label.set_fontsize(13)
ax_cell.text(-0.22, 1.05, 'B', fontsize=18, fontweight='bold', fontname='Arial', transform=ax_cell.transAxes)

# --- Positive vs Negative Optimal_Lag comparison (Row 2, Column 3) ---
ax_pos_neg = fig.add_subplot(outer_gs[1, 0])
if len(pos_r2) + len(neg_r2) > 0:
    bp_pn = ax_pos_neg.boxplot([pos_r2, neg_r2], showfliers=False, patch_artist=True)
    colors_pn = ['#2ca02c', '#2ca02c']  # green for positive, red for negative
    for patch, color in zip(bp_pn['boxes'], colors_pn):
        patch.set_facecolor('none')
        patch.set_edgecolor(color)
    # scatter individual points with slight jitter
    ax_pos_neg.scatter([1.2] * len(pos_r2), pos_r2, color=colors_pn[0], alpha=0.5)
    ax_pos_neg.scatter([1.8] * len(neg_r2), neg_r2, color=colors_pn[0], alpha=0.5)
    ax_pos_neg.set_ylim([-0.01, .4])
    ax_pos_neg.set_xticks([1, 2])
    # ax_pos_neg.set_xlim([])
    ax_pos_neg.set_xticklabels(['Positive\nspeed cells', 'Negative\nspeed cells'], fontname='Arial', fontsize=14)
    ax_pos_neg.set_ylabel('Decoding accuracy (R²)', fontname='Arial', fontweight='bold', fontsize=14)
    ax_pos_neg.spines['top'].set_visible(False)
    ax_pos_neg.spines['right'].set_visible(False)

    # annotate significance
    combined = np.concatenate([pos_r2, neg_r2]) if (len(pos_r2)+len(neg_r2))>0 else np.array([0.0])
    y_max_pn = np.nanmax(combined) if combined.size>0 else 0.0
    bar_h = y_max_pn + 0.02
    ax_pos_neg.plot([1, 2], [bar_h, bar_h], '-k', linewidth=1)
    ax_pos_neg.plot([1, 1], [bar_h - 0.01, bar_h], '-k', linewidth=1)
    ax_pos_neg.plot([2, 2], [bar_h - 0.01, bar_h], '-k', linewidth=1)
    if not np.isnan(pos_neg_p):
        sig_text = '*' if pos_neg_p < 0.05 else 'n.s.'
        ax_pos_neg.text(1.5, bar_h + 0.01, sig_text, ha='center', va='bottom', fontsize=14)
else:
    ax_pos_neg.text(0.5, 0.5, 'No data for positive vs negative lag comparison', ha='center', va='center', fontsize=14)

for label in ax_pos_neg.get_yticklabels():
    label.set_fontweight('bold')
    label.set_fontname('Arial')
    label.set_fontsize(13)

for label in ax_pos_neg.get_xticklabels():
    label.set_fontweight('bold')
    label.set_fontname('Arial')
    label.set_fontsize(13)

for median in bp_pn['medians']:
    median.set_color('black')
    
# Whiskers
bp_pn['whiskers'][0].set_color('#2ca02c')
bp_pn['whiskers'][1].set_color('#2ca02c')
bp_pn['whiskers'][2].set_color('#2ca02c')
bp_pn['whiskers'][3].set_color('#2ca02c')

# # Caps
bp_pn['caps'][0].set_color('#2ca02c')
bp_pn['caps'][1].set_color('#2ca02c')
bp_pn['caps'][2].set_color('#2ca02c')
bp_pn['caps'][3].set_color('#2ca02c')
ax_pos_neg.text(-0.23, 1.05, 'C', fontsize=18, fontweight='bold', fontname='Arial', transform=ax_pos_neg.transAxes)

# # --- Locomotion Cell Analysis (Row 2, Column 1) ---
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
ax_loco.set_ylabel('Decoding accuracy (R²)', fontname='Arial', fontweight='bold', fontsize=14)
ax_loco.set_xticks([1, 2, 3])
ax_loco.set_xticklabels(['Speed\ncells', 'Non-speed\ncells', 'Surrogate'], fontname='Arial', fontweight='bold', fontsize=14)
ax_loco.spines['top'].set_visible(False)
ax_loco.spines['right'].set_visible(False)
ax_loco.set_ylim(-0.01, 0.4)
for median in bp_loco['medians']:
    median.set_color('black')
for whisker in bp_loco['whiskers']:
    whisker.set_color(whisker.get_color())
for cap in bp_loco['caps']:
    cap.set_color(cap.get_color())
if loco_ttest_pvalue < 0.05:
    y_max_loco = max(np.max(other_data), np.max(speed_data))

ax_loco.plot([1, 3], [.395, .395], '-k', linewidth=1)
ax_loco.plot([1, 1], [.385, .395], '-k', linewidth=1)
ax_loco.plot([3, 3] ,[.385, .395], '-k', linewidth=1)
ax_loco.text(2, .38, '*', ha='center', va='bottom', fontsize=18)

ax_loco.plot([1, 2], [.37, .37], '-k', linewidth=1)
ax_loco.plot([1, 1], [.36, .37], '-k', linewidth=1)
ax_loco.plot([2, 2] ,[.36, .37], '-k', linewidth=1)
ax_loco.text(1.5, .353, '*', ha='center', va='bottom', fontsize=18)


for median in bp_loco['medians']:
    median.set_color('black')
    
# Whiskers
bp_loco['whiskers'][0].set_color('#2ca02c')
bp_loco['whiskers'][1].set_color('#2ca02c')
bp_loco['whiskers'][2].set_color('gray')
bp_loco['whiskers'][3].set_color('gray')
bp_loco['whiskers'][4].set_color('#CDCDCD')
bp_loco['whiskers'][5].set_color('#CDCDCD')


# # Caps
bp_loco['caps'][0].set_color('#2ca02c')
bp_loco['caps'][1].set_color('#2ca02c')
bp_loco['caps'][2].set_color('gray')
bp_loco['caps'][3].set_color('gray')
bp_loco['caps'][4].set_color('#CDCDCD')
bp_loco['caps'][5].set_color('#CDCDCD')

     
for label in ax_loco.get_yticklabels():
    label.set_fontweight('bold')
    label.set_fontname('Arial')
    label.set_fontsize(13)

for label in ax_loco.get_xticklabels():
    label.set_fontweight('bold')
    label.set_fontname('Arial')
    label.set_fontsize(13)
ax_loco.text(-0.24, 1.05, 'D', fontsize=18, fontweight='bold', fontname='Arial', transform=ax_loco.transAxes)



ax_multi = fig.add_subplot(outer_gs[1, 2])
# Prepare actual and surrogate data for boxplots
actual_box_data = [r2_data[i] for i in cell_numbers]
surrogate_box_data = [multi_cell_chance_r2[i] for i in cell_numbers]

# Plot actual data (green)
bp_actual = ax_multi.boxplot(actual_box_data, positions=positions_actual, widths=0.4,
                            showfliers=False, patch_artist=True, 
                            boxprops=dict(facecolor='none', edgecolor='#37a259'),
                            whiskerprops=dict(color='#37a259'),
                            capprops=dict(color='#37a259'))

# Plot surrogate data (gray)
bp_surrogate = ax_multi.boxplot(surrogate_box_data, positions=positions_surrogate, widths=0.4,
                               showfliers=False, patch_artist=True,
                               boxprops=dict(facecolor='none', edgecolor='gray'),
                               whiskerprops=dict(color='gray'),
                               capprops=dict(color='gray'))

# Set medians to black for both
for median in bp_actual['medians']:
    median.set_color('black')
for median in bp_surrogate['medians']:
    median.set_color('black')

################################################################################################
ax_multi.plot([1, 2.7], [.47, .47], '-k', linewidth=1)
ax_multi.plot([1, 1], [.46, .47], '-k', linewidth=1)
ax_multi.plot([2.7, 2.7], [.46, .47], '-k', linewidth=1)
ax_multi.text(1.75, .46, '*', ha='center', va='bottom', fontsize=13)

ax_multi.plot([2.7, 4.7], [.51, .51], '-k', linewidth=1)
ax_multi.plot([2.7, 2.7], [.5, .51], '-k', linewidth=1)
ax_multi.plot([4.7, 4.7], [.5, .51], '-k', linewidth=1)
ax_multi.text(3.75, .5, '*', ha='center', va='bottom', fontsize=13)

ax_multi.plot([4.7, 6.7], [.54, .54], '-k', linewidth=1)
ax_multi.plot([4.7, 4.7], [.53, .54], '-k', linewidth=1)
ax_multi.plot([6.7, 6.7], [.53, .54], '-k', linewidth=1)
ax_multi.text(5.75, .53, '*', ha='center', va='bottom', fontsize=13)

ax_multi.plot([6.7, 8.7], [.56, .56], '-k', linewidth=1)
ax_multi.plot([6.7, 6.7], [.55, .56], '-k', linewidth=1)
ax_multi.plot([8.7, 8.7], [.55, .56], '-k', linewidth=1)
ax_multi.text(7.75, .55, '*', ha='center', va='bottom', fontsize=13)

################################################################################################
ax_multi.plot([.7, 1.2], [.24, .24], '-k', linewidth=1)
ax_multi.plot([.7, .7], [.23, .24], '-k', linewidth=1)
ax_multi.plot([1.2, 1.2], [.23, .24], '-k', linewidth=1)
ax_multi.text(.95, .23, '*', ha='center', va='bottom', fontsize=13)

ax_multi.plot([2.7, 3.2], [.44, .44], '-k', linewidth=1)
ax_multi.plot([2.7, 2.7], [.43, .44], '-k', linewidth=1)
ax_multi.plot([3.2, 3.2], [.43, .44], '-k', linewidth=1)
ax_multi.text(2.95, .43, '*', ha='center', va='bottom', fontsize=13)

ax_multi.plot([4.7, 5.2], [.48, .48], '-k', linewidth=1)
ax_multi.plot([4.7, 4.7], [.47, .48], '-k', linewidth=1)
ax_multi.plot([5.2, 5.2], [.47, .48], '-k', linewidth=1)
ax_multi.text(4.95, .47, '*', ha='center', va='bottom', fontsize=13)

ax_multi.plot([6.7, 7.2], [.49, .49], '-k', linewidth=1)
ax_multi.plot([6.7, 6.7], [.48, .49], '-k', linewidth=1)
ax_multi.plot([7.2, 7.2], [.48, .49], '-k', linewidth=1)
ax_multi.text(6.95, .48, '*', ha='center', va='bottom', fontsize=13)

ax_multi.plot([8.7, 9.2], [.5, .5], '-k', linewidth=1)
ax_multi.plot([8.7, 8.7], [.49, .5], '-k', linewidth=1)
ax_multi.plot([9.2, 9.2], [.49, .5], '-k', linewidth=1)
ax_multi.text(8.95, .49, '*', ha='center', va='bottom', fontsize=13)

# Add legend
legend_elements = [Patch(facecolor='#37a259', edgecolor='#37a259', label='Speed cells'),
                   Patch(facecolor='gray', edgecolor='gray', label='Surrogate')]
ax_multi.legend(handles=legend_elements, bbox_to_anchor=(1, 1.15), loc='upper right', frameon=False,
 fontsize=12, ncol=2, prop={'weight': 'bold'})

ax_multi.set_xlabel('Number of cells', fontname='Arial', fontweight='bold', fontsize=14)
ax_multi.set_ylabel('Decoding accuracy (R²)', fontname='Arial', fontweight='bold', fontsize=14)
ax_multi.set_xticks(np.arange(1, 11, 2))
ax_multi.set_xticklabels(cell_numbers, fontsize=14)
ax_multi.spines['top'].set_visible(False)
ax_multi.spines['right'].set_visible(False)
ax_multi.set_xlim(0, 10)
ax_multi.set_ylim(-.01, .6)

for label in ax_multi.get_yticklabels():
    label.set_fontweight('bold')
    label.set_fontname('Arial')
    label.set_fontsize(13)

for label in ax_multi.get_xticklabels():
    label.set_fontweight('bold')
    label.set_fontname('Arial')
    label.set_fontsize(13)

ax_multi.text(-0.21, 1.05, 'E', fontsize=18, fontweight='bold', fontname='Arial', transform=ax_multi.transAxes)


plt.tight_layout()
plt.savefig('./figure7_decoding.png',dpi=300, bbox_inches='tight')
# plt.savefig('./figure7_decoding.tiff',dpi=300, compression='tiff_lzw')
plt.subplots_adjust(hspace=0.4)
plt.show()