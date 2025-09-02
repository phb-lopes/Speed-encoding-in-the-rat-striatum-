import matplotlib.pyplot as plt
from matplotlib_venn import venn3

# Positive FSI Speed cells | A = Cue | B = Trial outcome | C = Spatial choice 
a=4; b=3; c=5; ab=1; ac=2; bc=2; abc=1

sizefonte = 15
font = {'family' : 'DejaVu Sans',
        'weight' : 'bold',
        'size'   : sizefonte}
plt.rc('font', **font)
linhas = 1
colunas = 2
fig = plt.figure(figsize=(15.00,7.50))
fig.subplots_adjust(left=0.05, bottom=0.09,
    right=0.95, top=0.9, wspace=0.08, hspace=0.2)

ax1 = plt.subplot2grid((linhas, colunas), (0, 0))
venn_obj = venn3(subsets=(a, b, ab, c, ac, bc, abc),
             set_labels=('Cue', 'Trial\noutcome', 'Spatial choice'))

venn_obj.get_patch_by_id('100').set_facecolor('#09E054')  # A 
venn_obj.get_patch_by_id('010').set_facecolor('#E08C09')  # B 
venn_obj.get_patch_by_id('001').set_facecolor('#7909E0')  # C  
for patch in venn_obj.patches:
    if patch:
        patch.set_alpha(0.5)

for subset_id in ('100', '010', '110', '001', '101', '011', '111'):
    label = venn_obj.get_label_by_id(subset_id)
    if label:
        label.set_fontsize(20)  # Altere aqui para o tamanho desejado

# Negative FSI Speed cells | A = Cue | B = Trial outcome | C = Spatial choice 
a=1; b=3; c=4; ab=0; ac=1; bc=2; abc=0
ax2 = plt.subplot2grid((linhas, colunas), (0, 1))
venn_obj = venn3(subsets=(a, b, ab, c, ac, bc, abc),
             set_labels=('Cue', 'Trial\noutcome', 'Spatial choice'))

venn_obj.get_patch_by_id('100').set_facecolor('#09E054')  # A 
venn_obj.get_patch_by_id('010').set_facecolor('#E08C09')  # B 
venn_obj.get_patch_by_id('001').set_facecolor('#7909E0')  # C  
for patch in venn_obj.patches:
    if patch:
        patch.set_alpha(0.5)

for subset_id in ('100', '010', '110', '001', '101', '011', '111'):
    label = venn_obj.get_label_by_id(subset_id)
    if label:
        label.set_fontsize(20)  # Altere aqui para o tamanho desejado

save = 1
if save:
    fig.savefig('./figSup2_FRacrossCueOutcomeChoice_FSI.png', dpi=300)
else:
    plt.show()
#############################################################################
#############################################################################
##### MSN #####
#############################################################################
#############################################################################

# Positive MSN Speed cells | A = Cue | B = Trial outcome | C = Spatial choice 
a=9; b=3; c=7; ab=2; ac=1; bc=0; abc=1

sizefonte = 15
font = {'family' : 'DejaVu Sans',
        'weight' : 'bold',
        'size'   : sizefonte}
plt.rc('font', **font)
linhas = 1
colunas = 2
fig = plt.figure(figsize=(15.00,7.50))
fig.subplots_adjust(left=0.05, bottom=0.09,
    right=0.95, top=0.9, wspace=0.08, hspace=0.2)

ax1 = plt.subplot2grid((linhas, colunas), (0, 0))
venn_obj = venn3(subsets=(a, b, ab, c, ac, bc, abc),
             set_labels=('Cue', 'Trial\noutcome', 'Spatial choice'))

venn_obj.get_patch_by_id('100').set_facecolor('#09E054')  # A 
venn_obj.get_patch_by_id('010').set_facecolor('#E08C09')  # B 
venn_obj.get_patch_by_id('001').set_facecolor('#7909E0')  # C  
for patch in venn_obj.patches:
    if patch:
        patch.set_alpha(0.5)

for subset_id in ('100', '010', '110', '001', '101', '011', '111'):
    label = venn_obj.get_label_by_id(subset_id)
    if label:
        label.set_fontsize(20)  # Altere aqui para o tamanho desejado


# Negative MSN Speed cells | A = Cue | B = Trial outcome | C = Spatial choice 
a=3; b=2; c=5; ab=0; ac=1; bc=2; abc=0
ax2 = plt.subplot2grid((linhas, colunas), (0, 1))
venn_obj = venn3(subsets=(a, b, ab, c, ac, bc, abc),
             set_labels=('Cue', 'Trial\noutcome', 'Spatial choice'))

venn_obj.get_patch_by_id('100').set_facecolor('#09E054')  # A 
venn_obj.get_patch_by_id('010').set_facecolor('#E08C09')  # B 
venn_obj.get_patch_by_id('001').set_facecolor('#7909E0')  # C  
for patch in venn_obj.patches:
    if patch:
        patch.set_alpha(0.5)

for subset_id in ('100', '010', '110', '001', '101', '011', '111'):
    label = venn_obj.get_label_by_id(subset_id)
    if label:
        label.set_fontsize(20)  # Altere aqui para o tamanho desejado

save = 1
if save:
    fig.savefig('./figSup2_FRacrossCueOutcomeChoice_MSN.png', dpi=300)
else:
    plt.show()