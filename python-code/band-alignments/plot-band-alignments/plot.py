import matplotlib.pyplot as plt 
from matplotlib.ticker import MultipleLocator

materials = ['Ag$_3$SCl', 'Ag$_3$SBr', 'Ag$_3$SI', 'Ag$_3$SeCl', 'Ag$_3$SeBr', 'Ag$_3$SeI']
VB = [-5.721948101851854, -5.669427214564186, -5.245401969507099, -5.805546098927874, -5.716937898809525, -5.290616226190477]
CB = [-3.5499481018518537, -3.9034272145641857, -3.8434019695070987, -3.8105460989278743, -4.085937898809525, -4.017616226190476]

fig, axs = plt.subplots(1, figsize=(5,3))
axs.set_ylabel('Energy (eV)')
axs.set_xlim(0, 6)
axs.set_ylim(-6.5, -3)

for it in range(6):
    num_rows = 200
    alpha_val = 0.1
    y_ini = VB[it]
    y_end = -6.5
    interv = (y_end - y_ini)/num_rows
    for x in range(num_rows):
        axs.plot([it, it + .98], [y_ini, y_ini], color='coral', alpha=alpha_val)

        y_ini = y_ini + interv

        alpha_val = alpha_val + (.20)/num_rows

    num_rows = 200
    alpha_val = 0.1
    y_ini = -3
    y_end = CB[it]
    interv = (y_end - y_ini)/num_rows
    for x in range(num_rows):
        axs.plot([it, it + .98], [y_ini, y_ini], color='slateblue', alpha=alpha_val)

        y_ini = y_ini + interv

        alpha_val = alpha_val + (.20)/num_rows



for x in range(6):
    axs.plot([x, x+1], [VB[x], VB[x]], color='black', linewidth=2)
    axs.plot([x, x+1], [CB[x], CB[x]], color='black', linewidth=2)

    if x != 0:
        if VB[x] > VB[x-1]:
            axs.plot([x, x], [-8, VB[x]], color='black', linewidth=2)
        else:
            axs.plot([x, x], [-8, VB[x-1]], color='black', linewidth=2)

        if CB[x-1] > CB[x]:
            axs.plot([x, x], [CB[x], 0], color='black', linewidth=2)
        else:
            axs.plot([x, x], [CB[x-1], 0], color='black', linewidth=2)

axs.set_xticks([])

axs.text(0.5, -6.65, 'Ag$_3$SCl', ha='center', va='top', color='black')
axs.text(1.5, -6.65, 'Ag$_3$SBr', ha='center', va='top', color='black')
axs.text(2.5, -6.65, 'Ag$_3$SI', ha='center', va='top', color='black')
axs.text(3.5, -6.65, 'Ag$_3$SeCl', ha='center', va='top', color='black')
axs.text(4.5, -6.65, 'Ag$_3$SeBr', ha='center', va='top', color='black')
axs.text(5.5, -6.65, 'Ag$_3$SeI', ha='center', va='top', color='black')


axs.text(0.5, -5.82, '-5.72 eV', ha='center', va='top', color='black', fontsize=7.5)
axs.text(0.5, -3.30, '-3.55 eV', ha='center', va='top', color='black', fontsize=7.5)

axs.text(1.5, -5.77, '-5.67 eV', ha='center', va='top', color='black', fontsize=7.5)
axs.text(1.5, -3.65, '-3.90 eV', ha='center', va='top', color='black', fontsize=7.5)

axs.text(2.5, -5.35, '-5.25 eV', ha='center', va='top', color='black', fontsize=7.5)
axs.text(2.5, -3.59, '-3.84 eV', ha='center', va='top', color='black', fontsize=7.5)

axs.text(3.5, -5.91, '-5.81 eV', ha='center', va='top', color='black', fontsize=7.5)
axs.text(3.5, -3.57, '-3.81 eV', ha='center', va='top', color='black', fontsize=7.5)

axs.text(4.5, -5.82, '-5.72 eV', ha='center', va='top', color='black', fontsize=7.5)
axs.text(4.5, -3.84, '-4.09 eV', ha='center', va='top', color='black', fontsize=7.5)

axs.text(5.5, -5.39, '-5.29 eV', ha='center', va='top', color='black', fontsize=7.5)
axs.text(5.5, -3.77, '-4.02 eV', ha='center', va='top', color='black', fontsize=7.5)

axs.annotate("", xy=(0.15, -5.67), xytext=(0.15, -3.60), arrowprops=dict(arrowstyle="<->,head_length=.25,head_width=.25", mutation_scale=10, color='black'))
axs.annotate("", xy=(1.15, -5.62), xytext=(1.15, -3.95), arrowprops=dict(arrowstyle="<->,head_length=.25,head_width=.25", mutation_scale=10, color='black'))
axs.annotate("", xy=(2.15, -5.20), xytext=(2.15, -3.89), arrowprops=dict(arrowstyle="<->,head_length=.25,head_width=.25", mutation_scale=10, color='black'))
axs.annotate("", xy=(3.15, -5.76), xytext=(3.15, -3.86), arrowprops=dict(arrowstyle="<->,head_length=.25,head_width=.25", mutation_scale=10, color='black'))
axs.annotate("", xy=(4.15, -5.67), xytext=(4.15, -4.14), arrowprops=dict(arrowstyle="<->,head_length=.25,head_width=.25", mutation_scale=10, color='black'))
axs.annotate("", xy=(5.15, -5.24), xytext=(5.15, -4.07), arrowprops=dict(arrowstyle="<->,head_length=.25,head_width=.25", mutation_scale=10, color='black'))

for x in range(6):
    axs.text(x + 0.5, VB[x] +(CB[x]-VB[x])/1.8, str(format((CB[x]-VB[x]), '.2f')) + ' eV', ha='center', va='top', color='black', fontsize=7.5)

axs.axhline(-4.44, color='slategrey', linestyle='--')
axs.text(5.65, -4.2, 'H$^+$/H$_2$', ha='center', va='top', color='slategrey')
axs.axhline(-5.67, color='slategrey', linestyle='--')
axs.text(5.62, -5.72, 'H$_2$O/O$_2$', ha='center', va='top', color='slategrey')

plt.tight_layout()
plt.savefig('band-alignments.pdf')