import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator,FormatStrFormatter,MaxNLocator, LogFormatter, ScalarFormatter
import numpy as np
from matplotlib import ticker

from math import atan2,degrees
import numpy as np


mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times'] + plt.rcParams['font.serif']

plt.rcParams.update({'font.size': 26})

formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True) 
formatter.set_powerlimits((-1,1)) 

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))


### Data

labels = np.array([5,10,15,20,30,40,50,70,100,120,150,200,256])

t_pock = np.array([1.4,5.3,15.2,22.2,18.8,28.3,39.4,65.5,114.7,152.7,301.8,368.4,620.4])
t_moel = np.array([7.5,116.3,181.2,266.6,438.4,663.9,918.9,1260.4,1794.3,2195.4,2300,2500,3000])
t_gco = np.array([0.2,0.4,0.7,0.7,1.4,2,2.5,4.2,5.4,5,10.4,8.4,14.5])
t_ours_ql = np.array([8.1,7.8,7.8,8.4,8.3,9,9.3,11.2,11.9,11.7,17.1,15.1,21.1])

E_pock = np.array([6022,4820,4613,4545,4496,4481,4473,4467,4465,4464,4464,4465,4466])
E_moel = np.array([5079,4557,4509,4494,4479,4478,4480,4484,4490,4496,4525,4539,4554])
E_gco = np.array([6034,4823,4618,4550,4499,4484,4476,4470,4467,4466,4465,4464,4464])
E_ours_ql = np.array([5477,4551,4472,4466,4464,4464,4463,4464,4464,4464,4463,4463,4464])

# Plot 1. Energy vs time

p1 = ax1.plot(t_ours_ql,E_ours_ql, marker='*',  linestyle="-", zorder=3, label='Ours', 
              linewidth=3, markersize=11)
p2 = ax1.plot(t_gco, E_gco, marker='o',  linestyle="--", label='GCO', 
              linewidth=3, markersize=11)
p3 = ax1.plot(t_pock, E_pock, marker='s', linestyle=":", label='Pock', 
              linewidth=3, markersize=11)
p4 = ax1.plot(t_moel, E_moel, marker='^',  linestyle="-.",  label='Moellenhoff', 
              linewidth=3, markersize=11)

# ax1.legend(handles=[p1[0],p2[0], p3[0], p4[0]],labels=['Ours','GCO', 'Pock', 'Mollenhoff'])


ax1.set_xlabel('Time (s), log', fontsize=24)
ax1.set_ylabel('Energy', fontsize=24)
# ax.loglog()

ax1.set_xscale('log')
# ax1.set_yscale('log')

ax1.set_yticks([4500, 5000, 5500, 6000])

ax1.yaxis.set_major_formatter(formatter)
ax1.xaxis.set_major_formatter(ScalarFormatter())

ax1.set_xlim([0.19, 1000])
# ax1.set_ylim([4440, 5105])

# plt.savefig('energy_time.pdf', dpi=300,  bbox_inches='tight')

#########################################

# Plot 2. Time vs Labels

p1 = ax2.plot(labels,t_ours_ql,  marker='*',  linestyle="-", zorder=3, label='Ours', 
              linewidth=3, markersize=11)
p2 = ax2.plot(labels,t_gco,  marker='o',  linestyle="--", label='GCO', 
              linewidth=3, markersize=11)
p3 = ax2.plot(labels,t_pock,  marker='s', linestyle=":", label='Pock', 
              linewidth=3, markersize=11)
p4 = ax2.plot(labels,t_moel,  marker='^',  linestyle="-.", label='Moellenhoff', 
              linewidth=3, markersize=11)

# ax2.legend(handles=[p1[0],p2[0], p3[0], p4[0]],labels=['Ours','GCO', 'Pock', 'Mollenhoff'])

ax2.set_xlabel('Labels', fontsize=24)
ax2.set_ylabel('Time (s)', fontsize=24)
ax2.yaxis.set_label_position("right")
ax2.yaxis.tick_right()

ax2.yaxis.set_major_formatter(ScalarFormatter())
ax2.xaxis.set_major_formatter(ScalarFormatter())

ax2.set_ylim([0, 300])

ax2.ticklabel_format(style='sci')

# handles, labels = ax2.get_legend_handles_labels()
# fig.legend(handles, labels, loc=(0.27, 0.643), frameon=False, fontsize=20)

# ax.set_xlim([0, 50])
# ax.set_ylim([4460, 4560])

ax1.legend(prop={'size': 20})
ax2.legend(prop={'size': 20})

plt.savefig('results/energy_time_label.png', dpi=300,  bbox_inches='tight')