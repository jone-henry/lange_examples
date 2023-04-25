import numpy as np
import matplotlib.pyplot as plt

# Load data & prep
data = np.genfromtxt('bmif-Example.csv',
                     delimiter = ',')
# working data
wd = data[1:,2:]

# Count number of proteins and number of missing proteins
#  per sample. "NA" is taken to mean a protein is missing

prot_num = np.count_nonzero(~np.isnan(wd), axis=0)
prot_mis = np.count_nonzero(np.isnan(wd), axis=0)

# subset the data to plot
# present proteins
bm_dx_n = prot_num[0:16:2]
bm_29_n = prot_num[1:17:2]
pb_dx_n = prot_num[16::2]
pb_29_n = prot_num[17::2]

# missing proteins
bm_dx_m = prot_mis[0:16:2]
bm_29_m = prot_mis[1:17:2]
pb_dx_m = prot_mis[16::2]
pb_29_m = prot_mis[17::2]

# plotting
# cmap
colors_1=(
        (0.1215686275, 0.4666666667, 0.7058823529, 1.0),
        (1.0000000000, 0.4980392157, 0.0549019608, 1.0),
        (0.1725490196, 0.6274509804, 0.1725490196, 1.0),
        (0.8392156863, 0.1529411765, 0.1568627451, 1.0),
        (0.5803921569, 0.4039215686, 0.7411764706, 1.0),
        (0.5490196078, 0.3372549020, 0.2941176471, 1.0),
        (0.8901960784, 0.4666666667, 0.7607843137, 1.0),
        (0.4980392157, 0.4980392157, 0.4980392157, 1.0))

colors_2=(
        (0.6823529411764706, 0.7803921568627451, 0.90980392156862740, 1.0),
        (1.0000000000000000, 0.7333333333333333, 0.47058823529411764, 1.0),
        (0.5960784313725490, 0.8745098039215686, 0.54117647058823530, 1.0),
        (1.0000000000000000, 0.5960784313725490, 0.58823529411764710, 1.0),
        (0.7725490196078432, 0.6901960784313725, 0.83529411764705890, 1.0),
        (0.7686274509803922, 0.6117647058823530, 0.58039215686274510, 1.0),
        (0.9686274509803922, 0.7137254901960784, 0.82352941176470580, 1.0),
        (0.7803921568627451, 0.7803921568627451, 0.78039215686274510, 1.0))

# plot params
l, r, t, b  = 0.3, 0.7, 0.8, 0.4
fig, ax = plt.subplots(figsize=(8.5,11))
plt.subplots_adjust(left=l, right=r, top=t, bottom=b)

# y-axis placement
gap = 1
y1 = np.arange(0,   8, gap)
y2 = np.arange(9,  17, gap)
y3 = np.arange(18, 26, gap)
y4 = np.arange(27, 35, gap)

# gathering variables for plotting
ys = [y1,y2,y3,y4]
x1 = [bm_dx_n, bm_29_n, pb_dx_n, pb_29_n]
x2 = [bm_dx_m, bm_29_m, pb_dx_m, pb_29_m]

# make a h-bar chart
for i in range(4):
    ax.barh(ys[i], x1[i], color = colors_1[0])
    ax.barh(ys[i], x2[i], color = colors_2[0], left = x1[i])

# tick/spine/label formatting
ymax = 35.5
ymin = -1
plt.ylim(ymin, ymax)
plt.yticks([],[])
plt.tick_params(axis='x', width=1)
plt.xticks(np.arange(0, 501, step=100), fontsize=10, fontweight='demibold')
plt.xticks(ha='center')
for spine in ['top','right']:
    ax.spines[spine].set_visible(False)
for spine in ['bottom','left']:
    ax.spines[spine].set_linewidth(1)

# Y-axis labels
x = l - 0.075
lims = ymax - ymin
y = ax.get_position().bounds[1] + \
    (5/lims*ax.get_position().bounds[3])
labels = ['BMIF, Dx', 'BMIF, D29', 'PB, Dx', 'PB, D29']
for i in range(4):
    y1 = y + ((9*i)/lims*ax.get_position().bounds[3])
    plt.figtext(x, y1, labels[i], fontsize=13, fontweight='demibold', ha='center', va='center')

# X-axis label
x = ax.get_position().bounds[0] + ax.get_position().bounds[2]/2
y = ax.get_position().bounds[1]-0.05
plt.figtext(x, y, 'Unique Protein Count', fontsize=13, fontweight='demibold', ha='center', va='center')

# legend labels
x0 = l
y0 = ax.get_position().bounds[1] - 0.1
dx = 0.25
labels = ['Present', 'Missing']
for i in range(2):
    x1 = x0 + (dx*i)
    plt.figtext(x1, y0, labels[i], fontsize=13, fontweight='demibold', ha='center', va='center')

# legend patches
pw = 0.06
ph = 0.021
x0 = l + 0.065
y0 = ax.get_position().bounds[1] - 0.1 - (0.5*ph)

x1 = x0
fig.patches.extend([plt.Rectangle((x1, y0), pw, ph, fill=True,
                                  color=colors_1[0], transform=fig.transFigure,
                                  figure=fig)])
x1 = x0 + dx
fig.patches.extend([plt.Rectangle((x1, y0), pw, ph, fill=True,
                                  color=colors_2[0], transform=fig.transFigure,
                                  figure=fig)])

plt.savefig('protein_counts.svg')
#plt.show()

