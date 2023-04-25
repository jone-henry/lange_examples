import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# load data & prep
data = np.genfromtxt('bmif-Example.csv', delimiter = ',')

# get working data and subset
wd = data[1:,2:]
subset = ['bm_dx','bm_29','pb_dx','pb_dx']
subset[0] = wd[:,0:16:2]
subset[1] = wd[:,1:17:2]
subset[2] = wd[:,16::2]
subset[3] = wd[:,17::2]

# grab first column (gene names)
# remove secondary names for them and make into a list for later
h = np.genfromtxt('bmif-Example.csv',
                     delimiter = ',', dtype=str)[1:,0:1]
h2 = h.tolist()
flat_list = [item for sublist in h2 for item in sublist]
prot_names = []
for i in flat_list:
    temp = i.split(';')[0]
    prot_names.append(temp)

# making a place for the next output to go
# to make naming lists easier, a, b, c, d, = bm_dx, bm_29, pb_dx, pb_dx
out_arrays = ['a','b','b','d']

# count number of non-nans row-wise (across all samples
# if >6, change to a unique number (11) for counting later
# (the max value possible is 8)
for i in range(4):
    out_arrays[i] = np.count_nonzero(~np.isnan(subset[i]), axis=1)
    out_arrays[i][out_arrays[i] >= 6] = 11
    out_arrays[i] = out_arrays[i].tolist()

# if the element in the first list is 11, then
#  take the protein name from the second list and add
#  it to a new list
# this list is the unique proteins in each set of data
# here, a 'set' is one of the combinations of BMIF/PB/Dx/D29
unique_proteins = ['u1','u2','u3','u4']
for i in range(4):
    unique_proteins[i] = []
    for j in range(528):
        if out_arrays[i][j] == 11:
            unique_proteins[i].append(prot_names[j])

# convert to python sets
unique_sets = ['s1','s2','s3','s4']
for i in range(4):
    unique_sets[i] = set(unique_proteins[i])

# identify intersections between all combinations of 2, 3, 4 sets
# to make naming lists easier, a, b, c, d, = bm_dx, bm_29, pb_dx, pb_dx
# intersections between all 4 sets
abcd = unique_sets[0].intersection(unique_sets[1],unique_sets[2],unique_sets[3])

# intersections between 3 sets
abc = unique_sets[0].intersection(unique_sets[1],unique_sets[2])
bcd = unique_sets[1].intersection(unique_sets[2],unique_sets[3])
acd = unique_sets[0].intersection(unique_sets[2],unique_sets[3])
abd = unique_sets[0].intersection(unique_sets[1],unique_sets[3])

# intersections between 2 sets
ab = unique_sets[0].intersection(unique_sets[1])
ac = unique_sets[0].intersection(unique_sets[2])
ad = unique_sets[0].intersection(unique_sets[3])
bc = unique_sets[1].intersection(unique_sets[2])
bd = unique_sets[1].intersection(unique_sets[3])
cd = unique_sets[2].intersection(unique_sets[3])

# check if there are any unique proteins in a single set
# merge groups of sets to produce their union
u_abc = unique_sets[0].union(unique_sets[1],unique_sets[2])
u_bcd = unique_sets[1].union(unique_sets[2],unique_sets[3])
u_acd = unique_sets[0].union(unique_sets[2],unique_sets[3])
u_abd = unique_sets[0].union(unique_sets[1],unique_sets[3])

# now do the actual checking for unique proteins
a = []
for i in unique_sets[0]:
    if i not in u_bcd:
        a.append(i)
b = []
for i in unique_sets[1]:
    if i not in u_acd:
        b.append(i)
c = []
for i in unique_sets[2]:
    if i not in u_abd:
        c.append(i)
d = []
for i in unique_sets[3]:
    if i not in u_abc:
        d.append(i)

# get the plotting data - the size of each intersecting set
p_data = [len(abcd),
          len(abc), len(bcd), len(acd), len(abd),
          len(ab), len(ac), len(ad), len(bc), len(bd), len(cd),
          len(a), len(b), len(c), len(d)]

# more plotting data - the size of the base set
set_sizes = [len(unique_sets[0]), len(unique_sets[1]),
             len(unique_sets[2]), len(unique_sets[3])]

# plotting colors
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

# figure parameters
fig, ax = plt.subplots(figsize=(8.5,11))
gs = GridSpec(5, 5, figure=fig)
l, r, t, b  = 0.2, 0.8, 0.9, 0.3
plt.subplots_adjust(left=l, right=r, top=t, bottom=b, wspace=0.025)

# set up first axis; format the axis
# plot the vertical bar chart
ax1 = plt.subplot2grid((5,5), (0,2), colspan=4, rowspan=2)
for spine in ['top','right']:
    ax1.spines[spine].set_visible(False)
for spine in ['bottom','left']:
    ax1.spines[spine].set_linewidth(1)
sizes = plt.bar(range(15), p_data, color =colors_1[0], width=0.8)

# add labels above bars with the size of the intersections
for p in sizes:
   height = p.get_height()
   ax1.annotate('{}'.format(height),
      xy=(p.get_x() + p.get_width() / 2, height),
      xytext=(0, 3),
      textcoords="offset points",
      ha='center', va='bottom', fontsize=6)

# plot formatting
plt.tick_params(axis='y', width=1)
plt.ylim(0, 400)
plt.yticks(np.arange(0, 401, 100),fontsize=10, fontweight='demibold')
plt.xlim(-1, 15)
plt.xticks([],[])

# Y-axis label
x = ax1.get_position().bounds[0] - 0.15
y = ax1.get_position().bounds[1] + (ax1.get_position().bounds[3]/2)
plt.figtext(x, y, 'Intersection\nSize', fontsize=13, fontweight='demibold', ha='center', va='center')

# set up the second plot + some formatting
# this isn't a plot, but the colored dots under the primary figure
# used a scatter plot and manually made dark/light dots to
# represent which groups are being represented
ax2 = plt.subplot2grid((5,5), (2,2), colspan=4, rowspan=1)
for spine in ['top','right','bottom','left']:
    ax2.spines[spine].set_visible(False)
for i in range(4):
    plt.scatter(np.arange(0, 15, 1), np.full((15,), i, dtype=int), color=colors_2[0], s=50)
plt.xlim(-1, 15)
plt.ylim(-0.5, 3.5)
plt.xticks([],[])
plt.yticks([],[])

# add in dots for what was actually analysed
plt.scatter(np.full((4,), 0), [0,1,2,3], color=colors_1[0], s=50)

plt.scatter(np.full((3,), 1), [0,1,2], color=colors_1[0], s=50)
plt.scatter(np.full((3,), 2), [1,2,3], color=colors_1[0], s=50)
plt.scatter(np.full((3,), 3), [0,2,3], color=colors_1[0], s=50)
plt.scatter(np.full((3,), 4), [0,1,3], color=colors_1[0], s=50)

plt.scatter(np.full((2,), 5),  [0,1], color=colors_1[0], s=50)
plt.scatter(np.full((2,), 6),  [0,2], color=colors_1[0], s=50)
plt.scatter(np.full((2,), 7),  [0,3], color=colors_1[0], s=50)
plt.scatter(np.full((2,), 8),  [1,2], color=colors_1[0], s=50)
plt.scatter(np.full((2,), 9),  [1,3], color=colors_1[0], s=50)
plt.scatter(np.full((2,), 10), [2,3], color=colors_1[0], s=50)

plt.scatter(11, 0, color=colors_1[0], s=50)
plt.scatter(12, 1, color=colors_1[0], s=50)
plt.scatter(13, 2, color=colors_1[0], s=50)
plt.scatter(14, 3, color=colors_1[0], s=50)

# Y-axis labels
y_labs = ['BMIF, Dx','BMIF, D29','PB, Dx','PB, D29']
x = ax2.get_position().bounds[0] - 0.05
y = ax2.get_position().bounds[1] + (0.5/4*ax2.get_position().bounds[3])
for i in range(4):
    y1 = y + (i*1/4*ax2.get_position().bounds[3])
    plt.figtext(x, y1, y_labs[i], fontsize=11, fontweight='demibold', ha='center', va='center')

# set up the 3rd plot - the left hort. bar chart
ax3 = plt.subplot2grid((5,5), (2,0), colspan=1, rowspan=1)
for spine in ['top','left']:
    ax3.spines[spine].set_visible(False)
for spine in ['bottom','right']:
    ax3.spines[spine].set_linewidth(1)

# annotate with the size of each set
sizes_2 = plt.barh(range(4), set_sizes, 0.5)
for i, v in enumerate(set_sizes):
    ax3.text(v + 100, i, str(v), fontsize=6, ha='left', va='center')

# format the plot
plt.yticks([],[])
plt.xticks(range(0,451, 100), ['0','','200','','400'],
           fontsize=10, fontweight='demibold')
plt.tick_params(axis='x', width=1)
ax3.invert_xaxis()

# x-axis label
x = ax3.get_position().bounds[0] + (ax3.get_position().bounds[2]/2)
y = ax3.get_position().bounds[1] + ax3.get_position().bounds[3] + 0.02
plt.figtext(x, y, 'Set size', fontsize=13, fontweight='demibold', ha='center', va='center')

#plt.show()
plt.savefig('upset.svg')
