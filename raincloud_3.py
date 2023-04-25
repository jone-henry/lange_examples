import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns

# Load data & prep
data = np.genfromtxt('bmif-Example.csv',
                     delimiter = ',')
# subset to get the working data
wd = data[1:,2:]

# make some variables to hold the data
bm_dx  = ['s1a','s2a','s3a','s4a','s5a','s6a','s7a','s8a']
bm_d29 = ['s1b','s2b','s3b','s4b','s5b','s6b','s7b','s8b']
pb_dx  = ['s1c','s2c','s3c','s4c','s5c','s6c','s7c','s8c']
pb_d29 = ['s1d','s2d','s3d','s4d','s5d','s6d','s7d','s8d']

# make some lists to slice up the data
starts_1 = np.arange(0,15, 2)
starts_2 = np.arange(1,16, 2)
starts_3 = np.arange(16,37, 2)
starts_4 = np.arange(17,38, 2)

# sample-wise, remove the nans (missing proteins)
for i in range(0,8):
    s = starts_1[i]
    e = s+1
    bm_dx[i] = wd[:,s:e]
    bm_dx[i] = np.log10(bm_dx[i][~np.isnan(bm_dx[i])])

    s = starts_2[i]
    e = s+1
    bm_d29[i] = wd[:,s:e]
    bm_d29[i] = np.log10(bm_d29[i][~np.isnan(bm_d29[i])])

    s = starts_3[i]
    e = s+1
    pb_dx[i] = wd[:,s:e]
    pb_dx[i] = np.log10(pb_dx[i][~np.isnan(pb_dx[i])])

    s = starts_4[i]
    e = s+1
    pb_d29[i] = wd[:,s:e]
    pb_d29[i] = np.log10(pb_d29[i][~np.isnan(pb_d29[i])])

# gather variables up to plot in sets
plotting_data = ['data_x1','data_x2','data_x3','data_x4']
plotting_data[0] = bm_dx
plotting_data[1] = bm_d29
plotting_data[2] = pb_dx
plotting_data[3] = pb_d29

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

# loop over the plotting data
for z in range(4):

    # set up the plot
    l, r, t, b  = 0.2, 0.6, 0.8, 0.3
    fig, ax = plt.subplots(figsize=(8.5,11))
    plt.subplots_adjust(left=l, right=r, top=t, bottom=b,
                        hspace=0.05, wspace=0.05)

    # boxplot + formatting
    bp = ax.boxplot(plotting_data[z], patch_artist = True,
                vert = False, showfliers=False,
                widths=0.2, positions=range(1,9))

    for median in bp['medians']:
        median.set_color('black')
    for box in bp['boxes']:
        box.set(color='k')
    for patch in bp['boxes']:
        patch.set_facecolor(colors_1[0])

    # violin plot
    vp = ax.violinplot(plotting_data[z], points=100,
                            showmeans=False, showextrema=False,
                            showmedians=False, vert=False)
    # mask bottom of violin plot + formatting
    for idx, b in enumerate(vp['bodies']):
        m = np.mean(b.get_paths()[0].vertices[:, 0])
        b.get_paths()[0].vertices[:, 1] = np.clip(b.get_paths()[0].vertices[:, 1], idx+1, idx+2)
        b.set_color(colors_2[0])
        b.set_alpha(1)

    # scatter plot
    for idx, features in enumerate(plotting_data[z]):
        y = np.full(len(features), idx + 0.8)
        idxs = np.arange(len(y))
        out = y.astype(float)
        out.flat[idxs] += np.random.uniform(low=-0.05, high=0.05, size=len(idxs))
        y = out
        ax.scatter(features, y, s=0.3, c='0.3')

    # format the plot
    plt.xlim(-0.5, 10)
    plt.tick_params(axis='both', width=1)
    plt.xticks(np.arange(0, 11, step=1), ['0','','2','','4','','6','','8','','10'],
               fontsize=10, fontweight='demibold')
    plt.yticks(np.arange(1, 9, step=1),
               fontsize=10, fontweight='demibold')
    plt.xticks(ha='center')
    for spine in ['top','right']:
        ax.spines[spine].set_visible(False)
    for spine in ['bottom','left']:
        ax.spines[spine].set_linewidth(1)

    # X-axis label
    x = ax.get_position().bounds[0] + ax.get_position().bounds[2]/2
    y = ax.get_position().bounds[1]-0.05
    plt.figtext(x, y, 'Log10(intensity)', fontsize=13, fontweight='demibold', ha='center', va='center')

    # Y-axis label
    x = ax.get_position().bounds[0] - 0.1
    y = ax.get_position().bounds[1] + ax.get_position().bounds[3]/2
    plt.figtext(x, y, 'Sample ID', fontsize=13, fontweight='demibold', ha='center', va='center')

    # Header label
    label = ['BMIF, Dx','BMIF, D29','PB, Dx','PB, D29']
    x = ax.get_position().bounds[0] + ax.get_position().bounds[2]/2
    y = ax.get_position().bounds[1] + ax.get_position().bounds[3] + 0.02
    plt.figtext(x, y, label[z], fontsize=13, fontweight='demibold', ha='center', va='center')



    save_label = ['bmif_dx','bmif_d29','pb_dx','pb_d29']
    plt.savefig('cloud_{}.svg'.format(save_label[z]))
    #plt.show()
