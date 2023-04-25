import numpy as np
import scipy.stats as sci_ss
import statsmodels.stats.multitest as ssm
import matplotlib.pyplot as plt

# Load data & prep
data = np.genfromtxt('bmif-Example.csv',
                     delimiter = ',')
# working data (=wd)
wd = data[1:,2:]

# filter for proteins presented in all of the samples (32) for simplicity
prot_num = np.count_nonzero(~np.isnan(wd), axis=1)
filtered = wd[np.where(prot_num ==32)]

# also make one list of the protein names filtered as well
# this is for plotting later on
h = np.genfromtxt('bmif-Example.csv',
                     delimiter = ',', dtype=str)[1:,0:1]
h2 = h.tolist()
flat_list = [item for sublist in h2 for item in sublist]
prot_names = []
for i in flat_list:
    temp = i.split(';')[0]
    prot_names.append(temp)

prot_names = np.array(prot_names)
prot_names = prot_names[np.where(prot_num == 32)]

# median normalize column-wise
medians = np.median(filtered, axis=0)
med_norm = np.subtract(filtered, medians)

# subset for hypothesis testing/plotting
subset = ['bm_dx','bm_29','pb_dx','pb_dx']
subset[0] = med_norm[:,0:16:2]
subset[1] = med_norm[:,1:17:2]
subset[2] = med_norm[:,16::2]
subset[3] = med_norm[:,17::2]

# do t-tests between each group for each protein
# put p-vals into array for FDR correction after
pvals_bm_dx_v_d29 = []
pvals_pb_dx_v_d29 = []
pvals_dx_bm_v_pb  = []
pvals_d29_bm_v_pb = []
variance = False
for i in range(len(subset[0])):
    s=i
    e=i+1
    # reminder for subsets: bm_dx=0; bm_29=1; pb_dx=2; pb_dx=3
    # compare dx [0] vs d29 [1] for bmif
    p_val_bm_dx_v_d29 = sci_ss.ttest_ind(a=subset[0][s:e,:].reshape(-1,1),
                                         b=subset[1][s:e,:].reshape(-1,1),
                                         equal_var=variance)
    # compare dx [2] vs d29 [3] for pb
    p_val_pb_dx_v_d29 = sci_ss.ttest_ind(a=subset[2][s:e,:].reshape(-1,1),
                                         b=subset[3][s:e,:].reshape(-1,1),
                                         equal_var=variance)
    # compare bmif [0] vs pb [2] for dx
    p_val_dx_bm_v_pb = sci_ss.ttest_ind(a=subset[0][s:e,:].reshape(-1,1),
                                        b=subset[2][s:e,:].reshape(-1,1),
                                        equal_var=variance)
    # compare bmif [1] vs pb [3] for d29
    p_val_d29_bm_v_pb = sci_ss.ttest_ind(a=subset[1][s:e,:].reshape(-1,1),
                                         b=subset[3][s:e,:].reshape(-1,1),
                                         equal_var=variance)

    pvals_bm_dx_v_d29.append(p_val_bm_dx_v_d29[1].tolist())
    pvals_pb_dx_v_d29.append(p_val_pb_dx_v_d29[1].tolist())
    pvals_dx_bm_v_pb.append(p_val_dx_bm_v_pb[1].tolist())
    pvals_d29_bm_v_pb.append(p_val_d29_bm_v_pb[1].tolist())

# the above makes a nested list
# the below flattens it
pvals_bm_dx_v_d29 = [item for sublist in pvals_bm_dx_v_d29 for item in sublist]
pvals_pb_dx_v_d29 = [item for sublist in pvals_pb_dx_v_d29 for item in sublist]
pvals_dx_bm_v_pb  = [item for sublist in pvals_dx_bm_v_pb  for item in sublist]
pvals_d29_bm_v_pb = [item for sublist in pvals_d29_bm_v_pb for item in sublist]

# correcting for FDR using Benjamini-Hochberg Procedure
# https://www.nature.com/articles/s41598-017-05949-y
p_vals_bh = ['corrected_bm_dx_v_d29', 'corrected_pb_dx_v_d29',
             'corrected_dx_bm_v_pb',  'corrected_d29_bm_v_pb']
alpha_ = 0.05
p_vals_bh[0] = ssm.multipletests(pvals_bm_dx_v_d29, alpha=alpha_, method='fdr_bh')[1]
p_vals_bh[1] = ssm.multipletests(pvals_pb_dx_v_d29, alpha=alpha_, method='fdr_bh')[1]
p_vals_bh[2] = ssm.multipletests(pvals_dx_bm_v_pb,  alpha=alpha_, method='fdr_bh')[1]
p_vals_bh[3] = ssm.multipletests(pvals_d29_bm_v_pb, alpha=alpha_, method='fdr_bh')[1]

# prepare logged corrected pvals (=lcp) for plotting
logged_corrected_pvals = ['lcp_bm_time', 'lcp_pb_time',
                          'lcp_dx_site', 'lcp_d29_site']
for i in range(4):
    logged_corrected_pvals[i] = np.log10(p_vals_bh[i])*-1

# prepare the intensities for plotting
# take means (=m) across the samples
mean = ['m_bm_dx','m_bm_29','m_pb_dx','m_pb_dx']
for i in range(4):
    mean[i] = np.mean(subset[i], axis=1)

# take the ratio (=r) of means
ratios = ['r_bm_time', 'r_pb_time',
          'r_dx_site', 'r_d29_site']
ratios[0] = mean[0]/mean[1]
ratios[1] = mean[2]/mean[3]
ratios[2] = mean[0]/mean[2]
ratios[3] = mean[1]/mean[3]

# log_2 transform (=l)
logged_ratios = ['l_bm_time', 'l_pb_time',
                 'l_dx_site', 'l_d29_site']
for i in range(4):
    logged_ratios[i] = np.log2(np.absolute(ratios[i]))

# plotting
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
l, r, t, b  = 0.3, 0.9, 0.9, 0.4
fig, ax = plt.subplots(nrows=2, ncols=2,
                       figsize=(8.5,11))
plt.subplots_adjust(left=l, right=r, top=t, bottom=b,
                    hspace=0.3)

# reminder
# xdata: logged_ratios[i]
# ydata: logged_corrected_pvals[i]

axes = ['ax1', 'ax2', 'ax3', 'ax4']
for i in range(4):
    # create the plotting axis
    pnum=i+1
    axes[i] = plt.subplot(2,2,pnum)

    # make a list of colors to use for each data point
    # depending on the value
    mapped_colors = []
    for j in range(len(logged_ratios[i])):
            if logged_corrected_pvals[i][j] > 2 and logged_ratios[i][j] > 1:
                mapped_colors.append((0.173, 0.627, 0.173, 1.0))
            elif logged_corrected_pvals[i][j] > 2 and logged_ratios[i][j] < -1:
                mapped_colors.append((0.122, 0.467, 0.706, 1.0))
            elif logged_corrected_pvals[i][j] < 2 and logged_ratios[i][j] > 1:
                mapped_colors.append((0.596, 0.875, 0.541, 1.0))
            elif logged_corrected_pvals[i][j] < 2 and logged_ratios[i][j] < -1:
                mapped_colors.append((0.682, 0.780, 0.910, 1.0))
            else:
                mapped_colors.append((1.000, 0.733, 0.471, 1.0))

    mapped_colors = np.asanyarray(mapped_colors, dtype=object)

    # plot the data, using the above colors
    size = 20
    plt.scatter(logged_ratios[i], logged_corrected_pvals[i],
                c=mapped_colors, s=size, edgecolors='0.3', linewidths=0.5)

    # add some dashed indicator lines
    plt.plot([-1,-1],[0,6.5], alpha=0.5,
             c='0.5', linestyle='--')
    plt.plot([1,1],[0,6.5], alpha=0.5,
             c='0.5', linestyle='--')
    plt.plot([-10,10],[2,2], alpha=0.5,
             c='0.5', linestyle='--')

    # formatting
    plt.xlim(-8, 8)
    plt.ylim(0, 6.5)

    # annotations for significant (> log10(pval)=2) log2(protein ratios)
    #  >1 or <-1 using the protein names
    # data will be indexed to facilitate labeling
    mask = []
    a_x  = []
    a_y  = []

    # loop over all the data to check
    for j in range(len(logged_ratios[i])):
        if logged_corrected_pvals[i][j] > 2 and logged_ratios[i][j] > 1:
            mask.append(1)
            a_x.append(logged_ratios[i][j])
            a_y.append(logged_corrected_pvals[i][j])
        elif logged_corrected_pvals[i][j] > 2 and logged_ratios[i][j] < -1:
            mask.append(-1)
            a_x.append(logged_ratios[i][j])
            a_y.append(logged_corrected_pvals[i][j])
        else:
            mask.append(0)
            a_x.append(0)
            a_y.append(0)

    # convert the list to an np array
    mask = np.array(mask)
    a_x = np.array(a_x)
    a_y = np.array(a_y)

    # isolate the names and xy coords that are significantly
    #  higher ("_h") or significantly lower ("_l")
    prot_names_h = prot_names[np.where(mask==1)]
    prot_names_l = prot_names[np.where(mask==-1)]
    axh = a_x[np.where(mask==1)]
    axl = a_x[np.where(mask==-1)]
    ayh = a_y[np.where(mask==1)]
    ayl = a_y[np.where(mask==-1)]

    # add annotations to the plot
    for j in range(len(axh)):
        axes[i].annotate(prot_names_h[j], (axh[j], ayh[j]),
                         fontsize=7, rotation=30)
    for j in range(len(axl)):
        axes[i].annotate(prot_names_l[j], (axl[j], ayl[j]),
                         fontsize=7, rotation=-30, ha='right')

# loop over axes and add ticks where appropriate
for i in (0,1,2,3):
    plt.sca(axes[i])
    plt.xticks(np.arange(-8, 8.1, 1),
               ['','','-6','','','-3','','','0','','','3','','','6','',''],
               fontsize=10, fontweight='demibold')
    plt.yticks(np.arange(0, 6.1, 1),
               ['0','','2','','4','','6'],
               fontsize=10, fontweight='demibold')
    plt.tick_params(axis='both', width=1)

# format axes
for a in axes:
    for spine in ['top','right']:
        a.spines[spine].set_visible(False)
    for spine in ['bottom','left']:
        a.spines[spine].set_linewidth(1)

# X-axis label
x = l + (0.5*(r-l))
y = b-0.05
plt.figtext(x, y, 'log2(ratio of means)', fontsize=13,
            fontweight='demibold', ha='center', va='center')

# Y-axis label
x = l -0.065
y = b + (0.5*(t-b))
plt.figtext(x, y, '-log10(p-value)', fontsize=13, fontweight='demibold',
            rotation='vertical', ha='center', va='center')

#Add titles to each plot
titles = ['BMIF: Dx/D29', 'PB: Dx/D29',
          'Dx: BMIF/PB','D29: BMIF/PB']
for i in range(4):
    x = axes[i].get_position().bounds[0] + axes[i].get_position().bounds[2]/2
    y = axes[i].get_position().bounds[1] + axes[i].get_position().bounds[2]*0.85
    plt.figtext(x, y, titles[i], fontsize=11.5, fontweight='demibold',
                ha='center', va='center')

# legend
x0 = l
y0 = b - 0.125
d  = 0.032
h = 0.022
colors = [(0.173, 0.627, 0.173, 1.0),
          (0.122, 0.467, 0.706, 1.0),
          (0.596, 0.875, 0.541, 1.0),
          (0.682, 0.780, 0.910, 1.0),
          (1.000, 0.733, 0.471, 1.0)]
text = ['log2(ratio) > 1 & log10(pval) > 2',
        'log2(ratio) < -1 & log10(pval) > 2',
        'log2(ratio) > 1 & log10(pval) < 2',
        'log2(ratio) < -1 & log10(pval) < 2',
        '-1 < log2(ratio) > 1']

for i in range(5):
    y1 = y0 - (i*d)
    fig.patches.extend([plt.Rectangle((x0, y1), h, (h*8.5/11), fill=True,
                                      color=colors[i], transform=fig.transFigure,
                                      figure=fig)])
    x2 = x0 + 2*h
    y2 = y1 + 0.4*h
    plt.figtext(x2, y2, text[i], fontsize=11.5, fontweight='demibold', ha='left', va='center')


#plt.show()
#plt.savefig('volcano.svg')