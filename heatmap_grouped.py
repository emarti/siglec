from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import seaborn as sns
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram

def heatmap_grouped(mat, column_dict, row_dict, vmin=None, vmax=None, fig=None,
                    dendrogram_linkage=None,
                    fontsize=8, dendrogram_colors=None,
                    cbar_labelsize=8,
                    heatmap_kwargs={}, fig_kwargs={}, gridspec_kwargs={},
                    cbar_kwargs={},dendrogram_kwargs={},
                    DEBUG=False):
    if vmin==None:
        vmin = mat.min()
    if vmax==None:
        vmax = mat.max()

    if fig is None:
        fig = plt.figure(**fig_kwargs)

    # Use gridspec to chunk rows and columns appropriately
    # First, we need to determine the relative widths and heights based on
    # the sizes of each dict

    column_len = [len(x) for x in column_dict.values()]
    row_len = [len(x) for x in row_dict.values()]

    if dendrogram_linkage is None:
        gs = GridSpec(len(row_dict), len(column_dict),
                    width_ratios=column_len,
                    height_ratios=row_len,
                    hspace=0.02,
                    wspace=0.02,
                    **gridspec_kwargs)
    else:
        gs_total = GridSpec(2, 1,
                            height_ratios=(0.1, 0.9),
                            hspace=0.01,
                            **gridspec_kwargs
                            )
        ax_dendrogram = fig.add_subplot(gs_total[0])
        # Change coloring of groups
        if dendrogram_colors is not None:
            hierarchy.set_link_color_palette(dendrogram_colors)
        
        dendrogram(dendrogram_linkage,
                   ax=ax_dendrogram,
                   **dendrogram_kwargs)
        ax_dendrogram.set_axis_off()

        gs = GridSpecFromSubplotSpec(
            len(row_dict), len(column_dict),
            subplot_spec=gs_total[1],
            width_ratios=column_len,
            height_ratios=row_len,
            hspace=0.02,
            wspace=0.02,)
    
    for i, (row_key, row_val) in enumerate(row_dict.items()):
        for j, (col_key, col_val) in enumerate(column_dict.items()):
            mat_ = mat.loc[row_val, col_val]

            ax_ = fig.add_subplot(gs[i, j])

            # 
            draw_colorbar = (i==0) and (j==0)
            if draw_colorbar:
                cax = fig.add_axes([0.92, 0.4, 0.01, 0.2]) 
                cax.tick_params(labelsize=8, width=0)
                heatmap_kwargs.update({'cbar_ax': cax})
                cax.tick_params(labelsize=cbar_labelsize)

            if not DEBUG:
                sns.heatmap(mat_, vmin=vmin, vmax=vmax,
                            xticklabels=(i==len(row_dict)-1),
                            yticklabels=False,
                            cbar=draw_colorbar,
                            cbar_kws=cbar_kwargs,
                            ax=ax_,
                            **heatmap_kwargs,)
            
            # Labels on left and bottom only
            if i==len(row_dict)-1:
                ax_.set_xlabel(col_key)
                ax_.tick_params(labelsize=fontsize, width=0)
            else:
                ax_.set_xlabel(None)
            ax_.set_xlabel(None)

            if j==0:
                ax_.set_ylabel(row_key,
                               rotation=0,
                               horizontalalignment='right',
                               fontsize=fontsize)
            else:
                ax_.set_ylabel(None)

