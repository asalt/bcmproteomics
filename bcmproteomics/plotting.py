"""
Heavily inspired and borrowed from
https://github.com/rougier/ten-rules/blob/master/figure-1.py
"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import seaborn as sb

sb.set_style('white')

__all__ = ['autolabel', 'plot_nested_bar', 'plot_summary']


def autolabel(ax):
    """
    Attach a text label above each bar displaying its height
    """
    # adapted from http://matplotlib.org/examples/api/barchart_demo.html
    rects = ax.patches
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                '%d' % int(height),
                ha='center', va='bottom')

def plot_nested_bar(df, outer, inner, ax=None, **kwargs):
    """`outer` is column name for outer data
    and `inner` is column name for inner data"""
    if ax is None:
        ax = plt.gca()
    W,w = 0.8, 0.55
    ax.set_xlim(0, len(df))
    ax.set_ylim(0, df.max().max()*1.1)
    ax.set_xticks([x+.5 for x in range(len(df))])
    ax.set_xticklabels(df.index)
    tolabel1, tolabel2 = None, None
    for i, (ix, row) in enumerate(df.iterrows()):
        if i == len(df)-1:
            tolabel1 = outer
            tolabel2 = inner
        value = row[outer]
        x, y = i+(1-W)/2, 0
        p = patches.Rectangle(
            (x, y), W, value, fill=True, transform=ax.transData,
            lw=0, alpha=0.4, label=tolabel1, **kwargs)
        ax.add_patch(p)
        value = row[inner]
        x, y = i+(1-w)/2, 0
        p = patches.Rectangle(
            (x, y), w, value, fill=True, transform=ax.transData,
            lw=0, alpha=0.8, label=tolabel2, **kwargs)
        ax.add_patch(p)
    ax.legend()
    return ax

def plot_summary(df, color=None, description=None):
    """Plot a nice summary of E2G data.
    df :: DataFrame with data from e2g instance.
          Can be a filtered subset
    color :: color to use for faces of bars
    description :: additional description in the title

    returns :: matplotlib.figure, axes
    """
    pyg_cols = ('SRA', ['IDGroup', 'IDGroup_u2g'], 'IDSet',
                ['PSMs', 'PSMs_u2g',],
                ['PeptideCount', 'PeptideCount_S'],
                ['PeptideCount_u2g', 'PeptideCount_S_u2g',]
    )
    if color is None:
        color = sb.color_palette()[0]

    fig, axes = plt.subplots(2, 3, figsize=(16, 12), sharey='row')
    alpha = .8
    for ix, col in enumerate(pyg_cols):
        ax = axes.flatten()[ix]
        data = df[col]
        if any(col == x for x in ('SRA','IDSet')) or col == ['IDGroup', 'IDGroup_u2g']:
            if col == ['IDGroup', 'IDGroup_u2g']:
                idg = data['IDGroup'].value_counts()
                idg_u2g = data['IDGroup_u2g'].value_counts()
                plot_data = pd.concat((idg, idg_u2g), axis=1)
                plot_data = plot_data[plot_data.index != 0]
                plot_data.index = plot_data.index.astype(int)
                plot_nested_bar(plot_data, outer='IDGroup',
                                inner='IDGroup_u2g', ax=ax)
                ax.set_xlabel('IDGroup')
                continue
            else:
                plot_data = data.value_counts(sort=False)
            if col != 'SRA':
                plot_data.index = plot_data.index.astype(int)
            plot_data.plot(kind='bar', ax=ax, logy=False, alpha=alpha, linewidth=0)
            autolabel(ax)
            ax.set_xlabel(col)
            ax.set_ylabel('Count')
            continue

        if hasattr(data, 'columns') and len(data.columns) == 2:
            outer, inner = data.columns
            data[outer].plot(kind='hist', ax=ax, logy=True, alpha=.4,
                                  color=color, linewidth=0, label=outer)
            data[inner].plot(kind='hist', ax=ax, logy=True, alpha=.8,
                                  color=color, linewidth=0, label=inner)
            ax.set_xlabel(outer)
            ax.legend()
        else:
            data.plot(kind='hist', ax=ax, logy=True, alpha=alpha, linewidth=0)
            ax.set_xlabel(col)
        min_, max_ = ax.get_ylim()
        ax.set_ylim(0, max_)
    desc = ''
    if description:
        desc += ' in {}'.format(description)
    title = 'Characteristics of the {} Gene Products{}'
    fig.suptitle(title.format(len(df), desc))
    fig.tight_layout(rect=(0, 0, 1, .95))
    return fig, axes
