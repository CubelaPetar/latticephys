"""Plot phonon band structures with matplotlib."""

import matplotlib.pyplot as plt


def plot_bands(distance, frequencies, corners, labels, ax=None):
    """Branches against the distance along a k-path, with the k-points as ticks."""
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4))

    ax.plot(distance, frequencies, color="tab:blue")
    for corner in corners:
        ax.axvline(corner, color="grey", lw=0.5)
    ax.axhline(0, color="grey", lw=0.5)
    ax.set_xticks(corners, labels)
    ax.set_xlim(corners[0], corners[-1])
    return ax
