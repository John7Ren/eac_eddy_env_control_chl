"""Grid and mask helper functions for the EAC eddy chlorophyll-a workflow."""

from __future__ import annotations

import numpy as np
from matplotlib import path
from matplotlib.collections import LineCollection
from scipy.ndimage import label, sum_labels, zoom


def make_2deg_bins(lat_min=-50, lat_max=-20, lon_min=145, lon_max=175, dlat=2, dlon=2):
    """Return latitude and longitude bin edges used for the 2° x 2° EAC grid."""
    lat_bins = np.arange(lat_max, lat_min - dlat, -dlat)
    lon_bins = np.arange(lon_min, lon_max + dlon, dlon)
    return lat_bins, lon_bins


def resample_spatial(data, new_lat_size, new_lon_size):
    """Resample a 3D array (time, lat, lon) to new spatial dimensions."""
    time_size, orig_lat_size, orig_lon_size = data.shape
    lat_zoom = new_lat_size / orig_lat_size
    lon_zoom = new_lon_size / orig_lon_size
    resampled_data = np.empty((time_size, new_lat_size, new_lon_size))
    for t in range(time_size):
        resampled_data[t, :, :] = zoom(data[t, :, :], (lat_zoom, lon_zoom), order=1)
    return resampled_data


def assign_prof2boxes(prof):
    """Assign profile lat/lon values to the 2° x 2° EAC grid."""
    lat_bins, lon_bins = make_2deg_bins()
    lx = len(lon_bins) - 1
    ly = len(lat_bins) - 1
    iy_ar = np.ones_like(prof["n_prof"].data, dtype=float) * np.nan
    ix_ar = np.ones_like(prof["n_prof"].data, dtype=float) * np.nan
    for iy in range(ly):
        for ix in range(lx):
            ilat_d = lat_bins[iy + 1]
            ilat_u = lat_bins[iy]
            ilon_l = lon_bins[ix]
            ilon_r = lon_bins[ix + 1]
            mask_lat = (prof["lat"].data >= ilat_d) & (prof["lat"].data < ilat_u)
            mask_lon = (prof["lon"].data >= ilon_l) & (prof["lon"].data < ilon_r)
            iloc = np.where(mask_lat & mask_lon)[0]
            if np.size(iloc) != 0:
                iy_ar[iloc] = iy
                ix_ar[iloc] = ix
    return iy_ar, ix_ar


def inpolygon(xq, yq, xv, yv):
    """MATLAB-like inpolygon based on matplotlib.path.Path."""
    shape = xq.shape
    xq = xq.reshape(-1)
    yq = yq.reshape(-1)
    xv = xv.reshape(-1)
    yv = yv.reshape(-1)
    q = [(xq[i], yq[i]) for i in range(xq.shape[0])]
    p = path.Path([(xv[i], yv[i]) for i in range(xv.shape[0])])
    return p.contains_points(q).reshape(shape)


def get_patch_edges(x, y, zmask, c, ls, **kwargs):
    """Return a LineCollection for the boundary edges of a masked grid."""
    lw = kwargs.get("lw", 3)
    alpha = kwargs.get("alpha", 1.0)
    zorder = kwargs.get("zorder", 1)
    X, Y = np.meshgrid(x, y)
    edges = []
    for i in range(zmask.shape[0]):
        for j in range(zmask.shape[1]):
            if zmask[i, j]:
                x0, x1 = X[i, j], X[i, j + 1]
                y0, y1 = Y[i, j], Y[i + 1, j]
                if j == 0 or not zmask[i, j - 1]:
                    edges.append([(x0, y0), (x0, y1)])
                if j == zmask.shape[1] - 1 or not zmask[i, j + 1]:
                    edges.append([(x1, y0), (x1, y1)])
                if i == 0 or not zmask[i - 1, j]:
                    edges.append([(x0, y0), (x1, y0)])
                if i == zmask.shape[0] - 1 or not zmask[i + 1, j]:
                    edges.append([(x0, y1), (x1, y1)])
    return LineCollection(edges, colors=c, linestyle=ls, linewidths=lw, alpha=alpha, zorder=zorder)


def delete_isolate_mask(map_array, size_threshold=5):
    """Remove small isolated True regions from a boolean mask."""
    labeled_array, num_features = label(map_array)
    component_sizes = sum_labels(map_array, labeled_array, index=np.arange(1, num_features + 1))
    large_components_mask = np.isin(labeled_array, np.where(component_sizes >= size_threshold)[0] + 1)
    filtered_map = np.where(large_components_mask, map_array, 0)
    return filtered_map.astype(bool)


def fix_mask_hole(map_array, size_threshold=2):
    """Fill small holes in a boolean mask."""
    unmap_array = (map_array == 0).astype(int)
    labeled_array, num_features = label(unmap_array)
    component_sizes = sum_labels(unmap_array, labeled_array, index=np.arange(1, num_features + 1))
    large_components_mask = np.isin(labeled_array, np.where(component_sizes >= size_threshold)[0] + 1)
    filtered_unmap = np.where(large_components_mask, unmap_array, 0)
    filtered_map = (filtered_unmap == 0).astype(int)
    return filtered_map.astype(bool)
