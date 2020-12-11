""" Script for making simple visualisations to compare with paper An et al. (2012) """

import os
import sys
import time

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from matplotlib.animation import FuncAnimation

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.regrid as fr
import functions.visualisation as fv


### Parameters
case = "00"
make_ani = True
show_figs = False


### Set paths
output_dir = os.path.dirname(os.path.realpath(__file__)) + f"\\output_repr_{case}\\"
figure_dir = os.path.dirname(os.path.realpath(
    __file__)) + f"\\figures_repr_{case}\\"

os.makedirs(figure_dir, exist_ok=True)

filename_output = output_dir + "FlowFM_map.nc"
filename_processed = output_dir + "processed_data.nc"

print(f"Making figures for case {case}")
print(f"\t * Folder with data is \n\t'{output_dir}'")
print(f"\t * File with raw data is \n\t'{filename_output}'")
print(f"\t * File with processed data is \n\t'{filename_processed}'")
print(f"\t * Folder that contains figures is \n\t'{figure_dir}'")

# Convert data if it didn't happen yet
if not os.path.exists(filename_processed):
    print(f"File '{filename_processed}' is not found. Processing data.")
    fr.extract_data(filename_output, savename=filename_processed)


### Get data
print("Read data")
data = xr.open_dataset(filename_processed, chunks={"t": 3})
x = data["x"]
y = data["y"]
t = data["t"]

b = data["b"]
wl = data["wl"]
u = data["u"]
v = data["v"]
p = data["p"]

# offset x-coordinate for nicer plotting
x = x - x.min()


print(f"Creating figures in \n'{figure_dir}'")
### Figure bathymetry
print("Figure: bathymetry")
plt.figure(figsize=(5, 5), dpi=150)

plt.contourf(x/1000, y/1000, b)
plt.colorbar()

plt.title("Water depth [m]")
plt.xlabel("$x$ [km]")
plt.ylabel("$y$ [km]")

plt.axhline(color="black", linewidth=1)
plt.axvline(color="black", linewidth=1)

# plt.axhline(-7.56e3, linestyle="--", color="gray", linewidth=1)
# plt.axvline(10, linestyle="--", color="gray", linewidth=1)
# plt.axvline(100, linestyle="--", color="gray", linewidth=1)
# plt.axvline(200, linestyle="--", color="gray", linewidth=1)
    
plt.savefig(figure_dir + "bathymetry_contour.jpg", bbox_inches="tight")


### Figure pressure
print("Figure: pressure")
plot_times = [4e4, 8e4, 12e4, 16e4]

fig, ax = plt.subplots(2, 2, sharex=True, sharey=True)
fig.set_size_inches(8, 8)
fig.set_dpi(150)
fig.suptitle("Pressure distribution [Pa]")
fig.set_tight_layout(True)

for i in range(4):
    _ax = ax[i//2, i%2]
    im = _ax.contourf(x/1000, y/1000, p.interp(t=fv.to_timestr(plot_times[i])), vmin=0, vmax=2000)
    _ax.set_xlim([0, 300])
    _ax.set_ylim([0, -y.max()/1000])
    _ax.set_title(f"$t = {plot_times[i] : 0.0f}$s")
    
    if i//2: _ax.set_xlabel("$x$ [km]")
    if not i%2: _ax.set_ylabel("$y$ [km]")
    
    fig.colorbar(im, ax=_ax)
    
#     _ax.axhline(-7.56e3, linestyle="--", color="gray", linewidth=1)
#     _ax.axvline(10, linestyle="--", color="gray", linewidth=1)
#     _ax.axvline(100, linestyle="--", color="gray", linewidth=1)
#     _ax.axvline(200, linestyle="--", color="gray", linewidth=1)
    
fig.savefig(figure_dir + "pressure_contours.jpg", bbox_inches="tight")


### Figure sse contours
print("Figure: sse contours")
# figure 1
plot_times = [4e4, 8e4, 12e4, 16e4]

fig, ax = plt.subplots(2, 2, sharex=True, sharey=True)
fig.set_size_inches(8, 8)
fig.set_dpi(150)
fig.suptitle("Sea Surface Elevation [m]")
for i in range(4):
    _ax = ax[i//2, i%2]
    im = _ax.contourf(x/1000, y/1000, wl.interp(t=fv.to_timestr(plot_times[i])), vmin=-0.8, vmax=0.8)
    _ax.set_xlim([0, 300])
    _ax.set_ylim([0, -y.max()/1000])
    _ax.set_title(f"$t = {plot_times[i] : 0.0f}$s")
    
    if i//2: _ax.set_xlabel("$x$ [km]")
    if not i%2: _ax.set_ylabel("$y$ [km]")
    
    fig.colorbar(im, ax=_ax)
    
#     _ax.axhline(-7.56e3, linestyle="--", color="gray", linewidth=1)
#     _ax.axvline(10, linestyle="--", color="gray", linewidth=1)
#     _ax.axvline(100, linestyle="--", color="gray", linewidth=1)
#     _ax.axvline(200, linestyle="--", color="gray", linewidth=1)
    
fig.savefig(figure_dir + "sse_contours.jpg", bbox_inches="tight")


### Figure alongshore profiles sse
# figure 2
print("Figure: sse alongshore profile")
plot_times = [4e4, 8e4, 12e4, 16e4]
plot_ylims = np.array([[0, 4], [1, 6], [2, 8], [3, 10]]) * -1e3

fig, ax = plt.subplots(2, 2, sharex=False, sharey=False)
fig.set_size_inches(8, 8)
fig.set_dpi(150)
fig.suptitle("Along-shore profile of Sea Surface Elevation [m]")
for i in range(4):
    _ax = ax[i//2, i%2]
    im = _ax.plot(y/1000, wl.interp(t=fv.to_timestr(plot_times[i]), x=10e3))
    _ax.set_xlim(plot_ylims[i])
    _ax.set_ylim([-0.65, 0.65])
    _ax.set_title(f"$t = {plot_times[i] : 0.0f}$s")
    
    if i//2: _ax.set_xlabel("$x$ [km]")
    if not i%2: _ax.set_ylabel("$SSE$ [m]")
        
    _ax.axvline(-7.56e3, linestyle="--", color="gray", linewidth=1)
    _ax.axhline(color="black", linewidth=1)
    _ax.axvline(color="black", linewidth=1)
    
fig.savefig(figure_dir + "sse_along.jpg", bbox_inches="tight")


### Figure cross-shore profile sse
# figure 3
print("Figure: sse cross-shore profile")
y_slices = np.array([7.56]) * -1e6

plt.figure(figsize=(5, 5), dpi=150)
for i in range(len(y_slices)):
    plt.plot(x/1000, wl.interp(t=fv.to_timestr(1.6e5), y=y_slices[i]), label=f"$y={y_slices[i]/1000 : 0.0f}$ km")
plt.legend()
plt.axhline(color="black", linewidth=1)
plt.axvline(color="black", linewidth=1)
plt.title("Cross-shore profile of Sea Surface Elevation [m]")
plt.xlabel("$x$ [km]")
plt.ylabel("$SSE$ [m]")
plt.xlim([0, 600])
plt.ylim([0, 0.7])

# plt.axvline(10, linestyle="--", color="gray", linewidth=1)
# plt.axvline(100, linestyle="--", color="gray", linewidth=1)
# plt.axvline(200, linestyle="--", color="gray", linewidth=1)
    
plt.savefig(figure_dir + "sse_cross.jpg", bbox_inches="tight")


### Figure alongshore profiles
# figure
print("Figure: sse more alongshore profiles")
x_slices = np.array([10, 100, 200]) * 1e3
linestyles = ["-", "--", "--"]

plt.figure(figsize=(5, 5), dpi=150)
plt.title("Along-shore profile of Sea Surface Elevation [m]")
for i in range(len(x_slices)):
    plt.plot(
        y/1000, 
        wl.interp(t=fv.to_timestr(1.6e5), x=x_slices[i]), 
        label=f"$x={x_slices[i]/1000 : 0.0f}$ km",
        linestyle=linestyles[i]
    )
plt.legend()
plt.axhline(color="black", linewidth=1)
plt.axvline(color="black", linewidth=1)
plt.xlim([0, -y.max()/1000])
plt.axvline(-7.56e3, linestyle="--", color="gray", linewidth=1)
plt.xlabel("$y$ [km]")
plt.ylabel("$SSE$ [m]")
    
fig.savefig(figure_dir + "sse_along_2.jpg", bbox_inches="tight")


### Animation
print("Animation: sse + pressure")
fig, ax = plt.subplots(2, 1, sharex=True)
fig.set_size_inches(14.4, 7.2)
fig.set_dpi(100)
fig.set_tight_layout(True)

##
lines = np.zeros(2, dtype=np.object)
texts = np.zeros(1, dtype=np.object)

lines[0] = ax[0].plot(y/1000, wl.isel(t=0).interp(x=10000))[0]
lines[1] = ax[1].plot(y/1000, p.isel(t=0).interp(x=10000), color="tab:red")[0]
texts[0] = ax[0].set_title(f"$t = {t.isel(t=0).values.tolist() / 1e9 / 3600 : 5.1f}$ hours since start")

##
def initfig():
    for i in range(2):
        ax[i].grid()
        ax[i].set_xlim(0, -y.max()/1000)
        ax[i].axhline(0, color="black", linewidth=1)
    ax[0].set_ylabel("$SSE$ [m]")
    ax[1].set_ylabel("$Pressure$ [Pa]")
    ax[1].set_xlabel("$y$ [km]")
    ax[0].set_ylim([-0.7, 0.7])
    ax[1].set_ylim([0, 2000])
    return tuple(lines.flatten()) + tuple(texts.flatten())

initfig()

##
def update(i):
    lines[0].set_ydata(wl.isel(t=i).interp(x=10000))
    lines[1].set_ydata(p.isel(t=i).interp(x=10000))
    texts[0].set_text(f"$t = {t.isel(t=i).values.tolist() / 1e9 / 3600 : 5.1f}$ hours since start")

    return tuple(lines.flatten()) + tuple(texts.flatten())

##
frames = (np.arange(t.size)).astype(np.int)
anim = FuncAnimation(
    fig,
    update,
    init_func=initfig,
    blit=True,
    frames=frames,
    interval=1000/20
    )

if make_ani:
    print("\tAnimating")
    t0 = time.perf_counter()
    anim.save(figure_dir + "animation_sse_p.mp4")
    t1 = time.perf_counter()
    print(f"\tFinished animation in {t1 - t0 : 0.1f} seconds")
else:
    print("\tNo animation")

### End
if show_figs: 
    plt.show()
