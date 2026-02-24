import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

print("Loading universe data...")
df = pd.read_csv('orbit.csv')
df_unique_steps = df['Step'].unique()
print(f"Loaded {len(df_unique_steps)} frames of animation.")

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')

bound = 40000000000.0  
ax.set_xlim([-bound, bound])
ax.set_ylim([-bound, bound])

z_bound = 4000000000.0  
ax.set_zlim([-z_bound, z_bound])

ax.set_facecolor('black')
ax.set_box_aspect((1.0, 1.0, 0.15))
fig.patch.set_facecolor('black')
ax.grid(False)
ax.axis('off')

scatter = ax.scatter([], [], [], c=[], cmap='RdYlBu_r', marker='o', s=1.5)

def update(frame_index):
    current_step = df_unique_steps[frame_index]
    filtered_df = df[df['Step'] == current_step]

    x = filtered_df['X'].to_numpy()
    y = filtered_df['Y'].to_numpy()
    z = filtered_df['Z'].to_numpy()
    raw_masses = filtered_df['Mass'].to_numpy()

    scaled_masses = np.log10(raw_masses)

    min_mass = np.min(scaled_masses)
    max_mass = np.max(scaled_masses)

    # value - min / max - min
    normalized_masses = (scaled_masses - min_mass) / (max_mass - min_mass)

    # translate into real ui pixels
    final_sizes = (normalized_masses * 150) + 1.0
    
    scatter._offsets3d = (x, y, z)
    scatter.set_sizes(final_sizes)
    scatter.set_array(normalized_masses)
    
    return scatter,

print("Starting simulation playback...")
# interval=30 means 30 milliseconds between frame
ani = animation.FuncAnimation(fig, update, frames=len(df_unique_steps), interval=30, blit=False)

plt.show()