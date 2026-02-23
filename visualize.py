import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation

print("Loading universe data...")
df = pd.read_csv('orbit.csv')
df_unique_steps = df['Step'].unique()
print(f"Loaded {len(df_unique_steps)} frames of animation.")

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')

bound = 200000000.0
ax.set_xlim([-bound, bound])
ax.set_ylim([-bound, bound])

z_bound = 20000000.0
ax.set_zlim([-z_bound, z_bound])

ax.set_facecolor('black')
ax.set_box_aspect((1, 1, 0.15))
fig.patch.set_facecolor('black')
ax.grid(False)
ax.axis('off')

scatter = ax.scatter([], [], [], c='white', marker='o', s=1.5)

def update(frame_index):
    current_step = df_unique_steps[frame_index]
    filtered_df = df[df['Step'] == current_step]
    
    scatter._offsets3d = (filtered_df['X'], filtered_df['Y'], filtered_df['Z'])
    return scatter,

print("Starting simulation playback...")
# interval=30 means 30 milliseconds between frame
ani = animation.FuncAnimation(fig, update, frames=len(df_unique_steps), interval=30, blit=False)

plt.show()