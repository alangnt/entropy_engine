import csv
import matplotlib.pyplot as plt

trajectories = {}

with open('orbit.csv', 'r') as file:
  reader = csv.reader(file)
  next(reader)

  for row in reader:
    p_id = int(row[1])
    x = float(row[2])
    y = float(row[3])

    if p_id not in trajectories:
      trajectories[p_id] = {'x': [], 'y': []}

    trajectories[p_id]['x'].append(x)
    trajectories[p_id]['y'].append(y)

plt.figure(figsize=(10, 10))

labels = {0: "Earth", 1: "Moon", 2: "Low Earth Orbit Satellite"}
colors = {0: "green", 1: "blue", 2: "red"}

for p_id, data in trajectories.items():
    
    name = labels.get(p_id, f"Unknown Particle {p_id}")
    color = colors.get(p_id, "black")
    
    plt.plot(data['x'], data['y'], label=name, color=color)

plt.title("N-Body Dynamic Universe (3-Body System)")
plt.xlabel("X Position (m)")
plt.ylabel("Y Position (m)")
plt.legend()
plt.grid(True)
plt.axis('equal')

plt.xlim(-10000000, 10000000)
plt.ylim(-10000000, 10000000)

plt.show()