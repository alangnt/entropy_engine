# Entropy Engine

A 3D N-Body gravitational simulation that models the orbital dynamics of 1000 celestial bodies using the **Barnes-Hut algorithm** for efficient O(N log N) force computation, with OpenMP multithreading and real-time 3D visualization.

## How It Works

The simulation initializes 1000 Earth-mass particles at random positions in a 400-million-unit cube, then iteratively computes gravitational interactions over 2.36 million time steps. An **octree** spatially partitions the particles so that distant clusters can be approximated as single massive bodies, reducing the per-step complexity from O(N²) to O(N log N).

Each time step:
1. Build an octree from all particle positions
2. Compute gravitational forces via Barnes-Hut tree traversal (parallelized with OpenMP)
3. Update accelerations, velocities, and positions using Euler integration
4. Record particle states every 3600 steps to CSV

A softening parameter (`ε = 100,000`) prevents force singularities at close range.

## Building & Running

### Prerequisites

- **C++ compiler** with C++17 support (e.g. `g++`)
- **OpenMP** - on macOS, install via `brew install libomp`
- **Python 3** with packages listed in `requirements.txt`

### Compile

```bash
g++ -O3 -Xpreprocessor -fopenmp \
    -I/opt/homebrew/opt/libomp/include \
    -L/opt/homebrew/opt/libomp/lib -lomp \
    main.cpp -o universe
```

### Run the simulation

```bash
./universe
```

This produces `orbit.csv` containing trajectory data for all particles.

### Visualize

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python visualize.py
```

Opens an interactive 3D animation of the orbital trajectories using matplotlib and PyQt6.

## Project Structure

```
├── main.cpp            # Core simulation engine (Barnes-Hut + OpenMP)
├── visualize.py        # 3D scatter plot animation
├── requirements.txt    # Python dependencies (PyQt6, pandas, matplotlib)
└── INFOS.md            # Physics formulas and notes
```

## Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `G` | 6.674×10⁻¹¹ | Gravitational constant |
| `DT` | 1.0 s | Integration time step |
| `THETA` | 0.5 | Barnes-Hut opening angle criterion |
| `EPSILON` | 100,000 | Softening parameter |
| Particles | 1,000 | Number of simulated bodies |
| Steps | 2,360,000 | Total simulation iterations |
