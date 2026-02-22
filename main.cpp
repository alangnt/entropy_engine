#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

struct Vector3 {
  double x;
  double y;
  double z;
};

struct Particle {
  Vector3 position;
  Vector3 acceleration;
  Vector3 velocity;
  double mass;
};

struct BoundingBox {
  double x;
  double y;
  double halfDimension;

  bool contains(Particle* p) {
    double rightBoundary = x + halfDimension;
    double leftBoundary = x - halfDimension;
    double topBoundary = y + halfDimension;
    double bottomBoundary = y - halfDimension;

    if (p->position.x >= leftBoundary && p->position.x < rightBoundary && p->position.y >= bottomBoundary && p->position.y < topBoundary) {
      return true;
    }
    return false;
  }
};

struct QuadtreeNode {
  BoundingBox boundary;
  Particle* particle;
  
  double totalMass;
  Vector3 centerOfMass;

  QuadtreeNode* northWest;
  QuadtreeNode* northEast;
  QuadtreeNode* southWest;
  QuadtreeNode* southEast;

  QuadtreeNode(BoundingBox b) {
    boundary = b;
    particle = nullptr;
    totalMass = 0.0;
    centerOfMass.x = 0.0; centerOfMass.y = 0.0; centerOfMass.z = 0.0;

    northWest = nullptr; northEast = nullptr;
    southWest = nullptr; southEast = nullptr;
  }

  void subdivide() {
    // We cut the size of the box in half for the children
    double newHalf = boundary.halfDimension / 2.0;

    // 1. Create the NorthWest box
    BoundingBox nwBoundary;
    nwBoundary.x = boundary.x - newHalf; // Move Left
    nwBoundary.y = boundary.y + newHalf; // Move Up
    nwBoundary.halfDimension = newHalf;
    // Claim the RAM on the Heap!
    northWest = new QuadtreeNode(nwBoundary);

    // 2. Create the NorthEast box
    BoundingBox neBoundary;
    neBoundary.x = boundary.x + newHalf; // Move Right
    neBoundary.y = boundary.y + newHalf; // Move Up
    neBoundary.halfDimension = newHalf;
    northEast = new QuadtreeNode(neBoundary);

    // 3. Create the SouthWest box
    BoundingBox swBoundary;
    swBoundary.x = boundary.x - newHalf; // Move Left
    swBoundary.y = boundary.y - newHalf; // Move Down
    swBoundary.halfDimension = newHalf;
    southWest = new QuadtreeNode(swBoundary);

    // 4. Create the SouthEast box
    BoundingBox seBoundary;
    seBoundary.x = boundary.x + newHalf; // Move Right
    seBoundary.y = boundary.y - newHalf; // Move Down
    seBoundary.halfDimension = newHalf;
    southEast = new QuadtreeNode(seBoundary);
  }
};

// Define G as a constant double
const double GRAVITATIONAL_CONSTANT = 6.67430e-11;

// Define dt as a constant double
const double DT = 1.0;

// Passing a const ref to the original memory blocks
Vector3 calculateGravitationalForceVector(const Particle& pA, const Particle& pB) {
  // 1. Find the difference in each axis
  double dx = pB.position.x - pA.position.x;
  double dy = pB.position.y - pA.position.y;
  double dz = pB.position.z - pA.position.z;

  // 2. Calculate the distance squared
  // (no need to square because we need to use r^2 so they cancel each other)
  double distanceSquared = (dx * dx) + (dy * dy) + (dz * dz);

  // 3. Safety check to avoid dividing by zero
  // plus two particles can't be at the exact same point
  if (distanceSquared == 0.0) {
    Vector3 zeroForce;
    zeroForce.x = 0.0; zeroForce.y = 0.0; zeroForce.z = 0.0;
    return zeroForce;
  }

  // 4. Newton's Formula (F = G * (m1 * m2) / r^2)
  double force = GRAVITATIONAL_CONSTANT * (pA.mass * pB.mass) / distanceSquared;

  // 5. Find distance between pA and pB
  double distance = std::sqrt(distanceSquared);

  // 6. Normalize the vector and multiply by the force magnitude
  Vector3 forceVector;
  forceVector.x = (dx / distance) * force;
  forceVector.y = (dy / distance) * force;
  forceVector.z = (dz / distance) * force;

  return forceVector;
}

int main() {
  std::vector<Particle> universe;

  // Earth
  Particle earth;
  earth.position.x = 0.0; earth.position.y = 0.0; earth.position.z = 0.0;
  earth.velocity.x = 0.0; earth.velocity.y = 0.0; earth.velocity.z = 0.0;
  earth.mass = 5.972e24;
  universe.push_back(earth);

  // Moon
  Particle moon;
  moon.position.x = 384400000.0; moon.position.y = 0.0; moon.position.z = 0.0;
  moon.velocity.x = 0.0; moon.velocity.y = 1022.0; moon.velocity.z = 0.0;
  moon.mass = 7.34767309e22;
  universe.push_back(moon);

  // LEO Satellite
  Particle satellite;
  satellite.position.x = 6771000.0; satellite.position.y = 0.0; satellite.position.z = 0.0;
  satellite.velocity.x = 0.0; satellite.velocity.y = 7660.0; satellite.velocity.z = 0.0;
  satellite.mass = 420.0;
  universe.push_back(satellite);

  std::ofstream trajectoryFile("orbit.csv");
  trajectoryFile << "Step,X,Y,Z\n";

  int step = 0;
  while (step < 2360000) {

    std::vector<Vector3> totalForces(universe.size());
    for (int i = 0; i < totalForces.size(); i++) {
      totalForces[i].x = 0.0; totalForces[i].y = 0.0; totalForces[i].z = 0.0;
    }

    for (int i = 0; i < universe.size(); i++) {
      for (int j = 0; j < universe.size(); j++) {

        if (i == j) { 
          continue; 
        }

        Vector3 force = calculateGravitationalForceVector(universe[i], universe[j]);

        totalForces[i].x = totalForces[i].x + force.x;
        totalForces[i].y = totalForces[i].y + force.y;
        totalForces[i].z = totalForces[i].z + force.z;
      }
    }

    for (int i = 0; i < universe.size(); i++) {
            
      // Calculate Acceleration (a = F / m)
      double accelX = totalForces[i].x / universe[i].mass;
      double accelY = totalForces[i].y / universe[i].mass;
      double accelZ = totalForces[i].z / universe[i].mass;

      // Update Velocity (v = v + a * dt)
      universe[i].velocity.x = universe[i].velocity.x + (accelX * DT);
      universe[i].velocity.y = universe[i].velocity.y + (accelY * DT);
      universe[i].velocity.z = universe[i].velocity.z + (accelZ * DT);

      // Update Position (p = p + v * dt)
      universe[i].position.x = universe[i].position.x + (universe[i].velocity.x * DT);
      universe[i].position.y = universe[i].position.y + (universe[i].velocity.y * DT);
      universe[i].position.z = universe[i].position.z + (universe[i].velocity.z * DT);
    }

    // Print the position every hour of the simulation time
    if (step % 60 == 0) {
      for (int i = 0; i < universe.size(); i++) {
        trajectoryFile << step << "," 
          << i << "," // 'i' is the ParticleID (0 for Earth, 1 for Moon)
          << universe[i].position.x << "," 
          << universe[i].position.y << "," 
          << universe[i].position.z << "\n";
      }  
    }

    step++;
  }

  trajectoryFile.close();
  std::cout << "Simulation complete. Data saved to orbit.csv" << std::endl;

  return 0;
}