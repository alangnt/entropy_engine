#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>

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

  void free() {
    if (northWest != nullptr) {
      northWest->free();
      northEast->free();
      southWest->free();
      southEast->free();
    }

    delete northWest;
    delete northEast;
    delete southWest;
    delete southEast;
  }

  bool insert(Particle* p) {

    // 1. First we check if the current Quadtree contains the particle
    if (!boundary.contains(p)) {
      return false;
    }

    // 2. Then, we update the Quadtree
    update(p);

    // 3. We check if no children and no particle first
    if (northWest == nullptr && particle == nullptr) {
      particle = p;
      return true;
    }

    // 4. We check if no children and particle and push the old particle down
    if (northWest == nullptr && particle != nullptr) {
      subdivide();
      if (northWest->insert(particle)) { /* Success, do nothing else */ }
      else if (northEast->insert(particle)) { }
      else if (southWest->insert(particle)) { }
      else if (southEast->insert(particle)) { }
      particle = nullptr;
    }

    // 5. Push the new particle down
    if (northWest->insert(p)) {
      return true;
    }
    if (northEast->insert(p)) {
      return true;
    }
    if (southWest->insert(p)) {
      return true;
    }
    if (southEast->insert(p)) {
      return true;
    }

    if (p == particle) {}

    return false;
  }

  void update(Particle* p) {
    // Insert the new particle and calculate new center of mass and total mass
    centerOfMass.x = ((totalMass * centerOfMass.x) + (p->mass * p->position.x)) / (totalMass + p->mass);
    centerOfMass.y = ((totalMass * centerOfMass.y) + (p->mass * p->position.y)) / (totalMass + p->mass);
    centerOfMass.z = ((totalMass * centerOfMass.z) + (p->mass * p->position.z)) / (totalMass + p->mass);

    totalMass = totalMass + p->mass;
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

// Theta (MAC ratio)
const double THETA = 0.5;

// Softening parameter
const double EPSILON = 100000.0;

double calculateDistance(const Vector3& pA, const Vector3& pB) {
  double dx = pB.x - pA.x;
  double dy = pB.y - pA.y;
  double dz = pB.z - pA.z;

  // d = sqrt(dx^2 + dy^2 + dz^2)
  double distance = std::sqrt((dx * dx) + (dy * dy) + (dz * dz));

  return distance;
}

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

  double epsilonSquared = EPSILON * EPSILON;

  // 4. Newton's Formula (F = G * (m1 * m2) / r^2)
  double force = GRAVITATIONAL_CONSTANT * (pA.mass * pB.mass) / (distanceSquared + epsilonSquared);

  // 5. Find distance between pA and pB
  double distance = std::sqrt(distanceSquared);

  // 6. Normalize the vector and multiply by the force magnitude
  Vector3 forceVector;
  forceVector.x = (dx / distance) * force;
  forceVector.y = (dy / distance) * force;
  forceVector.z = (dz / distance) * force;

  return forceVector;
}

Vector3 calculateTreeForce(Particle* p, QuadtreeNode* node) {
  Vector3 totalForce = {0.0, 0.0, 0.0};

  if (node->totalMass == 0.0) {
    return totalForce;
  }

  if (node->northWest == nullptr && node->particle != nullptr) {
    if (p == node->particle) {
      return totalForce;
    }
    return calculateGravitationalForceVector(*p, *(node->particle));
  }

  if (node->northWest != nullptr) {
    // find d
    double distance = calculateDistance(p->position, node->centerOfMass);

    // find s
    double size = node->boundary.halfDimension * 2;
    
    // check if s/d is less than theta (0.5)
    // if true, far away
    if (size / distance < THETA) {
      Particle superParticle;
      superParticle.position = node->centerOfMass;
      superParticle.mass = node->totalMass;
      Vector3 force = calculateGravitationalForceVector(*p, superParticle);

      totalForce = force;
    } else {
      Vector3 nwForce = calculateTreeForce(p, node->northWest);
      Vector3 neForce = calculateTreeForce(p, node->northEast);
      Vector3 swForce = calculateTreeForce(p, node->southWest);
      Vector3 seForce = calculateTreeForce(p, node->southEast);

      double totalForceX = nwForce.x + neForce.x + swForce.x + seForce.x;
      double totalForceY = nwForce.y + neForce.y + swForce.y + seForce.y;
      double totalForceZ = nwForce.z + neForce.z + swForce.z + seForce.z;

      totalForce = {totalForceX, totalForceY, totalForceZ};
    }
  }

  return totalForce;
}

double randomDouble(double min, double max) {
  double fraction = (double)rand() / RAND_MAX;

  return min + fraction * (max - min);
}

int main() {
  srand(time(NULL));

  std::vector<Particle> universe;

  // Generate 100 planets
  for (int i = 0; i < 100; i++) {
    Particle planet;
    double randomPositionX = randomDouble(-200000000.0, 200000000.0);
    double randomPositionY = randomDouble(-200000000.0, 200000000.0);

    planet.position.x = randomPositionX; planet.position.y = randomPositionY;
    planet.position.z = 0.0;
    planet.mass = 5.972e24;
    planet.velocity.x = 0.0; planet.velocity.y = 0.0; planet.velocity.z = 0.0;

    universe.push_back(planet);
  }

  std::ofstream trajectoryFile("orbit.csv");
  trajectoryFile << "Step,X,Y,Z\n";

  int step = 0;
  while (step < 2360000) {

    std::vector<Vector3> totalForces(universe.size());
    for (int i = 0; i < totalForces.size(); i++) {
      totalForces[i].x = 0.0; totalForces[i].y = 0.0; totalForces[i].z = 0.0;
    }

    BoundingBox baseBoundary;
    baseBoundary.x = 0.0; baseBoundary.y = 0.0;
    baseBoundary.halfDimension = 400000000.0;

    QuadtreeNode root(baseBoundary);

    for (int i = 0; i < universe.size(); i++) {
      root.insert(&universe[i]);
    }

    for (int i = 0; i < universe.size(); i++) {
      totalForces[i] = calculateTreeForce(&universe[i], &root);
    }

    // free the memory
    root.free();

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
    if (step % 3600 == 0) {
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