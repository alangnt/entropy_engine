#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <omp.h>
#include <random>

// Define G as a constant double
const double GRAVITATIONAL_CONSTANT = 6.67430e-11;

// Define dt as a constant double
const double DT = 1.0;

// Theta (MAC ratio)
const double THETA = 0.5;

// Softening parameter
const double EPSILON = 100000.0;

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
  double z;
  double halfDimension;

  bool contains(Particle* p) {
    double rightBoundary = x + halfDimension;
    double leftBoundary = x - halfDimension;
    double topBoundary = y + halfDimension;
    double bottomBoundary = y - halfDimension;
    double frontBoundary = z + halfDimension;
    double backBoundary = z - halfDimension;

    if (p->position.x >= leftBoundary && p->position.x < rightBoundary 
      && p->position.y >= bottomBoundary && p->position.y < topBoundary
      && p->position.z >= backBoundary && p->position.z < frontBoundary) {
      return true;
    }
    return false;
  }
};

struct OctreeNode {
  BoundingBox boundary;
  Particle* particle;

  double totalMass;
  Vector3 centerOfMass;

  OctreeNode* topNorthWest; OctreeNode* bottomNorthWest;
  OctreeNode* topNorthEast; OctreeNode* bottomNorthEast;
  OctreeNode* topSouthWest; OctreeNode* bottomSouthWest;
  OctreeNode* topSouthEast; OctreeNode* bottomSouthEast;

  OctreeNode(BoundingBox b) {
    boundary = b;
    particle = nullptr;
    totalMass = 0.0;
    centerOfMass.x = 0.0; centerOfMass.y = 0.0; centerOfMass.z = 0.0;

    topNorthWest = nullptr; topNorthEast = nullptr;
    topSouthWest = nullptr; topSouthEast = nullptr;
    bottomNorthWest = nullptr; bottomNorthEast = nullptr;
    bottomSouthWest = nullptr; bottomSouthEast = nullptr;
  }

  void free() {
    if (topNorthWest != nullptr) {
      topNorthWest->free(); topNorthEast->free();
      topSouthWest->free(); topSouthEast->free();
      bottomNorthWest->free(); bottomNorthEast->free();
      bottomSouthWest->free(); bottomSouthEast->free();
    }

    delete topNorthWest; delete topNorthEast;
    delete topSouthWest; delete topSouthEast;
    delete bottomNorthWest; delete bottomNorthEast;
    delete bottomSouthWest; delete bottomSouthEast;
  }

  void update(Particle* p) {
    // Insert the new particle and calculate new center of mass and total mass
    centerOfMass.x = ((totalMass * centerOfMass.x) + (p->mass * p->position.x)) / (totalMass + p->mass);
    centerOfMass.y = ((totalMass * centerOfMass.y) + (p->mass * p->position.y)) / (totalMass + p->mass);
    centerOfMass.z = ((totalMass * centerOfMass.z) + (p->mass * p->position.z)) / (totalMass + p->mass);

    totalMass = totalMass + p->mass;
  }

  bool insert(Particle* p) {

    // 1. First we check if the current Quadtree contains the particle
    if (!boundary.contains(p)) {
      return false;
    }

    // 2. Then, we update the Quadtree
    update(p);

    // 3. We check if no children and no particle first
    if (topNorthWest == nullptr && particle == nullptr) {
      particle = p;
      return true;
    }

    // 4. We check if no children and particle and push the old particle down
    if (topNorthWest == nullptr && particle != nullptr) {
      subdivide();
      if (topNorthWest->insert(particle)) { /* Success, do nothing else */ }
      else if (topNorthEast->insert(particle)) { }
      else if (topSouthWest->insert(particle)) { }
      else if (topSouthEast->insert(particle)) { }
      else if (bottomNorthWest->insert(particle)) { }
      else if (bottomNorthEast->insert(particle)) { }
      else if (bottomSouthWest->insert(particle)) { }
      else if (bottomSouthEast->insert(particle)) { }
      particle = nullptr;
    }

    // 5. Push the new particle down
    if (topNorthWest->insert(p)) return true;
    if (topNorthEast->insert(p)) return true;
    if (topSouthWest->insert(p)) return true;
    if (topSouthEast->insert(p)) return true;
    if (bottomNorthWest->insert(p)) return true;
    if (bottomNorthEast->insert(p)) return true;
    if (bottomSouthWest->insert(p)) return true;
    if (bottomSouthEast->insert(p)) return true;

    return false;
  }

  void subdivide() {
    // We cut the size of the box in half for the children
    double newHalf = boundary.halfDimension / 2.0;

    BoundingBox tnw;
    tnw.x = boundary.x - newHalf; tnw.y = boundary.y + newHalf; tnw.z = boundary.z + newHalf;
    tnw.halfDimension = newHalf;
    topNorthWest = new OctreeNode(tnw);

    BoundingBox tne;
    tne.x = boundary.x + newHalf; tne.y = boundary.y + newHalf; tne.z = boundary.z + newHalf;
    tne.halfDimension = newHalf;
    topNorthEast = new OctreeNode(tne);

    BoundingBox tsw;
    tsw.x = boundary.x - newHalf; tsw.y = boundary.y - newHalf; tsw.z = boundary.z + newHalf;
    tsw.halfDimension = newHalf;
    topSouthWest = new OctreeNode(tsw);

    BoundingBox tse;
    tse.x = boundary.x + newHalf; tse.y = boundary.y - newHalf; tse.z = boundary.z + newHalf;
    tse.halfDimension = newHalf;
    topSouthEast = new OctreeNode(tse);

    BoundingBox bnw;
    bnw.x = boundary.x - newHalf; bnw.y = boundary.y + newHalf; bnw.z = boundary.z - newHalf;
    bnw.halfDimension = newHalf;
    bottomNorthWest = new OctreeNode(bnw);

    BoundingBox bne;
    bne.x = boundary.x + newHalf; bne.y = boundary.y + newHalf; bne.z = boundary.z - newHalf;
    bne.halfDimension = newHalf;
    bottomNorthEast = new OctreeNode(bne);

    BoundingBox bsw;
    bsw.x = boundary.x - newHalf; bsw.y = boundary.y - newHalf; bsw.z = boundary.z - newHalf;
    bsw.halfDimension = newHalf;
    bottomSouthWest = new OctreeNode(bsw);

    BoundingBox bse;
    bse.x = boundary.x + newHalf; bse.y = boundary.y - newHalf; bse.z = boundary.z - newHalf;
    bse.halfDimension = newHalf;
    bottomSouthEast = new OctreeNode(bse);
  }
};

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

Vector3 calculateTreeForce(Particle* p, OctreeNode* node) {
  Vector3 totalForce = {0.0, 0.0, 0.0};

  if (node->totalMass == 0.0) {
    return totalForce;
  }

  if (node->topNorthWest == nullptr && node->particle != nullptr) {
    if (p == node->particle) {
      return totalForce;
    }
    return calculateGravitationalForceVector(*p, *(node->particle));
  }

  if (node->topNorthWest != nullptr) {
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
      Vector3 tnwForce = calculateTreeForce(p, node->topNorthWest);
      Vector3 tneForce = calculateTreeForce(p, node->topNorthEast);
      Vector3 tswForce = calculateTreeForce(p, node->topSouthWest);
      Vector3 tseForce = calculateTreeForce(p, node->topSouthEast);
      Vector3 bnwForce = calculateTreeForce(p, node->bottomNorthWest);
      Vector3 bneForce = calculateTreeForce(p, node->bottomNorthEast);
      Vector3 bswForce = calculateTreeForce(p, node->bottomSouthWest);
      Vector3 bseForce = calculateTreeForce(p, node->bottomSouthEast);

      double totalForceX = tnwForce.x + tneForce.x + tswForce.x + tseForce.x + bnwForce.x + bneForce.x + bswForce.x + bseForce.x;
      double totalForceY = tnwForce.y + tneForce.y + tswForce.y + tseForce.y + bnwForce.y + bneForce.y + bswForce.y + bseForce.y;
      double totalForceZ = tnwForce.z + tneForce.z + tswForce.z + tseForce.z + bnwForce.z + bneForce.z + bswForce.z + bseForce.z;

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

  // Generate our first particle
  // the "Black Hole"
  Particle blackHole;
  blackHole.position.x = 0.0; blackHole.position.y = 0.0; blackHole.position.z = 0.0;
  blackHole.velocity.x = 0.0; blackHole.velocity.y = 0.0; blackHole.velocity.z = 0.0;
  blackHole.mass = 8.54e36;
  universe.push_back(blackHole);

  std::random_device rd;
  // Mersenne Twister prime-number
  std::mt19937 gen(rd());
  // then we define the Log-Normal probability curve
  std::lognormal_distribution<double> massDistribution(0.0, 1.0);


  // Generate 2,000 planets
  for (int i = 0; i < 2000; i++) {
    Particle planet;
    double randomPositionX = (rand() % 400000000) - 200000000;
    double randomPositionY = (rand() % 400000000) - 200000000;
    double randomPositionZ = (rand() % 40000000) - 20000000;

    planet.position.x = randomPositionX; 
    planet.position.y = randomPositionY;
    planet.position.z = randomPositionZ;

    double massMultiplier = massDistribution(gen);
    planet.mass = 5.972e24 * massMultiplier;

    // Calculate the exact distance (r) from the Black Hole (0,0,0)
    // Using the Pythagorean theorem: r = sqrt(x^2 + y^2)
    double r = std::sqrt(randomPositionX * randomPositionX + randomPositionY * randomPositionY);

    // Prevent division by zero just in case a planet spawns exactly at 0,0
    if (r == 0) r = 1.0; 

    // Calculate the exact Orbital Velocity (v)
    // Formula: v = sqrt((G * M) / r)
    double blackHoleMass = 8.54e36;
    double v = std::sqrt((GRAVITATIONAL_CONSTANT * blackHoleMass) / r);

    // Calculate the Tangential Unit Vector (-y/r, x/r) and multiply by v
    planet.velocity.x = v * (-randomPositionY / r);
    planet.velocity.y = v * (randomPositionX / r);
    planet.velocity.z = 0.0; // Keep the disk flat

    universe.push_back(planet);
  }

  std::ofstream trajectoryFile("orbit.csv");
  trajectoryFile << "Step,ParticleID,X,Y,Z,Mass\n";

  int step = 0;
  while (step < 2360000) {

    std::vector<Vector3> totalForces(universe.size());
    for (int i = 0; i < totalForces.size(); i++) {
      totalForces[i].x = 0.0; totalForces[i].y = 0.0; totalForces[i].z = 0.0;
    }

    BoundingBox baseBoundary;
    baseBoundary.x = 0.0; baseBoundary.y = 0.0;
    baseBoundary.halfDimension = 400000000.0;

    OctreeNode root(baseBoundary);

    for (int i = 0; i < universe.size(); i++) {
      root.insert(&universe[i]);
    }

    #pragma omp parallel for
    for (int i = 0; i < universe.size(); i++) {
      totalForces[i] = calculateTreeForce(&universe[i], &root);
    }

    // free the memory
    root.free();

    #pragma omp parallel for
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
          << universe[i].position.z << ","
          << universe[i].mass << "\n";
      }  
    }

    step++;
  }

  trajectoryFile.close();
  std::cout << "Simulation complete. Data saved to orbit.csv" << std::endl;

  return 0;
}