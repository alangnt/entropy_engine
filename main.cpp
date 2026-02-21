#include <iostream>
#include <cmath>

struct Vector3 {
  double x;
  double y;
  double z;
};

// Passing a const ref to the original memory blocks
double calculateDistance(const Vector3& particleA, const Vector3& particleB) {
  // 1. Find the difference in each axis
  double dx = particleB.x - particleA.x;
  double dy = particleB.y - particleA.y;
  double dz = particleB.z - particleA.z;

  // 2. Square the differences and find the square root
  double distance = std::sqrt((dx * dx) + (dy * dy) + (dz * dz));

  // 3. Return the final number
  return distance;
}

int main() {
  // Earth
  Vector3 earth;
  earth.x = 0.0;
  earth.y = 0.0;
  earth.z = 0.0;

  // Simple Satellite
  Vector3 satellite;
  satellite.x = 100.0;
  satellite.y = 0.0;
  satellite.z = 0.0;

  double distanceToSatellite = calculateDistance(earth, satellite);

  std::cout << "The distance is: " << distanceToSatellite << std::endl;

  return 0;
}