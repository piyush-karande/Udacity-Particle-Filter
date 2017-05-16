/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
//   x, y, theta and their uncertainties from GPS) and all weights to 1.
// Add random Gaussian noise to each particle.
// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

  // Set number of particles
  num_particles = 50;


  Particle particle;

  // Sample data from guassian distributions around initial GPS coordinates
  default_random_engine gen;
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  for (int i = 0; i < num_particles; i++) {
    particle.id = i;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1.0;

    weights.push_back(1.0);
    particles.push_back(particle);

  }

  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  // TODO: Add measurements to each particle and add random Gaussian noise.
  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/
  double x_pred;
  double y_pred;
  double theta_pred;

  double x0;
  double y0;
  double theta0;

  Particle particle;

  default_random_engine gen;

  for (int i = 0; i < num_particles; i++) {

    particle = particles[i];

    x0 = particle.x;
    y0 = particle.y;
    theta0 = particle.theta;

    if (yaw_rate == 0) {
      x_pred = velocity * delta_t*cos(theta0);
      y_pred = velocity * delta_t*sin(theta0);
      theta_pred = 0.0;
    } else {
      x_pred = velocity * (sin(theta0 + yaw_rate * delta_t) - sin(theta0)) / yaw_rate;
      y_pred = velocity * (-cos(theta0 + yaw_rate * delta_t) + cos(theta0)) / yaw_rate;;
      theta_pred = yaw_rate * delta_t;
    }

    normal_distribution<double> dist_x(x0 + x_pred, std_pos[0]);
    normal_distribution<double> dist_y(y0 + y_pred, std_pos[1]);
    normal_distribution<double> dist_theta(theta0 + theta_pred, std_pos[2]);

    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);

    particles[i] = particle;
  }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
  // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
  //   observed measurement to this particular landmark.
  // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
  //   implement this method and use it as a helper during the updateWeights phase.
  for (int i = 0; i < observations.size(); i++) {

    double min_distance = numeric_limits<double>::max();

    for (int j = 0; j < predicted.size(); j++) {
      LandmarkObs lmk = predicted[j];
      if (dist(observations[i].x, observations[i].y, lmk.x, lmk.y) < min_distance) {
        observations[i].id = lmk.id;
        min_distance = dist(observations[i].x, observations[i].y, lmk.x, lmk.y);
      }
    }


  }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   std::vector<LandmarkObs> observations, Map map_landmarks) {
  // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
  //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
  // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
  //   according to the MAP'S coordinate system. You will need to transform between the two systems.
  //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
  //   The following is a good resource for the theory:
  //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
  //   and the following is a good resource for the actual equation to implement (look at equation
  //   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account
  //   for the fact that the map's y-axis actually points downwards.)
  //   http://planning.cs.uiuc.edu/node99.html


  // Create a vector of all landmarks
  std::vector<LandmarkObs> map_list;
  for (int i = 0; i < map_landmarks.landmark_list.size(); i++) {
    LandmarkObs curr;
    curr.id = map_landmarks.landmark_list[i].id_i;
    curr.x = map_landmarks.landmark_list[i].x_f;
    curr.y = map_landmarks.landmark_list[i].y_f;
    map_list.push_back(curr);
  }

  // Loop through all particles to update weights
  for (int i = 0; i < num_particles; i++) {
    Particle particle = particles[i];
    std::vector<LandmarkObs> tobs;

    // Transform observations
    for (int j = 0; j < observations.size(); j++) {
      LandmarkObs curr_obs;
      curr_obs.x = observations[j].x * cos(particle.theta) - observations[j].y * sin(particle.theta) + particle.x;
      curr_obs.y = observations[j].x * sin(particle.theta) + observations[j].y * cos(particle.theta) + particle.y;
      tobs.push_back(curr_obs);
    }

    // Create a list of landmarks in the range of the sensor
    std::vector<LandmarkObs> map_list_inRange;

    for (int j = 0; j < map_list.size(); j++) {
      LandmarkObs curr = map_list[j];

      if (dist(curr.x, curr.y, particle.x, particle.y) <= sensor_range) {
        map_list_inRange.push_back(curr);
      }

    }

    // Update the ids of transformed observations
    dataAssociation(map_list_inRange, tobs);

    // Update weights using multivariate gaussian distribution
    particle.weight = 1.0;
    for (int j = 0; j < tobs.size(); j++) {
      double sig_x = std_landmark[0];
      double sig_y = std_landmark[1];
      double c = 1 / (2 * M_PI * sig_x * sig_y);
      double x_part = pow((tobs[j].x - map_list[tobs[j].id - 1].x), 2) / (2 * sig_x * sig_x);
      double y_part = pow((tobs[j].y - map_list[tobs[j].id - 1].y), 2) / (2 * sig_y * sig_y);
      particle.weight *= c * exp(-(x_part + y_part));
    }

    weights[i] = particle.weight;

  }

}

void ParticleFilter::resample() {
  // TODO: Resample particles with replacement with probability proportional to their weight.
  // NOTE: You may find std::discrete_distribution helpful here.
  //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  
  default_random_engine gen;
  discrete_distribution<int> d(weights.begin(), weights.end());
  std::vector<Particle> resampled_particles;
  int new_particle_ind;
  Particle new_particle;


  for (int i = 0; i < num_particles; i++) {

    new_particle_ind = d(gen);
    new_particle = particles[new_particle_ind];
    new_particle.weight = 1.0;

    resampled_particles.push_back(new_particle);
  }
  weights.assign(num_particles, 1.0);

  particles = resampled_particles;


}

void ParticleFilter::write(std::string filename) {
  // You don't need to modify this file.
  std::ofstream dataFile;
  dataFile.open(filename, std::ios::app);
  for (int i = 0; i < num_particles; ++i) {
    dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
  }
  dataFile.close();
}
