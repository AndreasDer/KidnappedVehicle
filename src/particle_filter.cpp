/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
    std::default_random_engine gen;
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
   // This lines create a normal (Gaussian) distribution for x,y and theta
    std::normal_distribution<double> dist_x(x, std[0]);
    std::normal_distribution<double> dist_y(y, std[1]);
    std::normal_distribution<double> dist_theta(theta, std[2]);
  num_particles = 10;  // TODO: Set the number of particles
  for (int i = 0; i < num_particles; i++) {
      Particle p = Particle();
      p.id = i;
      p.weight = 1;
      p.x= dist_x(gen);
      p.y = dist_y(gen);
      p.theta = dist_theta(gen);
      particles.push_back(p);
  }
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {

    std::default_random_engine gen;
    double x, y, theta;
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
    for (int i = 0; i < num_particles;i++) {
        Particle p = particles[i];
        theta = p.theta + yaw_rate * delta_t;
        x = p.x + velocity / yaw_rate*(sin(theta)-sin(p.theta));
        y = p.y + velocity / yaw_rate * (cos(theta) - cos(p.theta));
        std::normal_distribution<double> dist_x(x, std_pos[0]);
        std::normal_distribution<double> dist_y(y, std_pos[1]);
        std::normal_distribution<double> dist_theta(theta, std_pos[2]);
        p.theta = dist_theta(gen);
        p.x = dist_x(gen);
        p.y = dist_y(gen);
        particles[i] = p;
    }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations,
                                     Particle particle) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
    for (int i = 0; i < observations.size(); i++)
    {
        double final_dist = 1000.0;
        int id = -1;
        double x = 0.0;
        double y = 0.0;
        for (int j = 0; j < predicted.size(); j++) {
            double distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
            if (distance < final_dist) {
                final_dist = distance;
                id = predicted[j].id;
                x = observations[i].x;
                y = observations[i].y;
            }
        }
        particle.associations.push_back(id);
        particle.sense_x.push_back(x);
        particle.sense_y.push_back(y);
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
    for (int p = 0; p < particles.size(); p++) {
        vector<LandmarkObs> landmarksInSight;
        vector<LandmarkObs> translatedObs;
        //std::cout << "Landmarks in sight:" << std::endl;
        for (int l = 0; l < map_landmarks.landmark_list.size(); l++) {
            if (dist(particles[p].x, particles[p].y, map_landmarks.landmark_list[l].x_f, map_landmarks.landmark_list[l].y_f) <= sensor_range) {
                LandmarkObs map_landmark;
                map_landmark.id = map_landmarks.landmark_list[l].id_i;
                map_landmark.x = map_landmarks.landmark_list[l].x_f;
                map_landmark.y = map_landmarks.landmark_list[l].y_f;
                //std::cout << map_landmark.id << ", " << map_landmark.x << "," << map_landmark.y << std::endl;
                landmarksInSight.push_back(map_landmark);
            }
        }
        for (int o = 0; o < observations.size(); o++) {
            // transform to map x coordinate
            double x_map;
            x_map = particles[p].x + (cos(particles[p].theta) * observations[o].x) - (sin(particles[p].theta) * observations[o].y);

            // transform to map y coordinate
            double y_map;
            y_map = particles[p].x + (sin(particles[p].theta) * observations[o].x) + (cos(particles[p].theta) * observations[o].y);
            LandmarkObs translatedObservation;
            translatedObservation.id = 0;
            translatedObservation.x = x_map;
            translatedObservation.y = y_map;
            translatedObs.push_back(translatedObservation);
        }
        dataAssociation(landmarksInSight, translatedObs,particles[p]);
        double particle_weight = 0.0;
        for (int i = 0; i < particles[p].associations.size(); i++) {
            for (int j = 0; j < landmarksInSight.size(); j++) {
                if (landmarksInSight[j].id == particles[p].associations[i]) {
                    if (particle_weight == 0.0) {
                        particle_weight = multiv_prob(std_landmark[0], std_landmark[1], particles[p].sense_x[j], particles[p].sense_y[j], landmarksInSight[j].x, landmarksInSight[j].x);
                    }
                    else {
                        particle_weight *= multiv_prob(std_landmark[0], std_landmark[1], particles[p].sense_x[j], particles[p].sense_y[j], landmarksInSight[j].x, landmarksInSight[j].x);
                    }
                }
            }
        }
        /*std::cout << "Observations:" << std::endl;
        for (int i = 0; i < translatedObs.size(); i++) {
            std::cout << "ID: " << translatedObs[i].id << ", "<< translatedObs[i].x << ","<< translatedObs[i].y<<std::endl;
        }*/

    }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

}

double ParticleFilter::multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs,
    double mu_x, double mu_y) {
    // calculate normalization term
    double gauss_norm;
    gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

    // calculate exponent
    double exponent;
    exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
        + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));

    // calculate weight using normalization terms and exponent
    double weight;
    weight = gauss_norm * exp(-exponent);

    return weight;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

double ParticleFilter::getDistance(LandmarkObs landmark, LandmarkObs observation) {
    return std::sqrt(std::pow(landmark.x - observation.x, 2) + std::pow(landmark.y - observation.y, 2));
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}