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
#include <math.h>
#include <limits>

#include "particle_filter.h"

using namespace std;

#define PI 3.14159
#define DEBUG 0
void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    if (!std){
        perror("empty std array!");
        return;
    }
    this->num_particles = 10;
    particles.resize(num_particles);
    weights.resize(num_particles);

    default_random_engine gen;
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);
    cout << ""<<endl;
    for (int i = 0; i< particles.size(); i++){
        particles[i].x =  dist_x(gen);
        particles[i].y =  dist_y(gen);
        particles[i].theta =  dist_theta(gen);
        particles[i].weight = 1;
        weights[i] = particles[i].weight;
    }
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    default_random_engine gen;
    std::normal_distribution<double> dist_x(0.0, std_pos[0]);
    std::normal_distribution<double> dist_y(0.0, std_pos[1]);
    std::normal_distribution<double> dist_theta(0.0, std_pos[2]);

    for (int i = 0; i< particles.size(); i++)
    {
        if(fabs(yaw_rate) < 0.00001)
        {   // avoid nan
            particles[i].x += (velocity * delta_t * cos(particles[i].theta)) +dist_x(gen);
            particles[i].y += (velocity * delta_t * sin(particles[i].theta)) +dist_y(gen);
            particles[i].theta += dist_theta(gen);
        }
        else
        {
            particles[i].x += (velocity/yaw_rate) * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta)) +dist_x(gen);
            particles[i].y += (velocity/yaw_rate) * (cos(particles[i].theta) -  cos(particles[i].theta + yaw_rate*delta_t)) +dist_y(gen);
            particles[i].theta += yaw_rate * delta_t + dist_theta(gen);
        }
    }
}

void ParticleFilter::printout()
{
    cout << "Num of particles: " << num_particles<<endl;
    for (int i =0; i< num_particles;i++)
    {
        cout<<i<<") x: "<<particles[i].x<<" y: "<<particles[i].y<< " theta: "<<particles[i].theta<<endl;
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.


    // THIS METHOD IS NOT CALLED
    for(int i=0; i<observations.size(); i++)
    {
        double min_distance = std::numeric_limits<double>::max();
        unsigned int found_id = -1;
        for(int j=0; j<predicted.size(); j++)
        {
            double range = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);

            if(min_distance> range)
            {
                min_distance = range;
                found_id = predicted[i].id;
            }
        }
        if(found_id < predicted.size())
            observations[i] = predicted[found_id];
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
    // sensor_range: [m]
    // std_landmark: Landmark measurement uncertainty [x [m], y [m]]
    // noisy observations: vector of landmark observations from the vehicle sensor (x,y,id)
    // See Lesson 13: Particle Filter, Section 14 Importance Weight
/*    def measurement_prob(self, measurement):

        # calculates how likely a measurement should be

        prob = 1.0;
        for i in range(len(landmarks)):
            dist = sqrt((self.x - landmarks[i][0]) ** 2 + (self.y - landmarks[i][1]) ** 2)
            prob *= self.Gaussian(dist, self.sense_noise, measurement[i])
        return prob*/

    //1. Transform observations from the vehicle coordinate to map coordinate
    //2.
    for(int i=0; i<num_particles; i++)
    {
        double weight_p = 1.0;
        if (DEBUG)
            printf("P %d Weight: ", i);
        for(int j=0; j< observations.size(); j++)
        {
            // transform observation to map coordinates -> we need a Transformation Matrix 2x3
            // A =[cos(theta) -sin(theta) p.x,
            //     sin(theta)] cos(theta) p.y]
            // since the CS is not rotate, we only need to shift by the location of the particle
            // MatrixXd A = MatrixXd(2, 3);
            //A << cos(theta), -sin(theta), p.x,
            //     sin(theta), cos(theta), p.y;

            double obsx_map = cos(particles[i].theta)*observations[j].x - sin(particles[i].theta)*observations[j].y + particles[i].x;
            double obsy_map = sin(particles[i].theta)*observations[j].x + cos(particles[i].theta)*observations[j].y + particles[i].y;

            double min_distance = std::numeric_limits<double>::max();
            unsigned int found_id = -1;
            double landm_x, landm_y =0;
            for(int k=0; k< map_landmarks.landmark_list.size();k++)
            {
                double range = dist(map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f, obsx_map, obsy_map);
                if(min_distance > range)
                {
                    min_distance = range;
                    found_id = map_landmarks.landmark_list[k].id_i;
                    landm_x = map_landmarks.landmark_list[k].x_f;
                    landm_y = map_landmarks.landmark_list[k].y_f;
                }
            }

            // multivariate Gaussian
            // ignoring the factor since square root consumes more calc time


            double prob = exp(-0.5*( (pow(obsx_map - landm_x,2)/pow(std_landmark[0],2)) +
               (pow(obsy_map - landm_y,2)/pow(std_landmark[1],2)) )
            );
            if (DEBUG)
            {
                cout << "Found Landmark ID: "<<found_id<<endl;
                printf("obs x: %f; obs y: %f\n", obsx_map, obsy_map);
                printf("landm x: %f; landm y: %f\n", landm_x, landm_y);
                printf("%f ", prob);
            }

            weight_p *= prob;
        }
        if (DEBUG)
            printf("\n");
        particles[i].weight = weight_p;
        weights[i] = particles[i].weight;
    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    std::default_random_engine gen;
    std::discrete_distribution<int> d(weights.begin(), weights.end());

    std::vector<Particle> new_particles;

    int index=0;
    for(int n=0; n<num_particles; ++n) {
        index = d(gen);
        //copy particle, not call by reference
        Particle p = particles[index];
        new_particles.push_back(p);
    }
    this->particles = new_particles;
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
