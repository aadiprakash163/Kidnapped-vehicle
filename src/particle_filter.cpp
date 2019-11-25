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
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	
	num_particles = 100; // number of particles
	
	// noise distribution definition
	normal_distribution<double> uncertain_x(x , std[0]);
	normal_distribution<double> uncertain_y(y , std[1]);
	normal_distribution<double> uncertain_theta(theta , std[2]);

	// random number generator
	default_random_engine gen;

	cout<<"Initialization..............................."<<endl;
	cout<<"gps readings: "<<endl<<" x, y, theta: "<<x<<", "<<y<<", "<<theta<<endl; 

	// initialize particels
	for(int i = 0; i < num_particles; i++){
		Particle p;
		p.id = i;
		p.x = uncertain_x(gen);
		p.y = uncertain_y(gen);		
		p.theta = uncertain_theta(gen);
		cout<<"x: "<<p.x<<"\ty: "<<p.y<<"\ttheta: "<<p.theta<<endl;
		p.weight = 1;

		particles.push_back(p);
	}

	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	normal_distribution<double> std_pos_x(0 , std_pos[0]);
	normal_distribution<double> std_pos_y(0 , std_pos[1]);
	normal_distribution<double> std_pos_theta(0 , std_pos[2]);

	default_random_engine gen;
	
	// perform a movement step for all particles	
	for(int i = 0; i < num_particles; i++){

		if(fabs(yaw_rate) < 0.001f){
			particles[i].x += velocity * delta_t * cos(particles[i].theta);
			particles[i].y += velocity * delta_t * sin(particles[i].theta);
		} 
		else{
			particles[i].x += (velocity/yaw_rate) * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta)); 
			particles[i].y += (velocity/yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
			particles[i].theta += yaw_rate * delta_t; 
		}

		// add gaussian noise to all the particles
		particles[i].x += std_pos_x(gen);
		particles[i].y += std_pos_y(gen);
		particles[i].theta += std_pos_theta(gen);
		
	}


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> lm_in_range, std::vector<LandmarkObs>& obs_map_positions) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
	// make a vector to store closest predictions to observations
	// std::vector<int> closest_prediction;

	// iterate through all the measurements
	for(int i = 0; i < obs_map_positions.size(); i++){
		int temp_index = -1;
		double min_dist = numeric_limits<double>::max();

		for(int j = 0; j < lm_in_range.size(); j++){ //iterate through all the predictions			
			double distance = dist(obs_map_positions[i].x,obs_map_positions[i].y,lm_in_range[j].x,lm_in_range[j].y);
			if ( distance < min_dist){				
				min_dist = distance;
				temp_index = j;
			}
		}		
		obs_map_positions[i].id = lm_in_range[temp_index].id;		




	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// Argument details are as following:
	// 1). sensor_range: sensing range of sensor
	// 2). std_landmarks[]: standard deviations of landmark positions
	// 3). observations: noisy measurements of landmarks (vehicle coordinates)
	// 4). map_landmarks: landmarks in real world (Map coordinates)

	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution	
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	// get all the landmarks in sensor_range from particles

	// cout<<"\nObservation size: "<<observations.size()<<"\t Map landmark size: "<<map_landmarks.landmark_list.size()<<endl;

	// For each particle do following
	for(int i = 0; i < num_particles ; i++){ // Get pose
		double p_x = particles[i].x;
		double p_y = particles[i].y;
		double p_theta = particles[i].theta;

		std::vector<LandmarkObs> lm_in_range;

		// Get list of landmarks in the sensor range from the particle
		for(int j = 0; j < map_landmarks.landmark_list.size(); j++){
			double l_x = map_landmarks.landmark_list[j].x_f;
			double l_y = map_landmarks.landmark_list[j].y_f;
			int id = map_landmarks.landmark_list[j].id_i;
			if (dist(p_x,p_y,l_x,l_y) <= sensor_range){
				lm_in_range.push_back(LandmarkObs{id, l_x, l_y});				
			}
		}
		
		std::vector<LandmarkObs> obs_map_positions; 
		for(int j = 0; j < observations.size(); j++){
			double t_x = p_x + observations[j].x*cos(p_theta) - observations[j].y*sin(p_theta);
			double t_y = p_y + observations[j].x*sin(p_theta) + observations[j].y*cos(p_theta);
			int id = observations[j].id;			
			obs_map_positions.push_back(LandmarkObs{ id, t_x, t_y });
		}

		

		// Associate observations with nearest landmarks
		dataAssociation(lm_in_range, obs_map_positions);


		// Reinitialize weight
		particles[i].weight = 1.0;
		

		// Calculate weight of the particle using multi-variate gaussian
		for(int j = 0; j < obs_map_positions.size(); j++){
			double map_x = obs_map_positions[j].x;
			double map_y = obs_map_positions[j].y;

			int lm_id = obs_map_positions[j].id;
			double true_lm_x , true_lm_y;
			for(int k = 0; k < lm_in_range.size(); k++){			
				if (lm_in_range[k].id == lm_id){
					true_lm_x = lm_in_range[k].x;
					true_lm_y = lm_in_range[k].y;					
					break;
					
				}
			}	


			double std_x = std_landmark[0];
			double std_y = std_landmark[1];

			double weight = ( 1/(2*M_PI*std_x*std_y)) * exp( -( pow(true_lm_x-map_x,2)/(2*pow(std_x, 2)) + (pow(true_lm_y-map_y,2)/(2*pow(std_y, 2))) ) );
			
			particles[i].weight *= weight;
			
		}
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  	// get current weights
  	vector<double> current_weights;

  	// random number generator
  	default_random_engine gen;
  	for (int i = 0; i < num_particles; i++) {
    	current_weights.push_back(particles[i].weight);
  	}
  	
  	// generate random integer for index
  	uniform_int_distribution<int> rand_int(0, num_particles-1);
  	auto index = rand_int(gen);

  	
  	double max_weight = *max_element(current_weights.begin(), current_weights.end());

  	// generate random real number for beta
  	uniform_real_distribution<double> rand_real(0.0, max_weight);

  	double beta = 0.0;
  	vector<Particle> updated_particles;
  	// spin the resample wheel!
  	for (int i = 0; i < num_particles; i++) {
    	beta += rand_real(gen) * 2.0;
    	while (beta > current_weights[index]) {
      	beta -= current_weights[index];
      	index = (index + 1) % num_particles;
    	}
    
    	updated_particles.push_back(particles[index]);
  	
  	}


  	// Calculate sum of all weights for normalization
  	// double sum_of_weights = 0.0;
  	// for(int j = 0; j < num_particles; j++){
  	// 	sum_of_weights += updated_particles[j].weight;
  	// }


  	// // Normalize weights
  	// for(int j = 0; j < num_particles; j++){
  	// 	updated_particles[j].weight /= sum_of_weights ;
  	// }


  	particles = updated_particles;
}




Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

