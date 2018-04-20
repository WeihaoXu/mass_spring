#include "mass_spring.h"
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <iostream>
#include <glm/gtx/string_cast.hpp>

SpringNode::SpringNode(glm::vec3 init_pos, float init_mass, bool init_fixed)
				:position(init_pos), init_position(init_pos), mass(init_mass), fixed(init_fixed)
{
	force = glm::vec3(0.0, -mass * G, 0.0);
}

SpringNode::~SpringNode()
{
}

MassSpringSystem::MassSpringSystem(int init_x_size, int init_z_size)
									:x_size(init_x_size), z_size(init_z_size)
{
	// create nodes
	for(int x = 0; x < x_size; x++) {
		for(int z = 0; z < z_size; z++) {
			SpringNode curr_node(glm::vec3(x * grid_width_, 10.0, z * grid_width_), node_mass_, false);
			nodes_.push_back(curr_node);
		}
	}

	// connect nodes.
	for(int x = 0; x < x_size; x++) {
		for(int z = 0; z < z_size; z++) {
			SpringNode& curr_node = nodes_[getNodeIndex(x, z)];
			for(int delta_x = -1; delta_x <= 1; delta_x++) {
				for(int delta_z = -1; delta_z <= 1; delta_z++) {
					if(!isIndexValid(x + delta_x, z + delta_z) || delta_x == 0 && delta_z == 0) {
						continue;
					}
					curr_node.neighbors.push_back(&(nodes_[getNodeIndex(x + delta_x, z + delta_z)]));
					line_indices.push_back(glm::uvec2(getNodeIndex(x, z), getNodeIndex(x + delta_x, z + delta_z)));
				}
			}
		}
	}

	// set border fixed
	for(int x = 0; x < x_size; x++) {
		for(int z = 0; z < z_size; z++) {
			SpringNode& curr_node = nodes_[getNodeIndex(x, z)];
			// if(x == 0) {
			if(x == 0 || x == x_size - 1 || z == 0 || z == z_size - 1) {
				nodes_[getNodeIndex(x, z)].fixed = true;
			}
		}
	}


	std::cout << "system created, T = " << T_ << std::endl;
	

	srand (time(NULL));
	
	refreshCache();

}

bool MassSpringSystem::isIndexValid(int x, int z) {
	return x >= 0 && x < x_size && z >= 0 && z < z_size;
}

int MassSpringSystem::getNodeIndex(int x, int z) {
	return x * z_size + z;
}

MassSpringSystem::~MassSpringSystem() 
{

}

glm::vec3 MassSpringSystem::computeSingleForce(const SpringNode& curr_node, const SpringNode& nb_node) {
	float dist = glm::length(curr_node.position - nb_node.position);
	float init_dist = glm::length(curr_node.init_position - nb_node.init_position);
	glm::vec3 force = spring_k_ * (dist - init_dist) * glm::normalize(nb_node.position - curr_node.position);
	return force;
}



void MassSpringSystem::refreshCache() {
	node_positions.resize(nodes_.size());
	for(int i = 0; i < nodes_.size(); i++) {
		node_positions[i] = nodes_[i].position;
	}
}

const glm::vec3* MassSpringSystem::collectNodePositions() {
	return node_positions.data();
}

void MassSpringSystem::animate(float delta_t) {	// update system states and refresh cache.
	// update force
	for(int i = 0; i < nodes_.size(); i++) {
		SpringNode& curr_node = nodes_[i];
		if(curr_node.fixed) continue;	// anchor node, not movable
		curr_node.force = glm::vec3(0.0, -curr_node.mass * G, 0.0);
		for(int j = 0; j < curr_node.neighbors.size(); j++) {
			const SpringNode& nb_node = *(curr_node.neighbors[j]);
			curr_node.force += computeSingleForce(curr_node, nb_node);
		}
	}

	// update velocity and position (semi-Implicit Euler)
	for(auto& curr_node : nodes_) {
		if(curr_node.fixed) continue;
		glm::vec3 damper_force = - damper_ * curr_node.velocity;
		// std::cout << "damper force: " << glm::to_string(damper_force) 
		// 	<< ", spring force: " << glm::to_string(curr_node.force) << std::endl;
		curr_node.velocity += (curr_node.force + damper_force) / curr_node.mass * delta_t;
		curr_node.position += curr_node.velocity * delta_t;


		// curr_node.velocity += (curr_node.force / curr_node.mass) * delta_t;
		// curr_node.position += curr_node.velocity * delta_t;
		// curr_node.velocity -= damper_ * curr_node.velocity / curr_node.mass * delta_t;
		// curr_node.velocity = (1 - damper_ * delta_t) * curr_node.velocity;
	}

	refreshCache();
}

void MassSpringSystem::resetSystem() {	// update system states and refresh cache.
	// update force
	for(int i = 0; i < nodes_.size(); i++) {
		SpringNode& curr_node = nodes_[i];
		if(curr_node.fixed) continue;	// anchor node, not movable
		curr_node.position = curr_node.init_position;
		curr_node.velocity = glm::vec3(0.0);
		curr_node.force = glm::vec3(0.0, -curr_node.mass * G, 0.0);
	}
	
	refreshCache();
}

void MassSpringSystem::randomDisturb() {	// update system states and refresh cache.
	int x = std::floor((double)rand() / RAND_MAX * x_size);
	int z = std::floor((double)rand() / RAND_MAX * z_size);
	std::cout << "random change x: " << x << ", z: " << z << std::endl;
	SpringNode& curr_node = nodes_[getNodeIndex(x, z)];
	if(curr_node.fixed) {
		randomDisturb();
		return;
	}
	curr_node.position += (double)rand() / RAND_MAX * grid_width_ * 5.0;
	
	refreshCache();
}
