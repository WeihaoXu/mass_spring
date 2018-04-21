#include "mass_spring.h"
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <iostream>
#include <glm/gtx/string_cast.hpp>

SpringNode::SpringNode(int init_index, glm::vec3 init_pos, glm::vec3 curr_pos, float init_mass, bool init_fixed)
				:index(init_index), init_position(init_pos), position(curr_pos), mass(init_mass), fixed(init_fixed)
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
			glm::vec3 node_pos(x * grid_width_, 10.0, z * grid_width_);
			SpringNode* curr_node = new SpringNode(getNodeIndex(x, z), node_pos, node_pos, node_mass_, false);
			nodes_.push_back(curr_node);
		}
	}

	// connect nodes.
	for(int x = 0; x < x_size; x++) {
		for(int z = 0; z < z_size; z++) {
			SpringNode* curr_node = nodes_[getNodeIndex(x, z)];
			for(int delta_x = -1; delta_x <= 1; delta_x++) {
				for(int delta_z = -1; delta_z <= 1; delta_z++) {
					if(!isIndexValid(x + delta_x, z + delta_z) || delta_x == 0 && delta_z == 0) {
						continue;
					}
					curr_node->neighbors.push_back(nodes_[getNodeIndex(x + delta_x, z + delta_z)]);
					line_indices.push_back(glm::uvec2(getNodeIndex(x, z), getNodeIndex(x + delta_x, z + delta_z)));
				}
			}
		}
	}

	// set border fixed
	for(int x = 0; x < x_size; x++) {
		for(int z = 0; z < z_size; z++) {
			SpringNode* curr_node = nodes_[getNodeIndex(x, z)];
			// if(x == 0) {
			if(x == 0 || x == x_size - 1 || z == 0 || z == z_size - 1) {
				nodes_[getNodeIndex(x, z)]->fixed = true;
			}
			std::cout << "x = " << x << ", z = " << z << ", fixed? " << nodes_[getNodeIndex(x, z)]->fixed << std::endl;
		}
	}


	std::cout << "system created, T = " << T_ << std::endl;
	

	srand (time(NULL));
	
	refreshCache();

}

MassSpringSystem::~MassSpringSystem() 
{
	for(SpringNode* node : nodes_) {
		delete node;
	}

}

bool MassSpringSystem::isIndexValid(int x, int z) {
	return x >= 0 && x < x_size && z >= 0 && z < z_size;
}

int MassSpringSystem::getNodeIndex(int x, int z) {
	return x * z_size + z;
}



glm::vec3 MassSpringSystem::computeSingleForce(const SpringNode* curr_node, const SpringNode* nb_node) {
	float dist = glm::length(curr_node->position - nb_node->position);
	float init_dist = glm::length(curr_node->init_position - nb_node->init_position);
	glm::vec3 force = spring_k_ * (dist - init_dist) * glm::normalize(nb_node->position - curr_node->position);
	return force;
}



void MassSpringSystem::refreshCache() {
	std::cout << "to refresh cache. current nodes size: " << nodes_.size()  << std::endl;
	node_positions.resize(nodes_.size());
	line_indices.clear();

	// std::cout << "refresh cache" << std::endl;
	for(int i = 0; i < nodes_.size(); i++) {
		node_positions[i] = nodes_[i]->position;
		for(SpringNode* nb_node : nodes_[i]->neighbors) {
			// std::cout << "line index: (" << nodes_[i]->index << ", " << nb_node->index << ")" << std::endl;
			line_indices.push_back(glm::uvec2(nodes_[i]->index, nb_node->index));
		}
	}


}

const glm::vec3* MassSpringSystem::collectNodePositions() {
	return node_positions.data();
}

void MassSpringSystem::animate(float delta_t) {	// update system states and refresh cache.
	// update force
	for(int i = 0; i < nodes_.size(); i++) {
		SpringNode* curr_node = nodes_[i];
		if(curr_node->fixed) continue;	// anchor node, not movable
		curr_node->force = glm::vec3(0.0, -curr_node->mass * G, 0.0);
		for(int j = 0; j < curr_node->neighbors.size(); j++) {
			SpringNode* nb_node = curr_node->neighbors[j];
			curr_node->force += computeSingleForce(curr_node, nb_node);
		}
	}

	// https://stackoverflow.com/32776571/c-iterate-through-an-expanding-container/32776728
	// update velocity and position (semi-Implicit Euler)
	std::vector<SpringNode*> teared_new_nodes; 
	for(int i = 0; i < nodes_.size(); i++) {
		SpringNode* curr_node = nodes_[i];
		if(!curr_node->fixed) {
			glm::vec3 damper_force = - damper_ * curr_node->velocity;
			// std::cout << "damper force: " << glm::to_string(damper_force) 
			// 	<< ", spring force: " << glm::to_string(curr_node.force) << std::endl;
			curr_node->velocity += (curr_node->force + damper_force) / curr_node->mass * delta_t;
			curr_node->position += curr_node->velocity * delta_t;
		}
		
		// if(!curr_node->teared) {
		// 	for(int nb_idx = 0; nb_idx < curr_node->neighbors.size(); nb_idx++) {
		// 		SpringNode* nb_node = curr_node->neighbors[nb_idx];
		// 		checkTear(curr_node, nb_idx, 1.5, teared_new_nodes);
		// 	}
		// }
		
	}
	for(SpringNode* new_node : teared_new_nodes) {
		nodes_.push_back(new_node);
	}
	
	// std::cout << "before refresh cache" << std::endl;
	refreshCache();
	// std::cout << "done refresh cache" << std::endl;
}

void MassSpringSystem::checkTear(SpringNode* curr_node, int nb_idx, float max_deform_rate, 
									std::vector<SpringNode*>& teared_new_nodes) {
	SpringNode* nb_node = curr_node->neighbors[nb_idx];
	float curr_length = glm::length(curr_node->position - nb_node->position);
	float init_length = glm::length(curr_node->init_position - nb_node->init_position);
	if(std::abs(curr_length) > max_deform_rate * init_length) {
		glm::vec3 new_node_curr_pos = (curr_node->position + nb_node->position) / 2.0f;
		glm::vec3 new_node_init_pos = (curr_node->init_position + nb_node->init_position) / 2.0f;

		// SpringNode* new_node = new SpringNode(nodes_.size(), new_node_init_pos, new_node_curr_pos, curr_node->mass / 2.0, false);
		// new_node->teared = true;
		// std::cout << "push new node neighbor" << std::endl;
		
		// new_node->neighbors.push_back(curr_node);
		// std::cout << "done push new node neighbor" << std::endl;
	
		// teared_new_nodes.push_back(new_node);
		// std::cout << "done push new node" << std::endl;
	
		// curr_node->neighbors[nb_idx] = nodes_[nodes_.size() - 1];
		curr_node->neighbors.erase(curr_node->neighbors.begin() + nb_idx);
		// std::cout << "done update neighbor" << std::endl;

		// curr_node->neighbors.erase(curr_node->neighbors.begin() + nb_idx);
	}
}

void MassSpringSystem::resetSystem() {	// update system states and refresh cache.
	// update force
	for(int i = 0; i < nodes_.size(); i++) {
		SpringNode* curr_node = nodes_[i];
		if(curr_node->fixed) continue;	// anchor node, not movable
		curr_node->position = curr_node->init_position;
		curr_node->velocity = glm::vec3(0.0);
		curr_node->force = glm::vec3(0.0, -curr_node->mass * G, 0.0);
	}
	
	// refreshCache();
}

void MassSpringSystem::randomDisturb() {	// update system states and refresh cache.
	int x = std::floor((double)rand() / RAND_MAX * x_size);
	int z = std::floor((double)rand() / RAND_MAX * z_size);
	std::cout << "random change x: " << x << ", z: " << z << std::endl;
	SpringNode* curr_node = nodes_[getNodeIndex(x, z)];
	if(curr_node->fixed) {
		randomDisturb();
		return;
	}
	curr_node->position += (double)rand() / RAND_MAX * grid_width_ * 5.0;
	
	// refreshCache();
}
