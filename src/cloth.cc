#include "cloth.h"

SpringNode::SpringNode(glm::vec3 init_pos, float init_mass, bool init_fixed)
				:position(init_pos), mass(init_mass), fixed(init_fixed)
{
	force = glm::vec3(0.0, -mass * G, 0.0);
}

SpringNode::~SpringNode()
{
}

MassSpringSystem::MassSpringSystem(int width, int height) {

}

MassSpringSystem::~MassSpringSystem() 
{

}

glm::vec3 MassSpringSystem::computeSingleForce(const SpringNode& curr_node, const SpringNode& nb_node, int init_dist) {
	float dist = glm::length(curr_node.position - nb_node.position);
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
			curr_node.force += computeSingleForce(curr_node, nb_node, curr_node.init_nb_dists[j]);
		}
	}

	// update velocity and position
	for(auto& curr_node : nodes_) {
		if(curr_node.fixed) continue;
		curr_node.velocity += (curr_node.force / curr_node.mass) * delta_t;
		curr_node.position += curr_node.velocity * delta_t;
		curr_node.velocity *= energy_loss_;
	}

	refreshCache();
}
