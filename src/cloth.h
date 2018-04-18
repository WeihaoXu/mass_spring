#ifndef CLOTH_H
#define CLOTH_H

#include "gui.h"

#include <vector>

#define G 9.8f

struct SpringNode {
public:
	SpringNode(glm::vec3 init_pos, float mass, bool init_fixed = false);
	~SpringNode();

	glm::vec3 position;
	glm::vec3 init_position;
	glm::vec3 velocity = glm::vec3(0.0f, 0.0f, 0.0f);
	glm::vec3 force;

	std::vector<SpringNode*> neighbors;
	std::vector<float> init_nb_dists;

	float mass;
	bool fixed;
	int index;


};


// struct AnchorNode {
// public:
// 	AnchorNode();
// 	~AnchorNode();
// 	SpringNode* node;
// 	glm::vec3 position;	// anchor position
// };

class MassSpringSystem {
public:
	MassSpringSystem(int width, int height);
	~MassSpringSystem();

	glm::vec3 computeSingleForce(const SpringNode& curr_node, const SpringNode& nb_node, int init_dist);

	void setNodeFixed(int idx);
	void setNodeMovable(int idx);

	void refreshCache();	// copy node positions to opengl buffer
	const glm::vec3* collectNodePositions();

	void animate(float delta_t);

	std::vector<glm::vec3> node_positions;


private:
	const float spring_k_ = 1.0;
	const float energy_loss_ = 0.9;
	std::vector<SpringNode> nodes_;



};














#endif