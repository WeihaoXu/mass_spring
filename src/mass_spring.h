#ifndef MASS_SPRING_H
#define MASS_SPRING_H

#include <glm/glm.hpp>
#include <vector>

#define G 9.8f
#define PI 3.1416f

struct SpringNode {
public:
	SpringNode(int index, glm::vec3 init_pos, float mass, bool init_fixed = false);
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
	MassSpringSystem(int x_size, int z_size);
	~MassSpringSystem();

	glm::vec3 computeSingleForce(const SpringNode& curr_node, const SpringNode& nb_node);

	void setNodeFixed(int idx);
	void setNodeMovable(int idx);

	void refreshCache();	// copy node positions to opengl buffer
	const glm::vec3* collectNodePositions();

	void animate(float delta_t);
	void checkTear(SpringNode& curr_node, SpringNode& nb_node, float max_deform_rate, int nb_idx);

	void resetSystem();
	void randomDisturb();
	float getPeriod() {return T_;}

	std::vector<glm::vec3> node_positions;
	std::vector<glm::uvec2> line_indices;	// indices for line mesh
	std::vector<SpringNode> nodes_;


private:
	bool isIndexValid(int x, int z);
	int getNodeIndex(int x, int z);
	const float node_mass_ = 10.0;

	const float spring_k_ = 100.0;

	float T_ = 2 * PI * std::sqrt(node_mass_ / spring_k_);

	// const float energy_loss_ = 0.95;
	const float damper_ = (0.05) * (2 * std::sqrt(node_mass_ * spring_k_));
	const float grid_width_ = 1.0;
	int x_size, z_size;
	



};














#endif	// end define MASS_SPRING_H