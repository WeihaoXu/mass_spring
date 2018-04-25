#ifndef CLOTH_H
#define CLOTH_H

#include <vector>
#include <glm/glm.hpp>
#include <iostream>
#include <unordered_set>

#define G 9.8
#define PI 3.1416


struct Spring;
struct Triangle;


struct Particle {
	Particle(glm::vec3 init_position, float mass, glm::vec2 uv_coords, int grid_x, int grid_z);
	~Particle();
	
	void resetForce();
	void addForce(glm::vec3 f);

	void setFixed();
	void setMovable();

	glm::vec3 position_;
	glm::vec3 init_position_;
	glm::vec3 force_;
	glm::vec3 velocity_;

	glm::vec2 uv_coords_;
	glm::vec2 grid_coords;

	int grid_x_, grid_z_;

	std::vector<Spring*> springs_;
	float mass_;
	bool fixed_;

};

struct Triangle {
	// std::vector<Spring*> springs;
	std::vector<Particle*> particles_;

};

struct Spring {
	Spring(Particle* p1, Particle* p2, float k);
	~Spring();

	void computeForceQuantity();
	void applyForce();

	std::vector<Particle*> particles_;	// two particles
	std::vector<Triangle*> triangles_;	// two triangles

	Particle* p1_;
	Particle* p2_;
	
	float force_quantity_;

	Spring* bend_spring_ = nullptr;
	float init_length_;
	float k_;
};

class Cloth {

public:
	Cloth(int x_size, int z_size);
	~Cloth();
	void animate(float delta_t);




	std::vector<glm::vec3> vertices;
	std::vector<glm::vec2> cloth_uv_coords;
	// std::vector<glm::uvec3> faces;

	std::vector<glm::vec3> spring_vertices;	// used to linemesh springs. For debug use. 
	std::vector<glm::vec3> bend_spring_vertices;

private:
	int getParticleIdx(int x, int z);
	bool gridCoordValid(int x, int z);
	void refreshCache();



	std::vector<Particle*> particles_;
	std::unordered_set<Triangle*> triangles_;
	std::unordered_set<Spring*> springs_;
	// std::vector<Triangle*> triangles_;
	// std::vector<Spring*> springs_;

	int x_size_, z_size_;

	const float grid_width_ = 10.0;

	const float struct_k_ = 1.0;
	const float bend_k_ = 1.0;
	
	const float damper_ = 0.1;


	const float bend_sheer_k = 1.0;
	const float particle_mass_ = 0.1;
	const float init_height_ = 0.0;






};


#endif	// end define CLOTH_H