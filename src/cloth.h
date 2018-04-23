#ifndef CLOTH_H
#define CLOTH_H

#include <vector>
#include <glm/glm.hpp>

#define G 9.8
#define PI 3.1416

using namespace std;


struct Particle {
	
	glm::vec3 position;
	glm::vec3 force;
	glm::vec3 velocity;

	glm::vec2 uv_coords;
	vector<Spring*> springs;
	float mass;

};

struct Triangle {
	// vector<Spring*> springs;
	vector<Particle*> particles;

};

struct Spring {

	Particle* p1;
	Particle* p2;

	Particle* bend_p1;
	Particle* bend_p2;

	float init_length;
};

class Cloth {

public:
	vector<glm::vec3> vertices;
	vector<glm::uvec3> faces;


private:
	vector<Particle*> particles_;
	vector<Triangle*> triangles_;
	vector<Spring*> springs_;

	const float damper_;
	const float k_;

	const float partile_mass_;




};


#endif