#ifndef CLOTH_H
#define CLOTH_H

#include <vector>
#include <glm/glm.hpp>
#include <iostream>
#include <unordered_set>
#include <glm/gtx/string_cast.hpp>
#include <map>


#define G 9.8
#define PI 3.1416


struct Spring;
struct Triangle;

// Nodes in the spring.
struct Particle {
	Particle(glm::vec3 init_position, glm::vec3 curr_position, float mass, glm::vec2 uv_coords, int grid_x = -1, int grid_z = -1);
	~Particle();
	
	void resetForce();	// clear all forces except for gravity.
	void addForce(glm::vec3 f);
	void setFixed();
	void setMovable();

	glm::vec3 position_;
	glm::vec3 init_position_;
	glm::vec3 force_;
	glm::vec3 velocity_;

	glm::vec2 uv_coords_;
	glm::vec2 grid_coords;

	std::vector<Spring*> springs_;

	int grid_x_, grid_z_;
	float mass_;
	bool fixed_;


};

struct Triangle {
	Triangle(Particle* p1, Particle* p2, Particle* p3);
	Triangle();
	~Triangle();
	std::vector<Particle*> particles_;	// length == 3. Three particles.
};

struct Spring {

	Spring(Particle* p1, Particle* p2, float k);	// k is the spring constant
	~Spring();

	void computeForceQuantity();	// compute the force quantity, and store it in force_quantity_
	void applyForce();	// compute two force vectors, and apply them to two particles connected to the spring.
	void replaceTriangle(Triangle* t_old, Triangle* t_new);

	std::vector<Particle*> nb_particles_;	// two particles neighboring to but not owned by the spring.
	std::vector<Triangle*> triangles_;	// a spring is neighbor to either 1 or 2 triangles.

	Particle* p1_;
	Particle* p2_;
	Spring* bend_spring_ = nullptr;	// A bending spring (if there is one) related to this structural spring.
									// If this spring itself is a bending spring, this attribute will simply be nullptr.
	
	float force_quantity_;
	float init_length_;
	float k_;
	float max_deform_rate_ = 0.1f;

};

class Cloth {

public:
	Cloth(int x_size, int z_size);
	~Cloth();
	void animate(float delta_t);	// recalculate the forces, velocities and positions. Finally update cache
	void randomTear();

	// The following vectors are cache for GPU rendering.
	std::vector<glm::vec3> vertices;		// for rendering the cloth
	std::vector<glm::vec2> cloth_uv_coords;	// for texture mapping the the future.
	std::vector<glm::vec3> struct_spring_vertices;	// used to linemesh springs. For debug use. 
	std::vector<glm::vec3> bend_spring_vertices;	// used to linemesh springs. For debug use. 

private:
	int getParticleIdx(int x, int z);
	bool gridCoordValid(int x, int z);	
	void refreshCache();	// update the cache for rendiering
	void tear(Spring* s);
	Particle* getNeighborParticle(Triangle* t1, Spring* s);
	bool containsStructSpring(Particle* p1, Particle* p2);
	Spring* addStructSpring(Particle* p1, Particle* p2, float k);
	Spring* getStructSpring(Particle* p1, Particle* p2);
	void removeStructSpring(Spring* s);


	std::vector<Particle*> particles_;
	std::unordered_set<Triangle*> triangles_;	//stored in a hashset for constant time access, modify and delete
	std::unordered_set<Spring*> springs_;		//stored in a hashset for constant time access, modify and delete

	std::map<Particle*, std::map<Particle*, Spring*>> spring_map_; // key: particle pairs. Value: springs.


	int x_size_, z_size_;
	const float grid_width_ = 2.0;
	const float struct_k_ = 100.0;	// spring constant of bending springs
	const float bend_sheer_k_ = 10.0;	// spring constant of bending springs. (there bending springs also used as sheering springs)
	const float damper_ = 0.06;
	const float particle_mass_ = 0.1;	// init mass of every particle.
	const float init_height_ = 0.0;		// init height of the cloth. (i.e. init z position of all particles)

};



#endif	// end define CLOTH_H

