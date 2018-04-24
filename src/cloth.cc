#include "cloth.h"
#include <glm/gtx/string_cast.hpp>


Particle::Particle(glm::vec3 init_position, float mass, glm::vec2 uv_coords):
			init_position_(init_position), position_(init_position), mass_(mass), uv_coords_(uv_coords)
{

}

Particle::~Particle() {

}

Spring::Spring(Particle* p1, Particle* p2, Spring* bend_spring):
			p1_(p1), p2_(p2), bend_spring_(bend_spring_)
{
}


Spring::~Spring() 
{
}


Cloth::Cloth(int x_size, int z_size):
		x_size_(x_size), z_size_(z_size)
{
	// build grid
	for(int x = 0; x < x_size_; x++) {
		float z_offset = (x % 2 == 0)? 0.0 : (0.5 * grid_width_);
		for(int z = 0; z < z_size_; z++) {
			glm::vec3 position(x * grid_width_, init_height_, z * grid_width_ + z_offset);
			glm::vec2 uv_coords(1.0 * x / (x_size_ - 1), 1.0 * z / (z_size_ - 1));
			// std::cout << "uv = " << glm::to_string(uv_coords) << std::endl;
			Particle* particle = new Particle(position, particle_mass_, uv_coords);
			particles_.push_back(particle);
			// std::cout << "particles " << glm::to_string(particle->position_) << std::endl;
		}
	}
	// std::cout << "cloth built, particle number: " << particles_.size() << std::endl;

	// create triangles
	for(int x = 0; x < x_size_; x++) {
		for(int z = 0; z < z_size_; z++) {
			if(x % 2 == 0) {

				if(gridCoordValid(x, z + 1) && gridCoordValid(x + 1, z)) {
					Triangle* triangle = new Triangle();
					triangle->particles_.push_back(particles_[getParticleIdx(x, z)]);
					triangle->particles_.push_back(particles_[getParticleIdx(x, z + 1)]);
					triangle->particles_.push_back(particles_[getParticleIdx(x + 1, z)]);
					triangles_.push_back(triangle);
				}
				if(gridCoordValid(x, z + 1) && gridCoordValid(x - 1, z)) {
					Triangle* triangle = new Triangle();
					triangle->particles_.push_back(particles_[getParticleIdx(x, z)]);
					triangle->particles_.push_back(particles_[getParticleIdx(x - 1, z)]);
					triangle->particles_.push_back(particles_[getParticleIdx(x, z + 1)]);
					triangles_.push_back(triangle);
				}
			}
			else {
				if(gridCoordValid(x - 1, z + 1) && gridCoordValid(x, z + 1)) {
					Triangle* triangle = new Triangle();
					triangle->particles_.push_back(particles_[getParticleIdx(x, z)]);
					triangle->particles_.push_back(particles_[getParticleIdx(x - 1, z + 1)]);
					triangle->particles_.push_back(particles_[getParticleIdx(x, z + 1)]);
					triangles_.push_back(triangle);
				}

				if(gridCoordValid(x + 1, z + 1) && gridCoordValid(x, z + 1)) {
					Triangle* triangle = new Triangle();
					triangle->particles_.push_back(particles_[getParticleIdx(x, z)]);
					triangle->particles_.push_back(particles_[getParticleIdx(x, z + 1)]);
					triangle->particles_.push_back(particles_[getParticleIdx(x + 1, z + 1)]);
					triangles_.push_back(triangle);
				}
			}
		}
	}

	// structure springs

	// for(int x = 0; x < x_size_; x++) {
	// 	if(x % 2 == 0) {
	// 		for(int z = 0; z < z_size_; z++) {
	// 			if(gridCoordValid(x + 1, z - 1)) {
	// 				Spring* s = new Spring();
	// 				springs_.push_back(s);
	// 			}
	// 			if(gridCoordValid()) {
	// 				Spring* s = new Spring();
	// 				springs_.push_back(s);
	// 			}
	// 			if(gridCoordValid()) {
	// 				Spring* s = new Spring();
	// 				springs_.push_back(s);
	// 			}
	// 		}
	// 	}
	// 	else {
	// 		for(int z = 0; z < z_size_; z++) {
	// 			if(gridCoordValid(x + 1, z - 1)) {
	// 				Spring* s = new Spring();
	// 				springs_.push_back(s);
	// 			}
	// 			if(gridCoordValid()) {
	// 				Spring* s = new Spring();
	// 				springs_.push_back(s);
	// 			}
	// 			if(gridCoordValid()) {
	// 				Spring* s = new Spring();
	// 				springs_.push_back(s);
	// 			}
	// 		}
	// 	}
	// }

	refreshCache();

}


Cloth::~Cloth() {

}

void Cloth::refreshCache() {
	vertices.clear();
	cloth_uv_coords.clear();
	for(Triangle* triangle : triangles_) {
		for(Particle* p : triangle->particles_) {
			// std::cout << "particle pushed to cache: " << glm::to_string(p->position_) << std::endl;
			vertices.push_back(p->position_);
			cloth_uv_coords.push_back(p->uv_coords_);
		}
	}


}

void Cloth::animate(float delta_t) {
	refreshCache();
}


int Cloth::getParticleIdx(int x, int z) {
	return x * z_size_ + z;
}
bool Cloth::gridCoordValid(int x, int z) {
	return x >= 0 && x < x_size_ && z >= 0 && z < z_size_;
}






















