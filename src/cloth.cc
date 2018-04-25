#include "cloth.h"
#include <glm/gtx/string_cast.hpp>
#include <map>


Particle::Particle(glm::vec3 init_position, float mass, glm::vec2 uv_coords, int grid_x, int grid_z):
			init_position_(init_position), position_(init_position), mass_(mass), uv_coords_(uv_coords),
			grid_x_(grid_x), grid_z_(grid_z)
{
	resetForce();
	setMovable();
	// setFixed();
}

Particle::~Particle() {

}

void Particle::resetForce() {
	force_ = glm::vec3(0.0f, - 1.0 * mass_ * G, 0.0f);
}


void Particle::addForce(glm::vec3 f) {
	force_ += f;
}

void Particle::setFixed() {
	fixed_ = true;
}

void Particle::setMovable() {
	fixed_ = false;
}

Spring::Spring(Particle* p1, Particle* p2, float k):
			p1_(p1), p2_(p2), k_(k)
{
	init_length_ = glm::length(p1_->position_ - p2_->position_);
}


Spring::~Spring() 
{
}

void Spring::computeForceQuantity() {
	float curr_length = glm::length(p1_->position_ - p2_->position_);
	force_quantity_ = (init_length_ - curr_length) / init_length_ * k_;

}

void Spring::applyForce() {
	glm::vec3 force1 = glm::normalize(p1_->position_ - p2_->position_) * force_quantity_;
	glm::vec3 force2 = -1.0f * force1;
	p1_->addForce(force1);
	p2_->addForce(force2);
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
			Particle* particle = new Particle(position, particle_mass_, uv_coords, x, z);
			particles_.push_back(particle);
			// std::cout << "particles " << glm::to_string(particle->position_) << std::endl;
		}
	}
	// std::cout << "cloth built, particle number: " << particles_.size() << std::endl;

	// set two anchor nodes. For experiments.
	particles_[getParticleIdx(0, 0)]->setFixed();
	particles_[getParticleIdx(0, z_size_ - 1)]->setFixed();

	// create triangles
	for(int x = 0; x < x_size_; x++) {
		for(int z = 0; z < z_size_; z++) {
			if(x % 2 == 0) {

				if(gridCoordValid(x, z + 1) && gridCoordValid(x + 1, z)) {
					Triangle* triangle = new Triangle();
					triangle->particles_.push_back(particles_[getParticleIdx(x, z)]);
					triangle->particles_.push_back(particles_[getParticleIdx(x, z + 1)]);
					triangle->particles_.push_back(particles_[getParticleIdx(x + 1, z)]);
					triangles_.insert(triangle);
				}
				if(gridCoordValid(x, z + 1) && gridCoordValid(x - 1, z)) {
					Triangle* triangle = new Triangle();
					triangle->particles_.push_back(particles_[getParticleIdx(x, z)]);
					triangle->particles_.push_back(particles_[getParticleIdx(x - 1, z)]);
					triangle->particles_.push_back(particles_[getParticleIdx(x, z + 1)]);
					triangles_.insert(triangle);
				}
			}
			else {
				if(gridCoordValid(x - 1, z + 1) && gridCoordValid(x, z + 1)) {
					Triangle* triangle = new Triangle();
					triangle->particles_.push_back(particles_[getParticleIdx(x, z)]);
					triangle->particles_.push_back(particles_[getParticleIdx(x - 1, z + 1)]);
					triangle->particles_.push_back(particles_[getParticleIdx(x, z + 1)]);
					triangles_.insert(triangle);
				}

				if(gridCoordValid(x + 1, z + 1) && gridCoordValid(x, z + 1)) {
					Triangle* triangle = new Triangle();
					triangle->particles_.push_back(particles_[getParticleIdx(x, z)]);
					triangle->particles_.push_back(particles_[getParticleIdx(x, z + 1)]);
					triangle->particles_.push_back(particles_[getParticleIdx(x + 1, z + 1)]);
					triangles_.insert(triangle);
				}
			}
		}
	}

	// create structural springs
	std::map<Particle*, std::map<Particle*, Spring*>> springMap;
	for(Triangle* triangle : triangles_) {
		for(int idx = 0; idx < 3; idx++) {
			Particle* p1 = triangle->particles_[idx];
			Particle* p2 = triangle->particles_[(idx + 1) % 3];
			if(springMap[p1][p2] == nullptr && springMap[p2][p1] == nullptr) {
				Spring* s = new Spring(p1, p2, struct_k_);	// problem: how find bending spring?
				s->triangles_.push_back(triangle);
				p1->springs_.push_back(s);
				p2->springs_.push_back(s);
				springs_.insert(s);
				springMap[p1][p2] = s;
				springMap[p2][p1] = s;

			}
			else {
				springMap[p1][p2]->triangles_.push_back(triangle);
			}
		}
	}


	// add bending springs
	for(Spring* spring : springs_) {
		Particle* p1;
		Particle* p2;
		if(spring->p1_->grid_x_ < spring->p2_->grid_x_) {
			p1 = spring->p1_;
			p2 = spring->p2_;
		}
		else {
			p1 = spring->p2_;
			p2 = spring->p1_;
		}

		// create bending springs
		int bend_x1 = -1, bend_z1 = -1, bend_x2 = -1, bend_z2 = -1;

		if(p1->grid_x_ % 2 == 0) {
			if(p1->grid_x_ == p2->grid_x_) {
				bend_x1 = p1->grid_x_ - 1; 
				bend_z1 = std::min(p1->grid_z_, p2->grid_z_);
				bend_x2 = p1->grid_x_ + 1; 
				bend_z2 = std::min(p1->grid_z_, p2->grid_z_);
			}
			else if(p1->grid_z_ == p2->grid_z_) {
				bend_x1 = p1->grid_x_;
				bend_z1 = p1->grid_z_ + 1;
				bend_x2 = p2->grid_x_;
				bend_z2 = p2->grid_z_ - 1;
				
			}
			else {
				bend_x1 = p1->grid_x_;
				bend_z1 = p1->grid_z_ - 1;
				bend_x2 = p2->grid_x_;
				bend_z2 = p2->grid_z_ + 1;
				
			}
		}
		else {
			if(p1->grid_x_ == p2->grid_x_) {
				bend_x1 = p1->grid_x_ - 1;
				bend_z1 = std::max(p1->grid_z_, p2->grid_z_);
				bend_x2 = p1->grid_x_ + 1;
				bend_z2 = std::max(p1->grid_z_, p2->grid_z_);
				
			}
			else if(p1->grid_z_ == p2->grid_z_) {
				bend_x1 = p1->grid_x_;
				bend_z1 = p1->grid_z_ - 1;
				bend_x2 = p2->grid_x_;
				bend_z2 = p2->grid_z_ + 1;
			}
			else {
				bend_x1 = p1->grid_x_;
				bend_z1 = p1->grid_z_ + 1;
				bend_x2 = p2->grid_x_;
				bend_z2 = p2->grid_z_ - 1;
			}

		}
		if(gridCoordValid(bend_x1, bend_z1) && gridCoordValid(bend_x2, bend_z2)) {
			Spring* bend_spring = new Spring(particles_[getParticleIdx(bend_x1, bend_z1)], 
										particles_[getParticleIdx(bend_x2, bend_z2)],
										bend_sheer_k_);
			spring->bend_spring_ = bend_spring;
		}	
	}


	// update cache vertices
	refreshCache();

}


Cloth::~Cloth() {

}

void Cloth::refreshCache() {
	// vertices and uv_coords
	vertices.clear();
	cloth_uv_coords.clear();
	for(Triangle* triangle : triangles_) {
		for(Particle* p : triangle->particles_) {
			// std::cout << "particle pushed to cache: " << glm::to_string(p->position_) << std::endl;
			vertices.push_back(p->position_);
			cloth_uv_coords.push_back(p->uv_coords_);
		}
	}


	// spring linemesh
	struct_spring_vertices.clear();
	bend_spring_vertices.clear();
	for(Spring* s : springs_) {
		struct_spring_vertices.push_back(s->p1_->position_);
		struct_spring_vertices.push_back(s->p2_->position_);

		if(s->bend_spring_) {
			// std::cout << "push bend spring" << std::endl;
		
			bend_spring_vertices.push_back(s->bend_spring_->p1_->position_);
			bend_spring_vertices.push_back(s->bend_spring_->p2_->position_);

		}
	}
	// std::cout << "end push bend spring" << std::endl;




}

void Cloth::animate(float delta_t) {
	// clear all forces except for gravity
	for(Particle* particle : particles_) {
		particle->resetForce();
	}

	// update forces
	for(Spring* struct_s : springs_) {
		struct_s->computeForceQuantity();
		// TODO: if force quantity exceeds limit, break the spring
		struct_s->applyForce();
		
		if(struct_s->bend_spring_) {
			struct_s->bend_spring_->computeForceQuantity();
			struct_s->bend_spring_->applyForce();

		}
	}

	// update particle velocity and positions
	for(Particle* particle : particles_) {
		// Update velocity and positions
		if(!particle->fixed_) {
			glm::vec3 damper_force = -damper_ * particle->velocity_;
			glm::vec3 acceleration = (particle->force_ + damper_force) / particle->mass_;
			particle->velocity_ += acceleration * delta_t;
			particle->position_ += particle->velocity_ * delta_t;

		}
	}


	refreshCache();
}


int Cloth::getParticleIdx(int x, int z) {
	return x * z_size_ + z;
}
bool Cloth::gridCoordValid(int x, int z) {
	return x >= 0 && x < x_size_ && z >= 0 && z < z_size_;
}





















