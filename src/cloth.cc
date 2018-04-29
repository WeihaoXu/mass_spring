#include "cloth.h"



Particle::Particle(glm::vec3 init_position, glm::vec3 curr_position, float mass, glm::vec2 uv_coords, int grid_x, int grid_z):
			init_position_(init_position), position_(curr_position), mass_(mass), uv_coords_(uv_coords),
			grid_x_(grid_x), grid_z_(grid_z), is_secondary_(false)
{
	resetForce();
	setMovable();
	// setFixed();
}

Particle::Particle(glm::vec3 init_position, glm::vec3 curr_position, float mass, glm::vec2 uv_coords, bool is_secondary):
			init_position_(init_position), position_(curr_position), mass_(mass), uv_coords_(uv_coords),
			grid_x_(-1), grid_z_(-1), is_secondary_(is_secondary)
{
	resetForce();
	setMovable();
}

Particle::Particle (const Particle &old_obj):
			init_position_(old_obj.init_position_), position_(old_obj.position_), force_(old_obj.force_), velocity_(old_obj.velocity_), 
			uv_coords_(old_obj.uv_coords_), grid_x_(old_obj.grid_x_), grid_z_(old_obj.grid_z_), mass_(old_obj.mass_), 
			fixed_(old_obj.fixed_), is_secondary_(old_obj.is_secondary_) 

{

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

Triangle::Triangle(Particle* p1, Particle* p2, Particle* p3) {
	particles_.push_back(p1);
	particles_.push_back(p2);
	particles_.push_back(p3);
}

Triangle::Triangle() {

}

Triangle::~Triangle() {

}

Spring::Spring(Particle* p1, Particle* p2, float k, bool is_secondary):
			p1_(p1), p2_(p2), k_(k), is_secondary_(is_secondary)
{
	init_length_ = glm::length(p1_->position_ - p2_->position_);
}


Spring::~Spring() 
{
	if(bend_spring_) {
		delete bend_spring_;
	}
}

void Spring::computeForceQuantity() {
	float curr_length = glm::length(p1_->position_ - p2_->position_);
	float deform_rate = (init_length_ - curr_length) / init_length_;
	if(fabs(deform_rate) > max_deform_rate_) {	// constrains. Anti-superelastic.
		deform_rate = deform_rate * (fabs(deform_rate) / max_deform_rate_);
	}
	force_quantity_ = deform_rate * k_;

}

void Spring::applyForce() {
	glm::vec3 force1 = glm::normalize(p1_->position_ - p2_->position_) * force_quantity_;
	glm::vec3 force2 = -1.0f * force1;
	p1_->addForce(force1);
	p2_->addForce(force2);
}

void Spring::replaceTriangle(Triangle* t_old, Triangle* t_new) {
	for(int i = 0; i < triangles_.size(); i++) {
		if(triangles_[i] == t_old) {
			triangles_[i] = t_new;
			return;
		}
	}
}

void Spring::replaceParticle(Particle* p_old, Particle* p_new) {
	if(p1_ == p_old) {
		p1_ = p_new;
	}
	else if(p2_ == p_old) {
		p2_ = p_new;
	}
	else {
		std::cout << "the particle you want to replace in the spring doesn't exist" << std::endl;
		throw("the particle you want to replace in the spring doesn't exist");
	}
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
			Particle* particle = new Particle(position, position, particle_mass_, uv_coords, x, z);
			particles_.push_back(particle);
			// std::cout << "particles " << glm::to_string(particle->position_) << std::endl;
		}
	}
	// std::cout << "cloth built, particle number: " << particles_.size() << std::endl;

	// set two anchor nodes. For experiments.
	particles_[getParticleIdx(0, 0)]->setFixed();
	particles_[getParticleIdx(0, z_size_ - 1)]->setFixed();
	particles_[getParticleIdx(x_size_ - 1, 0)]->setFixed();
	particles_[getParticleIdx(x_size_ - 1, z_size_ - 1)]->setFixed();

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
	std::map<Particle*, std::map<Particle*, Spring*>> spring_map_;
	for(Triangle* triangle : triangles_) {
		for(int idx = 0; idx < 3; idx++) {
			Particle* p1 = triangle->particles_[idx];
			Particle* p2 = triangle->particles_[(idx + 1) % 3];
			if(!containsStructSpring(p1, p2)) {
				Spring* s = addStructSpring(p1, p2, struct_k_, false);	// problem: how find bending spring?
				s->triangles_.push_back(triangle);
				// std::cout << "add triangle when create spring" << std::endl;
			}
			else {
				getStructSpring(p1, p2) ->triangles_.push_back(triangle);
				// std::cout << "add triangle to existing spring" << std::endl;
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

		// push neighboring particles to this spring. Will be used in tearing. 
		// if(gridCoordValid(bend_x1, bend_z1)) {
		// 	spring->nb_particles_.push_back(particles_[getParticleIdx(bend_x1, bend_z1)]);
		// }

		// if(gridCoordValid(bend_x2, bend_z2)) {
		// 	spring->nb_particles_.push_back(particles_[getParticleIdx(bend_x2, bend_z2)]);
		// }	
	}

	for(Particle* p : particles_) {
		std::cout << "particle spring number: " << p->springs_.size() << std::endl;
	}

	// update cache vertices
	refreshCache();

}


Cloth::~Cloth() {

}

void Cloth::tear(Spring* s) {
	Particle *p1 = s->p1_, *p2 = s->p2_;	// particles of current springs.
	std::cout << "to remove spring at " << glm::to_string(glm::vec2(p1->grid_x_, p1->grid_z_)) 
				<< ", " << glm::to_string(glm::vec2(p2->grid_x_, p2->grid_z_)) << std::endl;
	// Particle *nb_p1 = nullptr, *nb_p2 = nullptr;	// neighboring particles.(if any)
	// if(s->nb_particles_.size() >= 1) {
	// 	nb_p1 = s->nb_particles_[0];
	// }
	// if(s->nb_particles_.size() >= 2) {
	// 	nb_p2 = s->nb_particles_[1];
	// }

	Triangle *t1 = nullptr, *t2 = nullptr;	// neighboring triangles. (if any)
	if(s->triangles_.size() >= 1) {
		t1 = s->triangles_[0];
	}	
	if(s->triangles_.size() >= 2) {
		t2 = s->triangles_[1];
	}


	// std::cout << "triangle number: " << s->triangles_.size() << std::endl;
	// std::cout << "deleted spring position_: " << glm::to_string(s->p1_->position_) << " " 
	// 			<< glm::to_string(s->p2_->position_) << std::endl;

	// center of the teared spring.
	glm::vec3 init_center_position = (p1->init_position_ + p2->init_position_) / 2.0f;
	glm::vec3 curr_center_position = (p1->position_ + p2->position_) / 2.0f;
	glm::vec2 center_uv_coords = (p1->uv_coords_ + p2->uv_coords_) / 2.0f;

	// two new particles created because of the tearing
	glm::vec3 pp1_init_pos = init_center_position;
	glm::vec3 pp1_curr_pos = p1->position_ + (curr_center_position - p1->position_) * 0.75f;
	glm::vec2 pp1_uv_coords = center_uv_coords;
	Particle* pp1 = new Particle(pp1_init_pos, pp1_curr_pos, p1->mass_ / 2.0, pp1_uv_coords, true);
	particles_.push_back(pp1);
	Spring* ss1 = addStructSpring(p1, pp1, struct_k_, true);



	glm::vec3 pp2_init_pos = init_center_position;
	glm::vec3 pp2_curr_pos = p2->position_ + (curr_center_position - p2->position_) * 0.75f;
	glm::vec2 pp2_uv_coords = center_uv_coords;
	Particle* pp2 = new Particle(pp2_init_pos, pp2_curr_pos, p2->mass_ / 2.0, pp2_uv_coords, true);
	particles_.push_back(pp2);
	Spring* ss2 = addStructSpring(p2, pp2, struct_k_, true);
	


	if(t1) {
		Particle* nb_p1 = getNeighborParticle(t1, s);
		Triangle* tt1 = new Triangle(nb_p1, p1, pp1);
		Triangle* tt2 = new Triangle(nb_p1, pp2, p2);

		triangles_.insert(tt1);
		triangles_.insert(tt2);

		ss1->triangles_.push_back(tt1);
		ss2->triangles_.push_back(tt2);

		// if(!nb_p1->is_secondary_) {
			std::cout << "nb_p1 not secondary" << std::endl;
			Spring* ss11 = addStructSpring(pp1, nb_p1, struct_k_, true);
			Spring* ss12 = addStructSpring(pp2, nb_p1, struct_k_, true);
			ss11->triangles_.push_back(tt1);
			ss12->triangles_.push_back(tt2);
		// }
	

		getStructSpring(p1, nb_p1)->replaceTriangle(t1, tt1);
		getStructSpring(p2, nb_p1)->replaceTriangle(t1, tt2);
		triangles_.erase(t1);
		delete t1;
	}
	if(t2) {
		Particle* nb_p2 = getNeighborParticle(t2, s);
		Triangle* tt1 = new Triangle(p1, nb_p2, pp1);
		Triangle* tt2 = new Triangle(pp2, nb_p2, p2);

		triangles_.insert(tt1);
		triangles_.insert(tt2);

		ss1->triangles_.push_back(tt1);
		ss2->triangles_.push_back(tt2);

		// if(!nb_p2->is_secondary_) {
			std::cout << "nb_p2 not secondary" << std::endl;
			Spring* ss21 = addStructSpring(pp1, nb_p2, struct_k_, true);
			Spring* ss22 = addStructSpring(pp2, nb_p2, struct_k_, true);
			ss21->triangles_.push_back(tt1);
			ss22->triangles_.push_back(tt2);	
		// }
		
		getStructSpring(p1, nb_p2)->replaceTriangle(t2, tt1);
		getStructSpring(p2, nb_p2)->replaceTriangle(t2, tt2);
		triangles_.erase(t2);
		delete t2;

	}
	removeStructSpring(s);

}

Particle* Cloth::getNeighborParticle(Triangle* t1, Spring* s) {
	for(Particle* p : t1->particles_) {
		if(p != s->p1_ && p != s->p2_) {
			return p;
		}
	}
	std::cout << "get neighbor particle function correctly!" << std::endl;
	throw("get neighbor particle function correctly!");
	return nullptr;
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
	// std::cout << "currently have " << particles_.size() << "particles" << std::endl;

	std::vector<Particle*> splitted_particles;
	for(auto itr = particles_.begin(); itr != particles_.end(); itr++) {
		std::map<int, std::unordered_set<Particle*>> particle_groups;
		groupNeighbors(*itr, particle_groups);
		duplicateParticles(*itr, particle_groups, splitted_particles);
	}

	for(Particle* splitted_p : splitted_particles) {
		particles_.push_back(splitted_p);
	}

	for(Particle* particle : particles_) {
		particle->resetForce();
	}

	
	// std::unordered_set<Spring*> springs_to_remove;
	// for(Triangle* t : triangles_) {
	// 	Particle *p0 = t->particles_[0], *p1 = t->particles_[1], *p2 = t->particles_[2];
	// 	Spring *s0 = spring_map_[p0][p1], *s1 = spring_map_[p1][p2], *s2 = spring_map_[p2][p0];
	// 	if(shareSecondaryParticle(s0, s1)) {
	// 		springs_to_remove.insert(s0);
	// 		springs_to_remove.insert(s1);
	// 	}	
	// 	if(shareSecondaryParticle(s1, s2)) {
	// 		springs_to_remove.insert(s1);
	// 		springs_to_remove.insert(s2);
	// 	}
	// 	if(shareSecondaryParticle(s2, s0)) {
	// 		springs_to_remove.insert(s2);
	// 		springs_to_remove.insert(s0);
	// 	}
	// }
	// for(Spring* s : springs_to_remove) {
	// 	springs_.erase(s);
	// }

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
	// for(auto itr = springs_.begin(); itr != springs_.end(); itr++) {
	// 	Spring* struct_s = *itr;
	// 	if(struct_s->p1_->is_secondary_ && struct_s->p2_->is_secondary_) {
	// 		springs_.erase(struct_s);
	// 		continue;
	// 	}

	// 	struct_s->computeForceQuantity();
	// 	// TODO: if force quantity exceeds limit, break the spring
	// 	struct_s->applyForce();
		
	// 	if(struct_s->bend_spring_) {
	// 		struct_s->bend_spring_->computeForceQuantity();
	// 		struct_s->bend_spring_->applyForce();

	// 	}
	// }


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

	setCurrentSpring();
	// std::cout << "pick ray start: " << glm::to_string(pick_ray_start) << std::endl;
	if(picked_spring) {
		// std::cout << "spring selected" << std::endl;
		if(to_tear && !picked_spring->is_secondary_) {
			tear(picked_spring);
		}
		
	}
	else {
		// std::cout << "no spring selected" << std::endl;
	}


	refreshCache();
}

void Cloth::groupNeighbors(Particle* p, std::map<int, std::unordered_set<Particle*>>& groups) {
	std::vector<Particle*> nb_particles;
	for(Spring* s : p->springs_) {
		if(s->p1_ == nullptr || s->p2_ == nullptr) {
			std::cout << "spring has null node" << std::endl;
			throw "spring has null node";

		}
		if(!containsStructSpring(s->p1_, s->p2_)) {
			std::cout << "spring not recorded in map!" << std::endl;
			throw("spring not recorded in map!");
		}
		Particle* nb_particle = nullptr;
		if(s->p1_ == p) {
			nb_particle = s->p2_;
		}
		else if(s->p2_ == p) {
			nb_particle = s->p1_;
		}
		else {
			std::cout << "spring doesn't belong to current particle" << std::endl;
			throw("spring doesn't belong to current particle");
		}
		nb_particles.push_back(nb_particle);
	}
	// we need union-find here!
	std::vector<int> uf;
	uf.resize(nb_particles.size());
	for(int i = 0; i < uf.size(); i++) {
		uf[i] = i;	// if -1, the root of the group
	}
	for(int i = 0; i < nb_particles.size(); i++) {
		Particle *p1 = nb_particles[i];
		for(int j = i + 1; j < nb_particles.size(); j++) {
			Particle *p2 = nb_particles[j];
			if(containsStructSpring(p1, p2)) {
				// std::cout << "spring " << i << " connected to spring " << j << std::endl;
				uf[findRoot(uf, j)] = findRoot(uf, i);
			}
		}
	}

	for(int i = 0; i < uf.size(); i++) {
		if(uf[i] == i) {
			// std::cout << "gorup " << i << " created" << std::endl;
			groups[i] = std::unordered_set<Particle*>();
		}
	}	

	for(int i = 0; i < uf.size(); i++) {
		if(containsStructSpring(p, nb_particles[i])) {
			int group_num = findRoot(uf, i);
			groups[group_num].insert(nb_particles[i]);
		}
		else {
			std::cout << "neighbor spring grouped doesn't exist" << std::endl;
			throw("neighbor spring grouped doesn't exist");
		}
		
	}
}

int Cloth::findRoot(std::vector<int>& uf, int idx) {
	while(idx != uf[idx]) {
		idx = uf[idx];
	}
	return idx;
}

void Cloth::duplicateParticles(Particle* p, std::map<int, std::unordered_set<Particle*>>& groups, std::vector<Particle*> new_particles) {
	if(groups.size() <= 1) {
		return;
	} 
	int group_count = groups.size(); 
	for(auto const& group : groups) {
		if(group_count == 1) break;	// if only one group, don't need to split the origin particle
		group_count--;
		const std::unordered_set<Particle*>& group_particles = group.second;
		Particle* p_copy = new Particle(*p);
		new_particles.push_back(p_copy);
		for(Particle* nb_particle : group_particles) {
			Spring* s = getStructSpring(p, nb_particle);
			// replace the particle in the original spring
			p->springs_.erase(s);
			spring_map_[p][nb_particle] = nullptr;
			spring_map_[nb_particle][p] = nullptr;

			s->replaceParticle(p, p_copy);
			p_copy->springs_.insert(s);
			spring_map_[p_copy][nb_particle] = s;
			spring_map_[nb_particle][p_copy] = s;
			

			for(Triangle* t : s->triangles_) {	// replace the particle in old triangles. At most two triangles
				for(int p_idx = 0; p_idx < t->particles_.size(); p_idx++) {
					if(t->particles_[p_idx] == p) {
						t->particles_[p_idx] = p_copy;
					}
				}
			}

		}
		
	}
}




void Cloth::setCurrentSpring() {
	picked_spring = nullptr;
	float min_distance = std::numeric_limits<float>::max();
	for(Spring* s : springs_) {	// iterate all springs, and find the one with min distance
		float curr_distance = line_segment_distance(pick_ray_start, pick_ray_end, s->p1_->position_, s->p2_->position_);
		if(curr_distance < SPRING_CYLINDER_RADIUS && curr_distance < min_distance) {
			min_distance = curr_distance;
			picked_spring = s;
		}
	}
}


int Cloth::getParticleIdx(int x, int z) {
	return x * z_size_ + z;
}
bool Cloth::gridCoordValid(int x, int z) {
	return x >= 0 && x < x_size_ && z >= 0 && z < z_size_;
}

bool Cloth::containsStructSpring(Particle* p1, Particle* p2) {
	return (spring_map_[p1][p2] != nullptr) || (spring_map_[p2][p1] != nullptr);
}

Spring* Cloth::addStructSpring(Particle* p1, Particle* p2, float k, bool is_secondary) {
	if(containsStructSpring(p1, p2)) {
		std::cout << "the sprinig you want to create already exists!" << std::endl;
		throw("the sprinig you want to create already exists!");
	}
	Spring* s = new Spring(p1, p2, k, is_secondary);
	spring_map_[p1][p2] = s;
	spring_map_[p2][p1] = s;
	p1->springs_.insert(s);
	p2->springs_.insert(s);

	springs_.insert(s);
	return s;
}
Spring* Cloth::getStructSpring(Particle* p1, Particle* p2) {
	if(!containsStructSpring(p1, p2)) {
		std::cout << "struct spring doesn't exist!" << std::endl;
		throw("struct spring doesn't exist!");
	}
	return spring_map_[p1][p2];
}

void Cloth::removeStructSpring(Spring* s) {
	s->p1_->springs_.erase(s);
	s->p2_->springs_.erase(s);
	springs_.erase(s);
	spring_map_[s->p1_][s->p2_] = nullptr;
	spring_map_[s->p2_][s->p1_] = nullptr;
	// delete s;
}

bool Cloth::shareSecondaryParticle(Spring* s1, Spring* s2) {
	if((!s1) || (!s2)) return false;
	// return ((s1->p1_ == s2->p1_ && s1->p1_->is_secondary_)
	// 		|| (s1->p2_ == s2->p2_ && s1->p2_->is_secondary_)
	// 		|| (s1->p1_ == s2->p2_ && s1->p1_->is_secondary_)
	// 		|| (s1->p2_ == s2->p1_ && s1->p2_->is_secondary_)
	// 		);
	return false;
}
























	