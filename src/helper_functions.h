#ifndef HELPER_FUNCTION_H
#define HELPER_FUNCTION_H

#define SMALL_NUM 0.00000001f   // used to avoid division overflow

#include <vector>
#include <glm/glm.hpp>


float line_segment_distance(const glm::vec3& line1_start, const glm::vec3& line1_end, 
							const glm::vec3& line2_start, const glm::vec3& line2_end);

void create_floor(std::vector<glm::vec4>& floor_vertices, std::vector<glm::uvec3>& floor_faces);
float line_point_distance(glm::vec3& line_start, glm::vec3& line_end, glm::vec3& point);

void create_sphere(std::vector<glm::vec3>& sphere_vertex, std::vector<glm::vec3>& sphere_normal, std::vector<glm::uvec3>& sphere_indices);

#endif