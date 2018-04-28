#ifndef HELPER_FUNCTION_H
#define HELPER_FUNCTION_H

#define SMALL_NUM 0.00000001f   // used to avoid division overflow

#include <glm/glm.hpp>


float line_segment_distance(const glm::vec3& line1_start, const glm::vec3& line1_end, 
							const glm::vec3& line2_start, const glm::vec3& line2_end);


#endif