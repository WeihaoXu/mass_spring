#include "helper_functions.h"
#include "config.h"
#include <glm/glm.hpp>
#include <math.h>
#include <iostream>
#include <glm/gtx/string_cast.hpp>
#include <GL/glew.h>

// Compute minimum distance between two line segments. Reference: http://geomalgorithms.com/a07-_distance.html
//    Input:  start and end points of two line segments
//    Return: the shortest distance between two lines.
float line_segment_distance(const glm::vec3& line1_start, const glm::vec3& line1_end, 
							const glm::vec3& line2_start, const glm::vec3& line2_end)
{
	glm::vec3 u = line1_end - line1_start;
	glm::vec3 v = line2_end - line2_start;
	glm::vec3 w = line1_start - line2_start;

    float    a = glm::dot(u,u);         // always >= 0
    float    b = glm::dot(u,v);
    float    c = glm::dot(v,v);         // always >= 0
    float    d = glm::dot(u,w);
    float    e = glm::dot(v,w);
    float    D = a*c - b*b;        // always >= 0
    float    sc, sN, sD = D;       // sc = sN / sD, default sD = D >= 0
    float    tc, tN, tD = D;       // tc = tN / tD, default tD = D >= 0

    // compute the line parameters of the two closest points
    if (D < SMALL_NUM) { // the lines are almost parallel
        sN = 0.0;         // force using point P0 on segment S1
        sD = 1.0;         // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    }
    else {                 // get the closest points on the infinite lines
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if (sN < 0.0) {        // sc < 0 => the s=0 edge is visible
            sN = 0.0;
            tN = e;
            tD = c;
        }
        else if (sN > sD) {  // sc > 1  => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }

    if (tN < 0.0) {            // tc < 0 => the t=0 edge is visible
        tN = 0.0;
        // recompute sc for this edge
        if (-d < 0.0)
            sN = 0.0;
        else if (-d > a)
            sN = sD;
        else {
            sN = -d;
            sD = a;
        }
    }
    else if (tN > tD) {      // tc > 1  => the t=1 edge is visible
        tN = tD;
        // recompute sc for this edge
        if ((-d + b) < 0.0)
            sN = 0;
        else if ((-d + b) > a)
            sN = sD;
        else {
            sN = (-d +  b);
            sD = a;
        }
    }
    // finally do the division to get sc and tc
    sc = (std::abs(sN) < SMALL_NUM ? 0.0 : sN / sD);
    tc = (std::abs(tN) < SMALL_NUM ? 0.0 : tN / tD);

    // get the difference of the two closest points
    glm::vec3   dP = w + (sc * u) - (tc * v);  // =  S1(sc) - S2(tc)

    return glm::length(dP);   // return the closest distance

}

void create_floor(std::vector<glm::vec4>& floor_vertices, std::vector<glm::uvec3>& floor_faces)
{
    floor_vertices.push_back(glm::vec4(kFloorXMin, kFloorY, kFloorZMax, 1.0f));
    floor_vertices.push_back(glm::vec4(kFloorXMax, kFloorY, kFloorZMax, 1.0f));
    floor_vertices.push_back(glm::vec4(kFloorXMax, kFloorY, kFloorZMin, 1.0f));
    floor_vertices.push_back(glm::vec4(kFloorXMin, kFloorY, kFloorZMin, 1.0f));
    floor_faces.push_back(glm::uvec3(0, 1, 2));
    floor_faces.push_back(glm::uvec3(2, 3, 0));
}

float line_point_distance(glm::vec3& line_start, glm::vec3& line_end, glm::vec3& point) {
    // std::cout << "start: " << glm::to_string(line_start) << ", end: " << glm::to_string(line_end) << std::endl;
    glm::vec3 sp = point - line_start;
    glm::vec3 se = line_end - line_start;
    float prj_len = glm::dot(glm::normalize(sp), glm::normalize(se)) * glm::length(sp);
    float sp_len = glm::length(sp);

    // std::cout << "len1: " << prj_len << ", len2: " << sp_len << std::endl;

    float res = sqrt(sp_len * sp_len - prj_len * prj_len);
    // std::cout << "line point distance: " << res << std::endl;
    return res;
}

void create_sphere(std::vector<glm::vec3>& sphere_vertex, std::vector<glm::vec3>& sphere_normal, std::vector<glm::uvec3>& sphere_indices){
    const int na=36;        // vertex grid size
    const int nb=18;
    GLfloat x,y,z,a,b,da,db,r=3.5;
    int ia,ib,ix,iy;
    da=2.0*M_PI/GLfloat(na);
    db=    M_PI/GLfloat(nb-1);
    // [Generate sphere point data]
    // spherical angles a,b covering whole sphere surface
    for (ix=0,b=-0.5*M_PI,ib=0;ib<nb;ib++,b+=db){
        for (a=0.0,ia=0;ia<na;ia++,a+=da,ix+=3){
            // unit sphere
            x=cos(b)*cos(a);
            y=cos(b)*sin(a);
            z=sin(b);
            glm::vec3 position = glm::vec3(x*r, y*r, z*r);
            glm::vec3 norm = glm::vec3(x, y, z);
            sphere_vertex.push_back(position);
            sphere_normal.push_back(norm);
        }
    }
    // [Generate GL_TRIANGLE indices]
    for (ix=0,iy=0,ib=1;ib<nb;ib++) {
        for (ia=1;ia<na;ia++,iy++){
            // first half of QUAD
            sphere_indices.push_back(glm::vec3(iy, iy+1,iy+na));
            // second half of QUAD
            sphere_indices.push_back(glm::vec3(iy+na, iy+1,iy+na+1));
        }
        // first half of QUAD
        sphere_indices.push_back(glm::vec3(iy, iy+1-na,iy+na));
        // second half of QUAD
        sphere_indices.push_back(glm::vec3(iy+na, iy+1-na,iy+1));
        iy++;
    }
}
