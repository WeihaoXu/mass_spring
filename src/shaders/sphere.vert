R"zzz(#version 330 core

in vec3 vertex_position;

uniform mat4 projection;
uniform mat4 model;
uniform mat4 view;
uniform vec3 sphere_position;

void main() {
	mat4 mvp = projection * view * model;
	vec3 tmp_vertex_position = vertex_position + sphere_position;
	gl_Position = mvp * vec4(tmp_vertex_position, 1.0);
}
)zzz"
