R"zzz(#version 330 core

in vec3 vertex_position;

uniform mat4 projection;
uniform mat4 model;
uniform mat4 view;

void main() {
	mat4 mvp = projection * view * model;
	vec3 tmp_vertex_position = vertex_position + vec3(10.0, -20.0, 5.0);
	gl_Position = mvp * vec4(tmp_vertex_position, 1.0);
}
)zzz"
