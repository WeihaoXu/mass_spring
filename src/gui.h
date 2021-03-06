#ifndef SKINNING_GUI_H
#define SKINNING_GUI_H

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <GLFW/glfw3.h>
#include "cloth.h"

#define PICK_RAY_LEN 200.0f

struct Mesh;


/*
 * Hint: call glUniformMatrix4fv on thest pointers
 */
struct MatrixPointers {
	const float *projection, *model, *view;
};

class GUI {
public:
	GUI(GLFWwindow*);
	~GUI();

	void keyCallback(int key, int scancode, int action, int mods);
	void mousePosCallback(double mouse_x, double mouse_y);
	void mouseButtonCallback(int button, int action, int mods);
	void updateMatrices();
	void assignCloth(Cloth* cloth) {cloth_ = cloth;}
	MatrixPointers getMatrixPointers() const;

	static void KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
	static void MousePosCallback(GLFWwindow* window, double mouse_x, double mouse_y);
	static void MouseButtonCallback(GLFWwindow* window, int button, int action, int mods);

	glm::vec3 getCenter() const { return center_; }
	const glm::vec3& getCamera() const { return eye_; }
	bool isPoseDirty() const { return pose_changed_; }
	void clearPose() { pose_changed_ = false; }

	
	bool toResetSystem() {return reset_ms_system_;}
	void clearResetFlag() {reset_ms_system_ = false;}

	bool drawClothEnabled() {return enable_draw_cloth_;}
	void toggleDrawCloth() {enable_draw_cloth_ = !enable_draw_cloth_;}

	bool drawSphereEnabled() {return cloth_->enable_sphere; }
	bool drawSpringEnabled() {return enable_draw_spring_;}
	void toggleDrawSpring() {enable_draw_spring_ = !enable_draw_spring_;}

	float getWindFactor() {return wind_factor_;}
	
	float getTimeSpeed() {return time_speed_;}
	
	bool toToggleWindDirect() {return to_toggle_wind_direct_;}
	void clearToggleWindDirectFlag() {to_toggle_wind_direct_ = false;}



	const float* getLightPositionPtr() const { return &light_position_[0]; }
	
	// int getCurrentBone() const { return current_bone_; }
	// const int* getCurrentBonePointer() const { return &current_bone_; }
	// bool setCurrentBone(int i);

	bool isTransparent() const { return transparent_; }
private:
	GLFWwindow* window_;
	// Mesh* mesh_;
	Cloth* cloth_;

	int window_width_, window_height_;

	bool to_toggle_wind_direct_ = false;
	bool reset_ms_system_ = false;
	bool enable_draw_cloth_ = true;
	bool enable_draw_spring_ = false;

	float wind_factor_ = 1.0f;

	bool control_pressed_ = false;
	bool drag_state_ = false;
	bool fps_mode_ = false;
	bool pose_changed_ = true;
	bool transparent_ = false;
	// int current_bone_ = -1;
	int current_button_ = -1;
	float roll_speed_ = M_PI / 64.0f;
	float last_x_ = 0.0f, last_y_ = 0.0f, current_x_ = 0.0f, current_y_ = 0.0f;
	float camera_distance_ = 20.0;
	float pan_speed_ = 0.1f * 3;
	float rotation_speed_ = 0.02f * 2;
	float zoom_speed_ = 0.1f * 3;
	float aspect_;

	float time_speed_ = 1.0;

	glm::vec3 eye_ = glm::vec3(0.0f, 0.0f, camera_distance_);
	glm::vec3 up_ = glm::vec3(0.0f, 1.0f, 0.0f);
	glm::vec3 look_ = glm::vec3(0.0f, 0.0f, -1.0f);
	glm::vec3 tangent_ = glm::cross(look_, up_);
	glm::vec3 center_ = eye_ - camera_distance_ * look_;
	glm::mat3 orientation_ = glm::mat3(tangent_, up_, look_);
	glm::vec4 light_position_;

	glm::mat4 view_matrix_ = glm::lookAt(eye_, center_, up_);
	glm::mat4 projection_matrix_;
	glm::mat4 model_matrix_ = glm::mat4(1.0f);

	bool captureWASDUPDOWN(int key, int action);
	
};
#endif
