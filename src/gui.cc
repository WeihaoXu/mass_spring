#include "gui.h"
#include "config.h"
#include <jpegio.h>

#include <iostream>
#include <debuggl.h>
#include <glm/gtc/matrix_access.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtx/string_cast.hpp>



GUI::GUI(GLFWwindow* window)
	:window_(window)
{
	glfwSetWindowUserPointer(window_, this);
	glfwSetKeyCallback(window_, KeyCallback);
	glfwSetCursorPosCallback(window_, MousePosCallback);
	glfwSetMouseButtonCallback(window_, MouseButtonCallback);

	glfwGetWindowSize(window_, &window_width_, &window_height_);
	float aspect_ = static_cast<float>(window_width_) / window_height_;
	projection_matrix_ = glm::perspective((float)(kFov * (M_PI / 180.0f)), aspect_, kNear, kFar);
}

GUI::~GUI()
{
}


void GUI::keyCallback(int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
		glfwSetWindowShouldClose(window_, GL_TRUE);
		return ;
	}
	if (key == GLFW_KEY_J && action == GLFW_RELEASE) {
		toggleDrawSpring();
	}

	if (key == GLFW_KEY_K && action == GLFW_RELEASE) {
		toggleDrawCloth();
	}
	if (key == GLFW_KEY_MINUS && action != GLFW_RELEASE) {
		wind_factor_ /= 1.1f;
		wind_factor_ = std::max(wind_factor_, 0.5f);
	}
	if (key == GLFW_KEY_EQUAL && action != GLFW_RELEASE) {
		wind_factor_ *= 1.1f;
		wind_factor_ = std::min(wind_factor_, 1.5f);
	} 
	if (key == GLFW_KEY_LEFT_CONTROL || key == GLFW_KEY_RIGHT_CONTROL) {
		if(action == GLFW_PRESS) {
			control_pressed_ = true;
		}
		else {
			control_pressed_ = false;
		}
	}


	if (captureWASDUPDOWN(key, action))
		return ;
	if (key == GLFW_KEY_LEFT || key == GLFW_KEY_RIGHT) {
		float roll_speed;
		if (key == GLFW_KEY_RIGHT)
			roll_speed = -roll_speed_;
		else
			roll_speed = roll_speed_;
		// FIXME: actually roll the bone here
	} else if (key == GLFW_KEY_C && action != GLFW_RELEASE) {
		fps_mode_ = !fps_mode_;
	} else if (key == GLFW_KEY_LEFT_BRACKET && action == GLFW_RELEASE) {
		time_speed_ /= 2.0;
	} else if (key == GLFW_KEY_RIGHT_BRACKET && action == GLFW_RELEASE) {
		time_speed_ *= 2.0;
	} else if (key == GLFW_KEY_T && action != GLFW_RELEASE) {
		cloth_->enable_wind = !cloth_->enable_wind;
	} else if (key == GLFW_KEY_O && action != GLFW_RELEASE) {
		cloth_->enable_sphere = !cloth_->enable_sphere;
	} else if (key == GLFW_KEY_I && action != GLFW_RELEASE) {
		cloth_->pause_sphere = !cloth_->pause_sphere;
	} else if (key == GLFW_KEY_R && action != GLFW_RELEASE) {
		reset_ms_system_ = true;
	} else if (key == GLFW_KEY_P && action != GLFW_RELEASE) {
		to_toggle_wind_direct_ = true;
	}
}

void GUI::mousePosCallback(double mouse_x, double mouse_y)
{
	last_x_ = current_x_;
	last_y_ = current_y_;
	current_x_ = mouse_x;
	current_y_ = window_height_ - mouse_y;
	float delta_x = current_x_ - last_x_;
	float delta_y = current_y_ - last_y_;
	if (sqrt(delta_x * delta_x + delta_y * delta_y) < 1e-15)
		return;
	glm::vec3 mouse_direction = glm::normalize(glm::vec3(delta_x, delta_y, 0.0f));
	glm::vec2 mouse_start = glm::vec2(last_x_, last_y_);
	glm::vec2 mouse_end = glm::vec2(current_x_, current_y_);
	glm::uvec4 viewport = glm::uvec4(0, 0, window_width_, window_height_);

	bool drag_camera = drag_state_ && current_button_ == GLFW_MOUSE_BUTTON_RIGHT;
	bool tear_particle = (!control_pressed_) && drag_state_ && current_button_ == GLFW_MOUSE_BUTTON_LEFT;
	bool drag_particle = control_pressed_ && drag_state_ && current_button_ == GLFW_MOUSE_BUTTON_LEFT;

	if (drag_camera) {
		glm::vec3 axis = glm::normalize(
				orientation_ *
				glm::vec3(mouse_direction.y, -mouse_direction.x, 0.0f)
				);
		orientation_ =
			glm::mat3(glm::rotate(rotation_speed_, axis) * glm::mat4(orientation_));
		tangent_ = glm::column(orientation_, 0);
		up_ = glm::column(orientation_, 1);
		look_ = glm::column(orientation_, 2);
	}

	// std::cout << "mouse position: " << glm::to_string(glm::vec2(mouse_x, mouse_y)) << std::endl;
	glm::vec3 mouse_pos = glm::unProject(glm::vec3(current_x_, current_y_, 1.0f),
											view_matrix_,
											projection_matrix_,
											viewport);
	// std::cout << "current position in world coords: (" << mouse_pos.x << ", " << mouse_pos.y << ", " << mouse_pos.z << ")" << std::endl;
	
	// pick spring, similar to bone picking in animation project.
	{
		glm::vec3 pick_ray_direct = glm::normalize(mouse_pos - eye_);
		glm::vec3 pick_ray_end = eye_ + PICK_RAY_LEN * pick_ray_direct;
		cloth_->pick_ray_start = eye_;
		cloth_->pick_ray_end = pick_ray_end;
	}


	if(tear_particle) {
		cloth_->to_tear = true;
	}
	else {
		cloth_->to_tear = false;
	}

	if(drag_particle) {
		glm::vec3 mouse_pos = glm::unProject(glm::vec3(current_x_, current_y_, 1.0f),
											view_matrix_,
											projection_matrix_,
											viewport);
		glm::vec3 last_mouse_pos_ = glm::unProject(glm::vec3(last_x_, last_y_, 1.0f),
											view_matrix_,
											projection_matrix_,
											viewport);	

		glm::vec3 darg_dist = mouse_pos - last_mouse_pos_;
		Particle* p = cloth_->getCurrentParticle();
		if(p) {
			p->move(darg_dist);
		}
	}

	
}

void GUI::mouseButtonCallback(int button, int action, int mods)
{
	drag_state_ = (action == GLFW_PRESS);
	current_button_ = button;
}

void GUI::updateMatrices()
{
	// Compute our view, and projection matrices.
	if (fps_mode_)
		center_ = eye_ + camera_distance_ * look_;
	else
		eye_ = center_ - camera_distance_ * look_;

	view_matrix_ = glm::lookAt(eye_, center_, up_);
	light_position_ = glm::vec4(eye_, 1.0f);

	aspect_ = static_cast<float>(window_width_) / window_height_;
	projection_matrix_ =
		glm::perspective((float)(kFov * (M_PI / 180.0f)), aspect_, kNear, kFar);
	model_matrix_ = glm::mat4(1.0f);
}

MatrixPointers GUI::getMatrixPointers() const
{
	MatrixPointers ret;
	ret.projection = &projection_matrix_[0][0];
	ret.model= &model_matrix_[0][0];
	ret.view = &view_matrix_[0][0];
	return ret;
}




bool GUI::captureWASDUPDOWN(int key, int action)
{
	if (key == GLFW_KEY_W) {
		if (fps_mode_)
			eye_ += zoom_speed_ * look_;
		else
			camera_distance_ -= zoom_speed_;
		return true;
	} else if (key == GLFW_KEY_S) {
		if (fps_mode_)
			eye_ -= zoom_speed_ * look_;
		else
			camera_distance_ += zoom_speed_;
		return true;
	} else if (key == GLFW_KEY_A) {
		if (fps_mode_)
			eye_ -= pan_speed_ * tangent_;
		else
			center_ -= pan_speed_ * tangent_;
		return true;
	} else if (key == GLFW_KEY_D) {
		if (fps_mode_)
			eye_ += pan_speed_ * tangent_;
		else
			center_ += pan_speed_ * tangent_;
		return true;
	} else if (key == GLFW_KEY_DOWN) {
		if (fps_mode_)
			eye_ -= pan_speed_ * up_;
		else
			center_ -= pan_speed_ * up_;
		return true;
	} else if (key == GLFW_KEY_UP) {
		if (fps_mode_)
			eye_ += pan_speed_ * up_;
		else
			center_ += pan_speed_ * up_;
		return true;
	}
	return false;
}

void dragParticle() {

}

// Delegrate to the actual GUI object.
void GUI::KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	GUI* gui = (GUI*)glfwGetWindowUserPointer(window);
	gui->keyCallback(key, scancode, action, mods);
}

void GUI::MousePosCallback(GLFWwindow* window, double mouse_x, double mouse_y)
{
	GUI* gui = (GUI*)glfwGetWindowUserPointer(window);
	gui->mousePosCallback(mouse_x, mouse_y);
}

void GUI::MouseButtonCallback(GLFWwindow* window, int button, int action, int mods)
{
	GUI* gui = (GUI*)glfwGetWindowUserPointer(window);
	gui->mouseButtonCallback(button, action, mods);
}
