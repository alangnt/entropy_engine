#include <GL/glew.h>
#include <GLFW/glfw3.h>

void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
    glfwSetWindowShouldClose(window, GLFW_TRUE);
  }
}

int main() {
  GLFWwindow* window;

  // Initialize the library
  if (!glfwInit())
    return 1;

  // Create a windowed mode window and its OpenGL context
  window = glfwCreateWindow(1280, 720, "Entropy Engine - GPU Core", NULL, NULL);
  if (!window) {
    glfwTerminate();
    return -1;
  }

  // Make the window's context current
  glfwMakeContextCurrent(window);

  // Init GLEW
  glewInit();

  // Background color
  glClearColor(20.0f/255.0f, 20.0f/255.0f, 20.0f/255.0f, 1.0f);

  // Listen to the keyboard input
  glfwSetKeyCallback(window, keyCallback);

  float vertices[] = { 
    0.2f, 0.2f, 
    0.4f, 0.2f, 
    0.3f, 0.4f 
  };

  // Loop until we close the window
  while (!glfwWindowShouldClose(window)) {
    // Render
    glClear(GL_COLOR_BUFFER_BIT);

    // Swap front and back buffers
    glfwSwapBuffers(window);

    // Poll for and process events
    glfwPollEvents();
  }

  glfwTerminate();
  return 0;
}