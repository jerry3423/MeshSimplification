#include<glad/glad.h>
#include<GLFW/glfw3.h>
#include<glm/glm.hpp>
#include<glm/gtc/matrix_transform.hpp>
#include<glm/gtc/type_ptr.hpp>
#include"imgui/imgui.h"
#include"imgui/imgui_impl_glfw.h"
#include"imgui/imgui_impl_opengl3.h"

#define STB_IMAGE_IMPLEMENTATION
#include<stb_image.h>

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<ctime>

#include"Mesh.h"

// Local parameters for OpenGL rendering
unsigned int VAO, VBO_pos, VBO_tex, VBO_norm, EBO;
unsigned int vertexShader, fragShader, myShader;
unsigned int texture;
glm::mat4 model = glm::mat4(1.0f);
glm::mat4 view = glm::mat4(1.0f);
glm::mat4 projection = glm::mat4(1.0f);

int WIDTH = 1920; // Window width
int HEIGHT = 1080; // Window height

// Camera settings
glm::vec3 cameraPos = glm::vec3(0.0f, 0.0f, 0.3f);// Camera position
glm::vec3 cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);// Camera up vector
float fov = 45.0f;// Field of view
glm::vec2 mouseStart, mouseEnd;
bool isWireframeMode = false;
bool isDragging = false;
glm::quat rotation = glm::quat(1.0, 0.0, 0.0, 0.0);// Quaternion for arcball rotation
double simplifyTime = 0.0f;


// Local functions
void framebuffer_size_callback(GLFWwindow* window, int width, int height);

void renderMesh(Mesh& mesh);

void key_callback(GLFWwindow* window);

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);

void mousebutton_callback(GLFWwindow* window, int button, int action, int mods);

void mouse_callback(GLFWwindow* window, double xpos, double ypos);

void getShader();

void sceneTransform(const Mesh& mesh);

void initImgui(GLFWwindow* window);

void renderImgui(Mesh& mesh);

void loadTexture(const char* filePath);

glm::vec3 convertScreenToArcball(glm::vec2 pos);

void calculateBoundingBox(Mesh& mesh);


int main()
{
    // Initialize GLFW
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    // Create a window
    GLFWwindow* window = glfwCreateWindow(WIDTH, HEIGHT, "Mesh Simplification", NULL, NULL);
    if (!window)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);

    // Set GLFW callbacks for mouse and scroll events
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetMouseButtonCallback(window, mousebutton_callback);
    glfwSetCursorPosCallback(window, mouse_callback);

    // Initialize GLAD
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // Initialize imgui
    initImgui(window);

    // Enable depth testing for 3D rendering
    glEnable(GL_DEPTH_TEST);

    // Set viewport size
    glViewport(0, 0, WIDTH, HEIGHT);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    // Load and compile shaders
    getShader();

    // Set texture uniforms in shaders
    glUniform1i(glGetUniformLocation(myShader, "modelTexture"), 0);
    glUniform1i(glGetUniformLocation(myShader, "enableTexture"), GL_FALSE);

    // Create an empty mesh
    Mesh mesh;

    // Main rendering loop
    while (!glfwWindowShouldClose(window))
    {
        // Handle keyboard input
        key_callback(window);

        // Clear color and depth buffer
        glClearColor(0.8f, 1.0f, 1.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Activate texture unit and bind texture
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, texture);

        // Activate the shader program
        glUseProgram(myShader);

        // Camera control
        sceneTransform(mesh);

        // Render the mesh
        glBindVertexArray(VAO);
        glDrawElements(GL_TRIANGLES, mesh.faceVertices.size() * 3, GL_UNSIGNED_INT, 0);

        // render imgui
        renderImgui(mesh);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // Clean up
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO_pos);
    glDeleteBuffers(1, &VBO_tex);
    glDeleteBuffers(1, &VBO_norm);
    glDeleteBuffers(1, &EBO);
    glDeleteProgram(myShader);
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwTerminate();
    return 0;
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

void renderMesh(Mesh& mesh) {
    std::vector<float> newVertices;
    std::vector<float> newNormals;
    std::vector<float> newTexCoords;
    std::vector<int> indices;

    int currentIndex = 0;
    for (const auto& face : mesh.faceVertices) {
        // Get positions of the vertices of the face
        Eigen::Vector3f v0 = mesh.vertices[face.v[0]].position;
        Eigen::Vector3f v1 = mesh.vertices[face.v[1]].position;
        Eigen::Vector3f v2 = mesh.vertices[face.v[2]].position;

        // Calculate normal for each face
        Eigen::Vector3f n = (v1 - v0).cross(v2 - v0).normalized();

        for (int i = 0; i < 3; i++) {
            Eigen::Vector3f pos = mesh.vertices[face.v[i]].position;
            newVertices.push_back(pos.x());
            newVertices.push_back(pos.y());
            newVertices.push_back(pos.z());
            newNormals.push_back(n.x());
            newNormals.push_back(n.y());
            newNormals.push_back(n.z());
            // Check if there texture coordinates has been loaded
            if (mesh.hasTexture) {
                glm::vec2 tex = mesh.texCoords[face.v[i]];
                newTexCoords.push_back(tex.x);
                newTexCoords.push_back(tex.y);
            }
            else {
                newTexCoords.push_back(0.0f);
                newTexCoords.push_back(0.0f);
            }

            indices.push_back(currentIndex++);
        }
    }
    // Calculate the bounding box of the mesh
    calculateBoundingBox(mesh);

    // Get center of bounding box and translate model
    glm::vec3 center = (mesh.minBBX + mesh.maxBBX) / 2.0f;
    model = glm::translate(model, -center);

    // Generate VAO, VBO, EBO
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO_pos);
    glGenBuffers(1, &VBO_tex);
    glGenBuffers(1, &VBO_norm);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    // For vertex position
    glBindBuffer(GL_ARRAY_BUFFER, VBO_pos);
    glBufferData(GL_ARRAY_BUFFER, newVertices.size() * sizeof(float), newVertices.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // For texture coordinates
    glBindBuffer(GL_ARRAY_BUFFER, VBO_tex);
    glBufferData(GL_ARRAY_BUFFER, newTexCoords.size() * sizeof(float), newTexCoords.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1);

    // For normals
    glBindBuffer(GL_ARRAY_BUFFER, VBO_norm);
    glBufferData(GL_ARRAY_BUFFER, newNormals.size() * sizeof(float), newNormals.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(2);

    // For indices
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(int), indices.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    
}

void key_callback(GLFWwindow* window)
{
    if (ImGui::IsAnyItemActive()) {
        return;
    }

    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
    {
        glfwSetWindowShouldClose(window, true);
    }
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    if (fov >= 1.0f && fov <= 45.0f)
    {
        fov -= yoffset;
    }
    if (fov <= 1.0f)
    {
        fov = 1.0f;
    }
    if (fov >= 45.0f)
    {
        fov = 45.0f;
    }
}

void mousebutton_callback(GLFWwindow* window, int button, int action, int mods)
{
    // Use mouse left button to control arcball
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    {
        double xpos, ypos;
        // Get mouse position
        glfwGetCursorPos(window, &xpos, &ypos);
        mouseStart = glm::vec2(xpos, ypos);
        isDragging = true;
    }
    else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE)
    {
        isDragging = false;
        rotation = glm::quat(1.0, 0.0, 0.0, 0.0);
    }
}

void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (ImGui::IsAnyItemActive()) {
        return;
    }
    if (isDragging)
    {
        // Update the end position of the mouse
        mouseEnd = glm::vec2(xpos, ypos);

        // Check if there has been any movement since the last frame
        if (mouseStart != mouseEnd)
        {
            // Convert start and end mouse positions to vectors on the arcball sphere
            glm::vec3 vStart = convertScreenToArcball(mouseStart);
            glm::vec3 vCurr = convertScreenToArcball(mouseEnd);

            // Calculate the angle of rotation based on the arcball vectors
            float angle = acos(std::min(1.0f, glm::dot(vStart, vCurr)));

            // Get the axis of rotation
            glm::vec3 axisInCameraCoord = glm::cross(vStart, vCurr);
            glm::vec3 axisInWorldCoord = glm::inverse(glm::mat3(view)) * axisInCameraCoord;

            // Apply the rotation to the current rotation quaternion
            rotation = glm::rotate(rotation, angle, axisInWorldCoord);

            mouseStart = mouseEnd;
        }
    }
}

void getShader()
{
    std::ifstream vertexFile, fragFile;
    std::stringstream vertexStream, fragStream;

    // load shaders from file
    vertexFile.open("../Mesh Simplification/base.vert");
    fragFile.open("../Mesh Simplification/base.frag");

    // Read the shader source code into string streams
    vertexStream << vertexFile.rdbuf();
    fragStream << fragFile.rdbuf();

    // Convert string streams into strings
    std::string vertexCode = vertexStream.str();
    std::string fragCode = fragStream.str();

    const char* v = vertexCode.c_str();
    const char* f = fragCode.c_str();

    vertexFile.close();
    fragFile.close();

    // compile shaders
    vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &v, NULL);
    glCompileShader(vertexShader);

    fragShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragShader, 1, &f, NULL);
    glCompileShader(fragShader);

    // Create a shader program and attach the vertex and fragment shaders
    myShader = glCreateProgram();
    glAttachShader(myShader, vertexShader);
    glAttachShader(myShader, fragShader);
    // Link the shader program
    glLinkProgram(myShader);

    glDeleteShader(vertexShader);
    glDeleteShader(fragShader);
}

void sceneTransform(const Mesh& mesh)
{
    // Apply arcball rotation to the model
    if (isDragging)
    {
        glm::mat4 rotMatrix = glm::mat4_cast(rotation);
        model = rotMatrix * model;

        // Reset the rotation quaternion for the next frame
        rotation = glm::quat(1.0, 0.0, 0.0, 0.0);
    }

    view = glm::lookAt(cameraPos, glm::vec3(0.0, 0.0, 0.0), cameraUp);
    projection = glm::perspective(glm::radians(fov), float(WIDTH) / float(HEIGHT), 0.01f, 100.0f);

    // Pass the transformation matrices to the shader
    glUniformMatrix4fv(glGetUniformLocation(myShader, "model"), 1, GL_FALSE, glm::value_ptr(model));
    glUniformMatrix4fv(glGetUniformLocation(myShader, "view"), 1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(glGetUniformLocation(myShader, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
}

void initImgui(GLFWwindow* window)
{
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.FontGlobalScale = 1.5f;
    ImGui::StyleColorsDark();// Use dark color scheme

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");
}

void renderImgui(Mesh& mesh)
{
    // Create new imgui frame
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    // Create imgui window
    ImGui::Begin("Editor");

    // Camera Control
    if (ImGui::CollapsingHeader("Camera"))
    {
        // Control camera position
        ImGui::InputFloat3("Camera Position", &cameraPos[0]);
    }

    // Model operations here: load and save
    if (ImGui::CollapsingHeader("Model Operation"))
    {
        // Input filename
        ImGui::InputText("FileName", mesh.fileName, IM_ARRAYSIZE(mesh.fileName));

        // Load Obj file with filename
        if (ImGui::Button("Load Model"))
        {
            // Clear mesh data first
            mesh.clear();
            simplifyTime = 0.0f;
            // Clear texture
            glDeleteTextures(1, &texture);
            texture = 0;
            glUniform1i(glGetUniformLocation(myShader, "enableTexture"), GL_FALSE);
            mesh.enableTexture = false;
            // Load file and render mesh
            mesh.readFile(mesh.fileName);
            renderMesh(mesh);
        }
        ImGui::SameLine();

        // Save Obj file using filename
        if (ImGui::Button("Save Model"))
        {
            mesh.writeFile(mesh.fileName);
        }
        ImGui::SameLine();

        // Load texture
        if (ImGui::Button("Load Texture"))
        {
            // Clear texture
            if (texture != 0)
            {
                glDeleteTextures(1, &texture);
                texture = 0;
            }

            // Load texture according to filename
            std::string baseFilePath = "../models/";
            std::string extension = ".jpg";
            std::string filePath = baseFilePath + mesh.fileName + extension;
            loadTexture(filePath.c_str());

            glUniform1i(glGetUniformLocation(myShader, "enableTexture"), GL_TRUE);
            mesh.enableTexture = true;
        }

        // Present texture or not
        ImGui::Checkbox("Texture", &mesh.enableTexture);
        glUniform1i(glGetUniformLocation(myShader, "enableTexture"), mesh.enableTexture ? GL_TRUE : GL_FALSE);


        // Polygon Mode for debug
        if (ImGui::Checkbox("Wireframe Mode", &isWireframeMode)) {
            if (isWireframeMode) {
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            }
            else {
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            }
        }

    }

    // Mesh simplification control
    if (ImGui::CollapsingHeader("Simplify Control"))
    {
        ImGui::Combo("Vertex Placement", &mesh.currentMode, mesh.SimplifyMode, IM_ARRAYSIZE(mesh.SimplifyMode));

        // Control simplify ratio
        ImGui::SliderFloat("Simplify Ratio", &mesh.simplifyRatio, 0.0f, 1.0f, "%.2f");

        // Aggregation control
        ImGui::Checkbox("Aggregation", &mesh.enableAggregation);

        if (mesh.enableAggregation)
        {
            // Threshold for aggregation
            ImGui::InputFloat("Distance Threshold", &mesh.t);
        }

        // Boundary preservation
        ImGui::Checkbox("Preserve Boundary", &mesh.preserveBoundary);

        // Prevent mesh inversion
        ImGui::Checkbox("Prevent Inversion", &mesh.preventInversion);

        // Start simplify
        if (ImGui::Button("Simplify")) {

            std::clock_t start = std::clock();
            mesh.preProcessData();
            // Do mesh simplification
            mesh.simplify();
            std::clock_t end = std::clock();
            simplifyTime = double(end - start) / CLOCKS_PER_SEC;
            renderMesh(mesh);
        }

        // Show simplify time
        ImGui::Text("Simplification time: %.3f seconds", simplifyTime);
        // Show vertex number
        ImGui::Text("Vertex number: %d ", mesh.vertices.size());
        // Show face number
        ImGui::Text("Face number: %d ", mesh.faceVertices.size());
    }
    ImGui::End();

    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

void loadTexture(const char* filePath)
{
    int width, height, channels;

    // Flip the image vertically
    stbi_set_flip_vertically_on_load(true);
    // Load the image from file path
    unsigned char* data = stbi_load(filePath, &width, &height, &channels, 0);

    // Check if the image data was loaded successfully
    if (data) {
        // Generate an OpenGL texture object
        glGenTextures(1, &texture);
        glBindTexture(GL_TEXTURE_2D, texture);

        // Set texture wrapping parameters
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

        // Set texture filtering parameters
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        // Load the image into the texture object
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);

        // Free image data
        stbi_image_free(data);
    }
    else {
        std::cerr << "Failed to load texture" << std::endl;
    }
}

// Converts 2D mouse coordinates into a 3D vector on an arcball
glm::vec3 convertScreenToArcball(glm::vec2 pos)
{
    // Normalize the 2D mouse position to [-1, 1] coordinates
    glm::vec3 a = glm::vec3(1.0 * pos.x / WIDTH * 2 - 1.0, 1.0 * pos.y / HEIGHT * 2 - 1.0, 0);
    a.y = -a.y;

    // Compute the square of the length of the vector
    float len = a.x * a.x + a.y * a.y;

    if (len <= 1)
    {
        // If the point is inside the arcball, calculate the z-coordinate
        a.z = sqrt(1 - len);
    }
    else
    {
        // If the point is outside the arcball, normalize the vector
        a = glm::normalize(a);
    }

    return a;
}

void calculateBoundingBox(Mesh& mesh)
{
    // Reset the model, view, and projection matrices
    model = glm::mat4(1.0f);
    view = glm::mat4(1.0f);
    projection = glm::mat4(1.0f);

    // Initialize the bounding box
    mesh.minBBX = glm::vec3(mesh.vertices[0].position.x(), mesh.vertices[0].position.y(), mesh.vertices[0].position.z());
    mesh.maxBBX = glm::vec3(mesh.vertices[0].position.x(), mesh.vertices[0].position.y(), mesh.vertices[0].position.z());

    // Iterate over each vertex in the mesh to expand the bounding box
    for (const Vertex& vertex : mesh.vertices) {
        mesh.minBBX.x = std::min(mesh.minBBX.x, vertex.position.x());
        mesh.minBBX.y = std::min(mesh.minBBX.y, vertex.position.y());
        mesh.minBBX.z = std::min(mesh.minBBX.z, vertex.position.z());

        mesh.maxBBX.x = std::max(mesh.maxBBX.x, vertex.position.x());
        mesh.maxBBX.y = std::max(mesh.maxBBX.y, vertex.position.y());
        mesh.maxBBX.z = std::max(mesh.maxBBX.z, vertex.position.z());
    }
}
