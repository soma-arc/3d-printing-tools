#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <memory>

#include "shape.h"
#include "args.hxx"

#define TINYOBJLOADER_IMPLEMENTATION
#include "./tiny_obj_loader.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

constexpr Object::Vertex rectangleVertex[] = {
    { -0.5f, -0.5f },
    { 0.5f, -0.5f },
    { 0.5f, 0.5f }
};

// LoadShader and LinkShader functions are from exrview
// https://github.com/syoyo/tinyexr/tree/master/examples/exrview
bool LoadShader(
    GLenum shaderType,  // GL_VERTEX_SHADER or GL_FRAGMENT_SHADER (or maybe GL_COMPUTE_SHADER)
    GLuint& shader,
    const char* shaderSourceFilename)
{
    GLint val = 0;

    // free old shader/program
    if (shader != 0) glDeleteShader(shader);

    static GLchar srcbuf[16384];
    FILE *fp = fopen(shaderSourceFilename, "rb");
    if (!fp) {
        fprintf(stderr, "failed to load shader: %s\n", shaderSourceFilename);
        return false;
    }
    fseek(fp, 0, SEEK_END);
    size_t len = ftell(fp);
    rewind(fp);
    len = fread(srcbuf, 1, len, fp);
    srcbuf[len] = 0;
    fclose(fp);

    static const GLchar *src = srcbuf;

    shader = glCreateShader(shaderType);
    glShaderSource(shader, 1, &src, NULL);
    glCompileShader(shader);
    glGetShaderiv(shader, GL_COMPILE_STATUS, &val);
    if (val != GL_TRUE) {
        char log[4096];
        GLsizei msglen;
        glGetShaderInfoLog(shader, 4096, &msglen, log);
        printf("%s\n", log);
        assert(val == GL_TRUE && "failed to compile shader");
    }

    printf("Load shader [ %s ] OK\n", shaderSourceFilename);
    return true;
}

bool LinkShader(
    GLuint& prog,
    GLuint& vertShader,
    GLuint& fragShader)
{
    GLint val = 0;

    if (prog != 0) {
        glDeleteProgram(prog);
    }

    prog = glCreateProgram();

    glAttachShader(prog, vertShader);
    glAttachShader(prog, fragShader);
    glLinkProgram(prog);

    glGetProgramiv(prog, GL_LINK_STATUS, &val);
    assert(val == GL_TRUE && "failed to link shader");

    printf("Link shader OK\n");

    return true;
}

int main(int argc, char** argv) {
    args::ArgumentParser parser("Generate mesh.");
    args::ValueFlag<float> sliceStep(parser, "step", "slice step", {'s', "sliceStep"}, 0.01f);
    args::ValueFlag<std::string> inputObj(parser, "obj", "Input obj file",
                                          {'i', "input"}, "scene.obj");
    args::ValueFlag<std::string> inputJson(parser, "json", "Input json file",
                                           {'j', "json"}, "scene.json");
    args::ValueFlag<std::string> outputBasenameArg(parser, "outputBasename", "Base name of output files (.vdb, .obj)",
                                                   {'o', "out"});
    args::Flag isFinite(parser, "isFinite", "Generate finite limit set", {'f', "finite"});
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});

    GLFWwindow* window;

    /* Initialize the library */
    if (!glfwInit())
        return -1;

    atexit(glfwTerminate);

    /* Create a windowed mode window and its OpenGL context */
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    // glfwWindowHint( GLFW_VISIBLE, 0 );
    window = glfwCreateWindow(640, 640, "Hello World", NULL, NULL);
    if (!window)
    {
        return -1;
    }

    /* Make the window's context current */
    glfwMakeContextCurrent(window);

    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK)
    {
        std::cerr << "Can't initialize GLEW" << std::endl;
        return 1;
    }


    GLuint vert_id = 0, frag_id = 0, prog_id = 0;
    {
        bool ret = LoadShader(GL_VERTEX_SHADER, vert_id, "./src/shader.vert");
        if (!ret) {
            fprintf(stderr, "Failed to compile vertex shader.\n");
            exit(-1);
        }
    }

    {
        bool ret = LoadShader(GL_FRAGMENT_SHADER, frag_id, "./src/shader.frag");
        if (!ret) {
            fprintf(stderr, "Failed to compile fragment shader.\n");
            exit(-1);
        }
    }

    {
        bool ret = LinkShader(prog_id, vert_id, frag_id);
        if (!ret) {
            fprintf(stderr, "Failed to link shader.\n");
            exit(-1);
        }
    }

    glBindAttribLocation(prog_id, 0, "position");
    glBindFragDataLocation(prog_id, 0, "fragment");
    glUseProgram(prog_id);

    glfwSwapInterval(1);
    glClearColor(0.f, 0.f, 0.f, 1.f);

    std::unique_ptr<const Shape> shape(new Shape(2, 4, rectangleVertex));

    /* Loop until the user closes the window */
    while (!glfwWindowShouldClose(window))
    {
        /* Render here */
        glClear(GL_COLOR_BUFFER_BIT);

        shape->draw();

        /* Swap front and back buffers */
        glfwSwapBuffers(window);

        /* Poll for and process events */
        glfwPollEvents();
    }

    return 0;
}
