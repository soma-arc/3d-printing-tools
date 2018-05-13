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
#include "vec3.h"

#define TINYOBJLOADER_IMPLEMENTATION
#include "./tiny_obj_loader.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

constexpr Object::Vertex rectangleVertex[] = {
    { -0.5f, -0.5f, 0.0f, 0.0f, 0.0f },
    { 0.5f, -0.5f, 0.0f, 0.0f, 0.8f },
    { 0.5f, 0.5f, 0.0f, 0.8f, 0.0f }
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

// void LoadObj(std::string objFilename) {
//     tinyobj::attrib_t attrib;
//     std::vector<tinyobj::shape_t> shapes;
//     std::vector<tinyobj::material_t> materials;

//     std::string err;

//     std::cout << "Load " << objFilename << std::endl;
    
//     bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, objFilename);

//     if (!err.empty()) { // `err` may contain warning message.
//         std::cerr << err << std::endl;
//     }

//     if (!ret) {
//         std::cerr << "Failed to load " << objFilename << std::endl;
//         exit(1);
//     }

//     printf("# of vertices  = %d\n", (int) (attrib.vertices.size()) / 3);
//     printf("# of normals   = %d\n", (int) (attrib.normals.size()) / 3);
//     printf("# of texcoords = %d\n", (int) (attrib.texcoords.size()) / 2);
//     printf("# of materials = %d\n", (int) materials.size());
//     printf("# of shapes    = %d\n", (int) shapes.size());
//     printf("\n");

//     unsigned char *textureData = new unsigned char[texSize * texSize * 3];
//     for (size_t s = 0; s < shapes.size(); s++) {
//         printf("Shape %d\n", (int) shapes.size());
//         printf("# of faces %d\n", (int) shapes[s].mesh.num_face_vertices.size());
//         size_t index_offset = 0;
//         for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
//             // The number of vertexes in the face. 3, 4 or 5?
//             int fv = shapes[s].mesh.num_face_vertices[f];
//             if(fv > 3) {
//                 printf("quad!!");
//                 continue;
//             }

//             std::vector<Vec3f> faceVert;
//             std::vector<Vec3f> faceUv;
//             Vec3f uvMin(2, 0, 2);
//             Vec3f uvMax(-1, 0, -1);
//             // Loop over vertices in the face.
//             for (int v = 0; v < fv; v++) {
//                 // access to vertex
//                 tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
//                 tinyobj::real_t vx = attrib.vertices[3*idx.vertex_index+0];
//                 tinyobj::real_t vy = attrib.vertices[3*idx.vertex_index+1];
//                 tinyobj::real_t vz = attrib.vertices[3*idx.vertex_index+2];
//                 tinyobj::real_t nx = attrib.normals[3*idx.normal_index+0];
//                 tinyobj::real_t ny = attrib.normals[3*idx.normal_index+1];
//                 tinyobj::real_t nz = attrib.normals[3*idx.normal_index+2];
//                 tinyobj::real_t tx = attrib.texcoords[2*idx.texcoord_index+0];
//                 tinyobj::real_t ty = attrib.texcoords[2*idx.texcoord_index+1];
//                 faceVert.push_back(Vec3f(vx, vy, vz));
//                 faceUv.push_back(Vec3f(tx, 0, ty));
//                 uvMin[0] = std::min(tx, uvMin[0]);
//                 uvMin[2] = std::min(ty, uvMin[2]);
//                 uvMax[0] = std::max(tx, uvMax[0]);
//                 uvMax[2] = std::max(ty, uvMax[2]);
//             }
//             index_offset += fv;

//             printf("(%f, %f), (%f, %f), (%f, %f)\n",
//                    faceUv[0].x(), faceUv[0].z(),
//                    faceUv[1].x(), faceUv[1].z(),
//                    faceUv[2].x(), faceUv[2].z());
//             printf("BBox (%f, %f) ~ (%f, %f)\n",
//                    uvMin.x(), uvMin.z(),
//                    uvMax.x(), uvMax.z());
//         }
//     }
// }

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

    glBindAttribLocation(prog_id, 0, "uv");
    glBindAttribLocation(prog_id, 1, "position");
    glBindFragDataLocation(prog_id, 0, "fragment");
    glUseProgram(prog_id);

    glfwSwapInterval(1);
    glClearColor(0.f, 0.f, 0.f, 1.f);

    std::unique_ptr<const Shape> shape(new Shape(2, 3, rectangleVertex));

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
