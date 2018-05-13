#pragma once
#include <GL/glew.h>

class Object {
    GLuint vao;
    GLuint vbo;

public:
    struct Vertex
    {
        GLfloat position[2];
        GLfloat color[3];
    };
    // コンストラクタ
    // size: 頂点の位置の次元
    // vertexcount: 頂点の数
    // vertex: 頂点属性を格納した配列
    Object (GLint size, GLsizei vertexcount, const Vertex *vertex) {

        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);

        glGenBuffers(1, &vbo);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER,
                     vertexcount * sizeof (Vertex), vertex, GL_STATIC_DRAW);

        glVertexAttribPointer(0, size, GL_FLOAT, GL_FALSE, sizeof (Vertex),
                              static_cast<Vertex *>(0)->position);
        glEnableVertexAttribArray(0);

        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof (Vertex),
                              static_cast<Vertex *>(0)->color);
        glEnableVertexAttribArray(1);
    }

    virtual ~Object() {
        glDeleteBuffers(1, &vao);
        glDeleteBuffers(1, &vbo);
    }

    void bind() const {
        glBindVertexArray(vao);
    }

private:
    // forbid to copy
    Object(const Object &o);
    Object &operator=(const Object &o);
};
