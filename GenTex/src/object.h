#pragma once
#include <GL/glew.h>

class Object {
    GLuint vao;
    GLuint vbo;

public:
    struct Vertex
    {
        GLfloat position[2];
    };
    // コンストラクタ
    // size: 頂点の位置の次元
    // vertexcount: 頂点の数
    // vertex: 頂点属性を格納した配列
    Object (GLint size, GLsizei vertexcount, const Vertex *vertex) {
        // 頂点配列オブジェクト
        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);
        // 頂点バッファオブジェクト
        glGenBuffers(1, &vbo);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER,
                     vertexcount * sizeof (Vertex), vertex, GL_STATIC_DRAW);
        // 結合されている頂点バッファオブジェクトを in 変数から参照できるようにする
        glVertexAttribPointer(0, size, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(0);
    }

    virtual ~Object() {
        // 頂点配列オブジェクトを削除する
        glDeleteBuffers(1, &vao);
        // 頂点バッファオブジェクトを削除する
        glDeleteBuffers(1, &vbo);
    }

    void bind() const {
        // 描画する頂点配列オブジェクトを指定する
        glBindVertexArray(vao);
    }

private:
    // forbid to copy
    Object(const Object &o);
    Object &operator=(const Object &o);
};
