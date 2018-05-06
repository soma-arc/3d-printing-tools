#pragma once
#include <memory>
#include "object.h"

class Shape {
    std::shared_ptr<const Object> object;

protected:
    const GLsizei vertexcount;

public:
    // コンストラクタ
    // size: 頂点の位置の次元
    // vertexcount: 頂点の数
    // vertex: 頂点属性を格納した配列
    Shape(GLint size, GLsizei vertexcount, const Object::Vertex *vertex)
        : object(new Object(size, vertexcount, vertex))
        , vertexcount(vertexcount)
        {}

    void draw() const {
        object->bind();
        execute();
    }

    virtual void execute() const {
        glDrawArrays(GL_TRIANGLES, 0, vertexcount);
    }
};
