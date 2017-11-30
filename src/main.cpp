#include <iostream>

#define TINYOBJLOADER_IMPLEMENTATION
#include "./tiny_obj_loader.h"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cout << "Needs input.obj\n" << std::endl;
        return 0;
    }

    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;

    std::string err;

    std::cout << "Load " << argv[1] << std::endl;
    bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, argv[1]);

    if (!err.empty()) { // `err` may contain warning message.
        std::cerr << err << std::endl;
    }

    if (!ret) {
        std::cerr << "Failed to load " << argv[1] << std::endl;
        exit(1);
    }

    printf("# of vertices  = %d\n", (int) (attrib.vertices.size()) / 3);
    printf("# of normals   = %d\n", (int) (attrib.normals.size()) / 3);
    printf("# of texcoords = %d\n", (int) (attrib.texcoords.size()) / 2);
    printf("# of materials = %d\n", (int) materials.size());
    printf("# of shapes    = %d\n", (int) shapes.size());
    printf("\n");

    for (size_t s = 0; s < shapes.size(); s++) {
        printf("Shape %d\n", (int) shapes.size());
        printf("# of faces %d\n", (int) shapes[s].mesh.num_face_vertices.size());
        size_t index_offset = 0;
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
            // The number of vertexes in the face. 3, 4 or 5?
            int fv = shapes[s].mesh.num_face_vertices[f];

            // Loop over vertices in the face.
            for (int v = 0; v < fv; v++) {
                // access to vertex
                tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
                tinyobj::real_t vx = attrib.vertices[3*idx.vertex_index+0];
                tinyobj::real_t vy = attrib.vertices[3*idx.vertex_index+1];
                tinyobj::real_t vz = attrib.vertices[3*idx.vertex_index+2];
                tinyobj::real_t nx = attrib.normals[3*idx.normal_index+0];
                tinyobj::real_t ny = attrib.normals[3*idx.normal_index+1];
                tinyobj::real_t nz = attrib.normals[3*idx.normal_index+2];
                tinyobj::real_t tx = attrib.texcoords[2*idx.texcoord_index+0];
                tinyobj::real_t ty = attrib.texcoords[2*idx.texcoord_index+1];
                // Optional: vertex colors
                // tinyobj::real_t red = attrib.colors[3*idx.vertex_index+0];
                // tinyobj::real_t green = attrib.colors[3*idx.vertex_index+1];
                // tinyobj::real_t blue = attrib.colors[3*idx.vertex_index+2];
            }
            index_offset += fv;

            // per-face material
            shapes[s].mesh.material_ids[f];
        }
    }


    return 0;
}
