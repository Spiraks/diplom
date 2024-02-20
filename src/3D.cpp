#include <SFML/Graphics.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <chrono>
#include "FDTD.h"

int main()
{

    int _x = 200;
    int _y = 200;
    int _z = 20;
    int px = _x / 2;
    int py = _y / 2;
    int pz = _z / 2;
    std::string path = {"/home/alex/src/diplom_FDTD/"};

    Mesh3D mesh{_x, _y, _z};
    FDTD grid(mesh);
    Field E;
    E = Field(mesh);
    float j = -1;
    for (int z = 0; z < _z; z++)
    {
        grid.J.setNode(px, py, z, X, j);
        grid.J.setNode(px, py, z, Y, j);
        grid.J.setNode(px, py, z, Z, j);
    }
    int tick = 50;
    grid.update(dt * tick);
    // for (int i = 0; i < 5; ++i)
    // {
    //     grid.update(dt * (tick += 10));
    grid.H.saveExToFileZ(path + "Ex.txt", X);
    grid.H.saveExToFileZ(path + "Ey.txt", Y);
    // }

    return 0;
}
