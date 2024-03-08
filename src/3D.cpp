#include <SFML/Graphics.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <chrono>
#include "FDTD.h"
#include "Obj.h"

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
    Obj obj(75, 75, 1, 50, 50, 18);
    FDTD grid(mesh, obj);

    float j = -1;
    px -= _x / 4;
    py -= _y / 4;

    for (int z = 0; z < _z; z++)
    {
        grid.J.setNode(px, py, z, X, j);
        grid.J.setNode(px, py, z, Y, j);
        grid.J.setNode(px, py, z, Z, j);
    }
    j = 1;
    px += _x / 2;
    py += _y / 2;
    for (int z = 0; z < _z; z++)
    {
        grid.J.setNode(px, py, z, X, j);
        grid.J.setNode(px, py, z, Y, j);
        grid.J.setNode(px, py, z, Z, j);
    }

    int tick = 150;
    grid.update(dt * tick);
    grid.H.saveExToFileZ(path + "Ex.txt", X);
    grid.H.saveExToFileZ(path + "Ey.txt", Y);

    return 0;
}
