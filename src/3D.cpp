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
    int size = 100;
    int _x = size;
    int _y = size;
    int _z = size;
    int px1 = _x / 2;
    int py1 = _y / 2;
    int pz = _z / 2;
    std::string path = {"/home/alex/src/diplom_FDTD/"};

    Mesh3D mesh{_x, _y, _z};
    // Obj obj(75, 75, 1, 50, 50, 18);
    Obj obj(0, 0, 0, 0, 0, 0);
    FDTD grid(mesh, obj);

    float j = 1;

    int tick = 250;
    int step = 1;
    std::string name = ".txt";
    for (int i = 0; i < tick; i++)
    {
        for (int z = 0; z < _z; z++)
        {
            grid.J.setNode(px1, py1, z, X, j);
            grid.J.setNode(px1, py1, z, Y, j);
            grid.J.setNode(px1, py1, z, Z, j);
        }
        grid.update(dt);
    }
    grid.H.saveExToFileZ(path + "mur_Ex" + std::to_string(0) + name, X);
    grid.H.saveExToFileZ(path + "mur_Ey" + std::to_string(0) + name, Y);
    for (int i = 1; i < 20; i++)
    {
        for (int z = 0; z < _z; z++)
        {
            grid.J.setNode(px1, py1, z, X, j);
            grid.J.setNode(px1, py1, z, Y, j);
            grid.J.setNode(px1, py1, z, Z, j);
        }
        grid.update(dt * step);
        grid.H.saveExToFileZ(path + "mur_Ex" + std::to_string(i) + name, X);
        grid.H.saveExToFileZ(path + "mur_Ey" + std::to_string(i) + name, Y);
    }
    std::cout << "\a";
    return 0;
}
