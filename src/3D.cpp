#include <SFML/Graphics.hpp>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <stdlib.h>
#include <signal.h>
#include <chrono>
#include "FDTD.h"
#include "Obj.h"

int main(void)
{
    int size = 100;
    int _x = size;
    int _y = size;
    int _z = size;
    int px1 = _x / 2;
    int py1 = _y / 2;
    std::string path = {"/home/alex/src/diplom_FDTD/"};

    Mesh3D mesh{_x, _y, _z};
    // Obj obj(75, 75, 1, 50, 50, 18);
    Obj obj(0, 0, 0, 0, 0, 0);
    FDTD grid(mesh, obj);

    float j = 1;

    int tick = 170;
    int step = 5;
    std::string name = ".txt";
    for (int i = 0; i < tick; i++)
    {
        for (int z = 0; z < _z; z++)
        {
            grid.J.setNode(px1, py1, z, X, 1);
            grid.J.setNode(px1, py1, z, Y, 1);
            grid.J.setNode(px1, py1, z, Z, 1);
        }
        grid.update(dt);
    }
    grid.H.saveExToFileZ(path + "mur_Ex" + std::to_string(0) + name, X);
    grid.H.saveExToFileZ(path + "mur_Ey" + std::to_string(0) + name, Y);
    // size_t num = 1;
    // for (int i = 0; i < 40; i++)
    // {
    //     for (int z = 0; z < _z; z++)
    //     {
    //         grid.J.setNode(px1, py1, z, X, 1);
    //         grid.J.setNode(px1, py1, z, Y, 1);
    //         grid.J.setNode(px1, py1, z, Z, 1);
    //     }
    //     grid.update(dt);
    //     if (i % step == 0)
    //     {

    //         grid.H.saveExToFileZ(path + "mur_Ex" + std::to_string(num) + name, X);
    //         grid.H.saveExToFileZ(path + "mur_Ey" + std::to_string(num) + name, Y);
    //         num++;
    //     }
    // }
    system("canberra-gtk-play -f /home/alex/src/diplom_FDTD/gtr-nylon22.wav");
    return 0;
}
