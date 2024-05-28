#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <stdlib.h>
#include <signal.h>
#include <chrono>
#include "FDTD.h"

#define setPoint(x,y,z)({ \
            grid.J._X[x][y][z]= j;\
            grid.J._Y[x][y][z]= j;\
            grid.J._Z[x][y][z]= j;})
#define setBase(x,y,z) ({\
            grid.Sigm._X[x][y][z]= 0;\
            grid.Sigm._Y[x][y][z]= 0;\
            grid.Sigm._Z[x][y][z]= 0;\
            grid.Mu._X[x][y][z]= 1;\
            grid.Mu._Y[x][y][z]= 1;\
            grid.Mu._Z[x][y][z]= 1;\
            grid.Eps._X[x][y][z]= 4;\
            grid.Eps._Y[x][y][z]= 4;\
            grid.Eps._Z[x][y][z]= 4;})
#define setPlot(x,y,z) ({\
            grid.Sigm._X[x][y][z]= 1000000;\
            grid.Sigm._Y[x][y][z]= 1000000;\
            grid.Sigm._Z[x][y][z]= 1000000;\
            grid.Mu._X[x][y][z]= 1;\
            grid.Mu._Y[x][y][z]= 1;\
            grid.Mu._Z[x][y][z]= 1;\
            grid.Eps._X[x][y][z]= 1;\
            grid.Eps._Y[x][y][z]= 1;\
            grid.Eps._Z[x][y][z]= 1;})

int main(void)
{
    int size = 150;
    int _x = 160;
    int _y = 160;
    int _z = 300;
    int px1 = _x / 2;
    int py1 = _y / 2;
    std::string path = {"/home/alex/src/diplom_FDTD/results/"};


    Mesh3D mesh{_x, _y, _z};
    FDTD grid(mesh);

    float j = -1;
    int tick = 300; //60
    int step = 5;
    std::string name = ".txt";
    auto start = std::chrono::system_clock::now();

    for (size_t i = 0; i < tick; i++)
    {

        for(int a = 0; a < 18; a++){
            for(int b = 0; b < 18; b++){
                setPoint(px1-9+a,py1-9+b,50);
                setPoint(px1-9+a,py1-9+b,51);
            }
        }
        
        grid.update(dt);
    }
    grid.H.saveExToFileZX(path + "mur_Ezx" + name);
    grid.H.saveExToFileZY(path + "mur_Ezy" + name);
    grid.H.saveExToFileXY(path + "mur_Exy" + name);
    grid.H.saveExToFileXZ(path + "mur_Exz" + name);
    grid.H.saveExToFileYX(path + "mur_Eyx" + name);
    grid.H.saveExToFileYZ(path + "mur_Eyz" + name);

    // size_t num = 1;
    // for (size_t i = 0; i < 150; i++)
    // {
    //     for (size_t z = 0; z < _z; z++)
    //     {
    //         grid.J.setNode(px1, py1, z, X, 1);
    //         grid.J.setNode(px1, py1, z, Y, 1);
    //         grid.J.setNode(px1, py1, z, Z, 1);
    //     }
    //     grid.update(dt);
    //     if (i % step == 0)
    //     {
    //         grid.H.saveExToFileZ(path + "mur_Ezx" + std::to_string(num) + name, X);
    //         grid.H.saveExToFileZ(path + "mur_Ezy" + std::to_string(num) + name, Y);
    //         num++;
    //     }
    // }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s"
              << std::endl;
    // system("canberra-gtk-play -f /home/alex/src/diplom_FDTD/boom.mp3");
    return 0;
}
