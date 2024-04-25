// #include <SFML/Graphics.hpp>

// #include <iostream>
// #include <fstream>
// #include <unistd.h>
// #include <chrono>
// #include "FDTD.h"

// int main()
// {

//     int _x = 200;
//     int _y = 200;
//     int _z = 30;
//     int px = _x / 2;
//     int py = _y / 2;
//     int pz = _z / 2;

//     sf::RenderWindow window(sf::VideoMode(_x, _y), "3D Wave Simulation");
//     Mesh3D mesh{_x, _y, _z};
//     FDTD grid(mesh);
//     My_Field E;
//     E = My_Field(mesh);
//     float j = 1;
//     for (int z = 0; z < _z; z++)
//     {
//         grid.J.setNode(px, py, z, X, j);
//         grid.J.setNode(px, py, z, Y, j);
//         grid.J.setNode(px, py, z, Z, j);
//     }

//     while (window.isOpen())
//     {
//         sf::Event event;
//         while (window.pollEvent(event))
//         {
//             if (event.type == sf::Event::Closed)
//             {
//                 window.close();
//             }
//         }
//         grid.update(dt * 50);

//         // Draw the wave
//         sf::Texture texture;
//         texture.create(_x, _y);
//         sf::Uint8 *pixels = new sf::Uint8[_x * _y * 4];
//         for (int x = 1; x < _x; ++x)
//         {
//             for (int y = 1; y < _y; y++)
//             {
//                 int index = (x * _y + y) * 4;

//                 float vx = abs(grid.H.getNodeE(x, y, X));
//                 float vy = abs(grid.H.getNodeE(x, y, Y));
//                 if (vx < 10 && vx > 0)
//                 {
//                     pixels[index] = (2 + vx) * 20;
//                     pixels[index + 1] = 0;
//                     pixels[index + 2] = (2 + vy) * 20;
//                     pixels[index + 3] = 255;
//                     continue;
//                 }
//                 // float v = (vx+vy)/2;
//                 pixels[index] = static_cast<sf::Uint8>(vx);
//                 pixels[index + 1] = 0;
//                 pixels[index + 2] = static_cast<sf::Uint8>(vy);
//                 pixels[index + 3] = 255;
//             }
//         }

//         texture.update(pixels);

//         grid.H.saveExToFileZ("/home/alex/src/diplom_FDTD/Ex.txt", X);
//         grid.H.saveExToFileZ("/home/alex/src/diplom_FDTD/Ey.txt", Y);
//         return 0;
//         delete[] pixels;
//         sf::Sprite sprite(texture);
//         window.clear();
//         window.draw(sprite);
//         window.display();
//     }

//     return 0;
// }
