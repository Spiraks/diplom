
#include <vector>
#include "FDTD.h"

#define murE(x, y, z, x1, y1, z1, cf, axis) ({ \
E.setNode(x, y, z, axis,    \
(E.gtO(x1, y1, z1, axis) + cf * (E.gtO(x, y, z, axis) - \
E.gt(x1, y1, z1, axis)))); })

#define murH(x, y, z, x1, y1, z1, cf, axis) ({ \
H.setNode(x, y, z, axis,    \
(H.gtO(x1, y1, z1, axis) + cf * (H.gtO(x, y, z, axis) - \
H.gt(x1, y1, z1, axis)))); })

void FDTD::MursFirstABC()
{
    // std::cout << "MursFirstABC\n";
    size_t s = 1;
    size_t s1 = s + 1;
    size_t end = 2;
    size_t x = len_x - end;
    size_t y = len_y - end;
    size_t z = len_z - end;
    size_t x1 = len_x - end - 1;
    size_t y1 = len_y - end - 1;
    size_t z1 = len_z - end - 1;

    // УСЛОВИЯ ПОГЛОЩЕНИЯ. Массивы вида EZX1, EZXN и т.п. - массивы сохраненных на предыдущем шаге по времени значений поля Е в приграничной области. Процедура сохранения приводится ниже.
    // коэффициенты Mur_Tx, Mur_Ty и Mur_Tz равны c*dT/dX, c*dT/dY и c*dT/dZ соответственно. c-скорость света, dT- шаг по времени; dX, dY, dZ - шаги по пространству.

    // ЛЕВАЯ И ПРАВАЯ ПЛОСКОСТИ
    //  Ez для двух плоскостей YxZ (левой и правой)
    // Ey для двух плоскостей YxZ [левой и правой]
    float _cfMur = (Mur_Tx - 1) / (Mur_Tx + 1);
    for (int J = s; J <= y; J++)
    {
        for (int K = s; K <= z; K++)
        {
            // std::cout << K << "-" << J << "\n";
            murE(s, J, K, s1, J, K, _cfMur, X);
            murE(s, J, K, s1, J, K, _cfMur, Y);
            murE(s, J, K, s1, J, K, _cfMur, Z);

            murE(x, J, K, x1, J, K, _cfMur, X);
            murE(x, J, K, x1, J, K, _cfMur, Y);
            murE(x, J, K, x1, J, K, _cfMur, Z);



            murH(s, J, K, s1, J, K, _cfMur, X);
            murH(s, J, K, s1, J, K, _cfMur, Y);
            murH(s, J, K, s1, J, K, _cfMur, Z);

            murH(x, J, K, x1, J, K, _cfMur, X);
            murH(x, J, K, x1, J, K, _cfMur, Y);
            murH(x, J, K, x1, J, K, _cfMur, Z);

        }
    }
    // ВЕРХНЯЯ И НИЖНЯЯ ПЛОСКОСТИ
    //  Ex для двух плоскостей XxZ [верхней и нижней]
    // Ez для двух плоскостей XxZ [верхней и нижней]

    for (int K = s; K <= z; K++)
    {
        for (int I = s; I <= x; I++)
        {
            murE(I, s, K, I, s1, K, _cfMur, X);
            murE(I, s, K, I, s1, K, _cfMur, Y);
            murE(I, s, K, I, s1, K, _cfMur, Z);

            murE(I, y, K, I, y1, K, _cfMur, X);
            murE(I, y, K, I, y1, K, _cfMur, Y);
            murE(I, y, K, I, y1, K, _cfMur, Z);



            murH(I, s, K, I, s1, K, _cfMur, X);
            murH(I, s, K, I, s1, K, _cfMur, Y);
            murH(I, s, K, I, s1, K, _cfMur, Z);

            murH(I, y, K, I, y1, K, _cfMur, X);
            murH(I, y, K, I, y1, K, _cfMur, Y);
            murH(I, y, K, I, y1, K, _cfMur, Z);
        }
    }

    // БЛИЖНЯЯ И ДАЛЬНЯЯ ПЛОСКОСТИ
    // Ey для двух плоскостей XxY [ближней и дальней]
    // Ex для двух плоскостей XxY [ближней и дальней]
    // std::cout << "5 Ey для двух плоскостей XxY [ближней и дальней]\n";
    for (int I = s; I <= x; I++)
    {
        for (int J = s; J <= y; J++)
        {
            murE(I, J, s, I, J, s1, _cfMur, X);
            murE(I, J, s, I, J, s1, _cfMur, Y);
            murE(I, J, s, I, J, s1, _cfMur, Z);

            murE(I, J, z, I, J, z1, _cfMur, X);
            murE(I, J, z, I, J, z1, _cfMur, Y);
            murE(I, J, z, I, J, z1, _cfMur, Z);



            murH(I, J, s, I, J, s1, _cfMur, X);
            murH(I, J, s, I, J, s1, _cfMur, Y);
            murH(I, J, s, I, J, s1, _cfMur, Z);

            murH(I, J, z, I, J, z1, _cfMur, X);
            murH(I, J, z, I, J, z1, _cfMur, Y);
            murH(I, J, z, I, J, z1, _cfMur, Z);
        }
    }

    // РЕБРА, параллельные оси Z
    // std::cout << "7 РЕБРА, параллельные оси Z\n";

        for (int K = s; K <= z; K++)
        {
            murE(s, s, K, s1, s1, K, _cfMur, X);
            murE(s, y, K, s1, y1, K, _cfMur, X);
            murE(x, s, K, x1, s1, K, _cfMur, X);
            murE(x, y, K, x1, y1, K, _cfMur, X);

            murE(s, s, K, s1, s1, K, _cfMur, Y);
            murE(s, y, K, s1, y1, K, _cfMur, Y);
            murE(x, s, K, x1, s1, K, _cfMur, Y);
            murE(x, y, K, x1, y1, K, _cfMur, Y);

            murE(s, s, K, s1, s1, K, _cfMur, Z);
            murE(s, y, K, s1, y1, K, _cfMur, Z);
            murE(x, s, K, x1, s1, K, _cfMur, Z);
            murE(x, y, K, x1, y1, K, _cfMur, Z);



            murH(s, s, K, s1, s1, K, _cfMur, X);
            murH(s, y, K, s1, y1, K, _cfMur, X);
            murH(x, s, K, x1, s1, K, _cfMur, X);
            murH(x, y, K, x1, y1, K, _cfMur, X);

            murH(s, s, K, s1, s1, K, _cfMur, Y);
            murH(s, y, K, s1, y1, K, _cfMur, Y);
            murH(x, s, K, x1, s1, K, _cfMur, Y);
            murH(x, y, K, x1, y1, K, _cfMur, Y);

            murH(s, s, K, s1, s1, K, _cfMur, Z);
            murH(s, y, K, s1, y1, K, _cfMur, Z);
            murH(x, s, K, x1, s1, K, _cfMur, Z);
            murH(x, y, K, x1, y1, K, _cfMur, Z);
        }

    // РЕБРА, параллельные оси Y
    // std::cout << "8 РЕБРА, параллельные оси Y\n";

        size_t cf = Y;
        for (int J = s; J <= y; J++)
        {
            murE(s, J, s, s1, J, s1, _cfMur, X);
            murE(s, J, z, s1, J, z1, _cfMur, X);
            murE(x, J, s, x1, J, s1, _cfMur, X);
            murE(x, J, z, x1, J, z1, _cfMur, X);

            murE(s, J, s, s1, J, s1, _cfMur, Y);
            murE(s, J, z, s1, J, z1, _cfMur, Y);
            murE(x, J, s, x1, J, s1, _cfMur, Y);
            murE(x, J, z, x1, J, z1, _cfMur, Y);

            murE(s, J, s, s1, J, s1, _cfMur, Z);
            murE(s, J, z, s1, J, z1, _cfMur, Z);
            murE(x, J, s, x1, J, s1, _cfMur, Z);
            murE(x, J, z, x1, J, z1, _cfMur, Z);



            murH(s, J, s, s1, J, s1, _cfMur, X);
            murH(s, J, z, s1, J, z1, _cfMur, X);
            murH(x, J, s, x1, J, s1, _cfMur, X);
            murH(x, J, z, x1, J, z1, _cfMur, X);

            murH(s, J, s, s1, J, s1, _cfMur, Y);
            murH(s, J, z, s1, J, z1, _cfMur, Y);
            murH(x, J, s, x1, J, s1, _cfMur, Y);
            murH(x, J, z, x1, J, z1, _cfMur, Y);

            murH(s, J, s, s1, J, s1, _cfMur, Z);
            murH(s, J, z, s1, J, z1, _cfMur, Z);
            murH(x, J, s, x1, J, s1, _cfMur, Z);
            murH(x, J, z, x1, J, z1, _cfMur, Z);
        }

    // РЕБРА, параллельные оси X
    // std::cout << "9 РЕБРА, параллельные оси X\n";
        for (int I = s; I <= x; I++)
        {
            // std::cout << I << "\n";
            murE(I, s, s, I, s1, s1, _cfMur, X);
            murE(I, s, z, I, s1, z1, _cfMur, X);
            murE(I, y, s, I, y1, s1, _cfMur, X);
            murE(I, y, z, I, y1, z1, _cfMur, X);

            murE(I, s, s, I, s1, s1, _cfMur, Y);
            murE(I, s, z, I, s1, z1, _cfMur, Y);
            murE(I, y, s, I, y1, s1, _cfMur, Y);
            murE(I, y, z, I, y1, z1, _cfMur, Y);

            murE(I, s, s, I, s1, s1, _cfMur, Z);
            murE(I, s, z, I, s1, z1, _cfMur, Z);
            murE(I, y, s, I, y1, s1, _cfMur, Z);
            murE(I, y, z, I, y1, z1, _cfMur, Z);



            murH(I, s, s, I, s1, s1, _cfMur, X);
            murH(I, s, z, I, s1, z1, _cfMur, X);
            murH(I, y, s, I, y1, s1, _cfMur, X);
            murH(I, y, z, I, y1, z1, _cfMur, X);

            murH(I, s, s, I, s1, s1, _cfMur, Y);
            murH(I, s, z, I, s1, z1, _cfMur, Y);
            murH(I, y, s, I, y1, s1, _cfMur, Y);
            murH(I, y, z, I, y1, z1, _cfMur, Y);
            
            murH(I, s, s, I, s1, s1, _cfMur, Z);
            murH(I, s, z, I, s1, z1, _cfMur, Z);
            murH(I, y, s, I, y1, s1, _cfMur, Z);
            murH(I, y, z, I, y1, z1, _cfMur, Z);
        }
}

void FDTD::GetBorderValues()
{
    E.copyNodes();
    H.copyNodes();
}
