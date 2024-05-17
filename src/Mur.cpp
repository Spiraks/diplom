
#include <vector>
#include "FDTD.h"

#define _murX(x, y, z, x1, y1, z1, cf,fiels) ({ \
fiels._X[x][y][z] = (fiels._X_1[x1][y1][z1] + cf * (fiels._X_1[x][y][z] - fiels._X[x1][y1][z1])); })

#define _murY(x, y, z, x1, y1, z1, cf,fiels) ({ \
fiels._Y[x][y][z] = (fiels._Y_1[x1][y1][z1] + cf * (fiels._Y_1[x][y][z] - fiels._Y[x1][y1][z1])); })

#define _murZ(x, y, z, x1, y1, z1, cf,fiels) ({ \
fiels._Z[x][y][z] = (fiels._Z_1[x1][y1][z1] + cf * (fiels._Z_1[x][y][z] - fiels._Z[x1][y1][z1])); })



#define mur(x, y, z, x1, y1, z1, cf) ({ \
_murX(x, y, z, x1, y1, z1, cf, H);\
_murY(x, y, z, x1, y1, z1, cf, H);\
_murZ(x, y, z, x1, y1, z1, cf, H);\
_murX(x, y, z, x1, y1, z1, cf, E);\
_murY(x, y, z, x1, y1, z1, cf, E);\
_murZ(x, y, z, x1, y1, z1, cf, E); })

void FDTD::Mur()
{
    // std::cout << "Mur\n";
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

    for (size_t J = s; J <= y; J++)
    {
        for (size_t K = s; K <= z; K++)
        {
            // std::cout << K << "-" << J << "\n";
            mur(s, J, K, s1, J, K, _cfMur);
            mur(x, J, K, x1, J, K, _cfMur);
        }
    }
    // ВЕРХНЯЯ И НИЖНЯЯ ПЛОСКОСТИ
    //  Ex для двух плоскостей XxZ [верхней и нижней]
    // Ez для двух плоскостей XxZ [верхней и нижней]

    for (size_t K = s; K <= z; K++)
    {
        for (size_t I = s; I <= x; I++)
        {
            mur(I, s, K, I, s1, K, _cfMur);
            mur(I, y, K, I, y1, K, _cfMur);
        }
    }

    // БЛИЖНЯЯ И ДАЛЬНЯЯ ПЛОСКОСТИ
    // Ey для двух плоскостей XxY [ближней и дальней]
    // Ex для двух плоскостей XxY [ближней и дальней]
    // std::cout << "5 Ey для двух плоскостей XxY [ближней и дальней]\n";
    for (size_t I = s; I <= x; I++)
    {
        for (size_t J = s; J <= y; J++)
        {
            mur(I, J, s, I, J, s1, _cfMur);
            mur(I, J, z, I, J, z1, _cfMur);
        }
    }

    // РЕБРА, параллельные оси Z
    // std::cout << "7 РЕБРА, параллельные оси Z\n";

    for (size_t K = s; K <= z; K++)
    {
        mur(s, s, K, s1, s1, K, _cfMur);
        mur(s, y, K, s1, y1, K, _cfMur);
        mur(x, s, K, x1, s1, K, _cfMur);
        mur(x, y, K, x1, y1, K, _cfMur);
    }

    // РЕБРА, параллельные оси Y
    // std::cout << "8 РЕБРА, параллельные оси Y\n";

    for (size_t J = s; J <= y; J++)
    {
        mur(s, J, s, s1, J, s1, _cfMur);
        mur(s, J, z, s1, J, z1, _cfMur);
        mur(x, J, s, x1, J, s1, _cfMur);
        mur(x, J, z, x1, J, z1, _cfMur);

    }

    // РЕБРА, параллельные оси X
    // std::cout << "9 РЕБРА, параллельные оси X\n";
    for (size_t I = s; I <= x; I++)
    {
        mur(I, s, s, I, s1, s1, _cfMur);
        mur(I, s, z, I, s1, z1, _cfMur);
        mur(I, y, s, I, y1, s1, _cfMur);
        mur(I, y, z, I, y1, z1, _cfMur);
    }
}

void FDTD::GetBorderValuesMur()
{
    E.copyNodes();
    H.copyNodes();
}
