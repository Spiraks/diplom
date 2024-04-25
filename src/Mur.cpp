
#include <vector>
#include "FDTD.h"

void FDTD::MursFirstABC()
{
    // std::cout << "MursFirstABC\n";
    int start = 2;

    // УСЛОВИЯ ПОГЛОЩЕНИЯ. Массивы вида EZX1, EZXN и т.п. - массивы сохраненных на предыдущем шаге по времени значений поля Е в приграничной области. Процедура сохранения приводится ниже.
    // коэффициенты Mur_Tx, Mur_Ty и Mur_Tz равны c*dT/dX, c*dT/dY и c*dT/dZ соответственно. c-скорость света, dT- шаг по времени; dX, dY, dZ - шаги по пространству.

    // ЛЕВАЯ И ПРАВАЯ ПЛОСКОСТИ
    //  Ez для двух плоскостей YxZ (левой и правой)
    // Ey для двух плоскостей YxZ [левой и правой]
    float val = 0;
    float _coef_mur = (Mur_Tx - 1) / (Mur_Tx + 1);
    for (int J = 2; J < len_y - 2; J++)
    {
        for (int K = 2; K < len_z - 2; K++)
        {
            // std::cout << K << "-" << J << "\n";

            val = E.getNodeOld(start + 1, J, K, Z) + _coef_mur * (E.getNodeOld(start, J, K, Z) - E.getNode(start + 1, J, K, Z));
            E.setNode(start, J, K, Z, val);

            val = E.getNodeOld(start + 1, J, K, Y) + _coef_mur * (E.getNodeOld(start, J, K, Y) - E.getNode(start + 1, J, K, Y));
            E.setNode(start, J, K, Y, val);

            val = E.getNodeOld(len_x - 4, J, K, Z) + _coef_mur * (E.getNodeOld(len_x - 3, J, K, Z) - E.getNode(len_x - 4, J, K, Z));
            E.setNode(len_x - 3, J, K, Z, val);

            val = E.getNodeOld(len_x - 4, J, K, Y) + _coef_mur * (E.getNodeOld(len_x - 3, J, K, Y) - E.getNode(len_x - 4, J, K, Y));
            E.setNode(len_x - 3, J, K, Y, val);
        }
    }

    // ВЕРХНЯЯ И НИЖНЯЯ ПЛОСКОСТИ
    //  Ex для двух плоскостей XxZ [верхней и нижней]
    // Ez для двух плоскостей XxZ [верхней и нижней]

    _coef_mur = (Mur_Ty - 1) / (Mur_Ty + 1);
    for (int K = 2; K < len_z - 2; K++)
    {
        for (int I = 2; I < len_x - 2; I++)
        {
            val = E.getNodeOld(I, start + 1, K, X) + _coef_mur * (E.getNodeOld(I, start, K, X) - E.getNode(I, start + 1, K, X));
            E.setNode(I, start, K, X, val);

            val = E.getNodeOld(I, start + 1, K, Z) + _coef_mur * (E.getNodeOld(I, start, K, Z) - E.getNode(I, start + 1, K, Z));
            E.setNode(I, start, K, Z, val);

            val = E.getNodeOld(I, len_y - 4, K, X) + _coef_mur * (E.getNodeOld(I, len_y - 3, K, X) - E.getNode(I, len_y - 4, K, X));
            E.setNode(I, len_y - 3, K, X, val);

            val = E.getNodeOld(I, len_y - 4, K, Z) + _coef_mur * (E.getNodeOld(I, len_y - 3, K, Z) - E.getNode(I, len_y - 4, K, Z));
            E.setNode(I, len_y - 3, K, Z, val);
        }
    }

    // БЛИЖНЯЯ И ДАЛЬНЯЯ ПЛОСКОСТИ
    // Ey для двух плоскостей XxY [ближней и дальней]
    // Ex для двух плоскостей XxY [ближней и дальней]
    // std::cout << "5 Ey для двух плоскостей XxY [ближней и дальней]\n";
    _coef_mur = (Mur_Tz - 1) / (Mur_Tz + 1);
    for (int I = 2; I < len_x - 2; I++)
    {
        for (int J = 2; J < len_y - 2; J++)
        {
            val = E.getNodeOld(I, J, start + 1, Y) + _coef_mur * (E.getNodeOld(I, J, start, Y) - E.getNode(I, J, start + 1, Y));
            E.setNode(I, J, start, Y, val);

            val = E.getNodeOld(I, J, start + 1, X) + _coef_mur * (E.getNodeOld(I, J, start, X) - E.getNode(I, J, start + 1, X));
            E.setNode(I, J, 1, X, val);

            val = E.getNodeOld(I, J, len_z - 4, Y) + _coef_mur * (E.getNodeOld(I, J, len_z - 3, Y) - E.getNode(I, J, len_z - 4, Y));
            E.setNode(I, J, len_z - 3, Y, val);

            val = E.getNodeOld(I, J, len_z - 4, X) + _coef_mur * (E.getNodeOld(I, J, len_z - 3, X) - E.getNode(I, J, len_z - 4, X));
            E.setNode(I, J, len_z - 3, X, val);
        }
    }

    // РЕБРА, параллельные оси Z
    // std::cout << "7 РЕБРА, параллельные оси Z\n";

    for (int K = 2; K < len_z - 2; K++)
    {
        val = E.getNodeOld(start + 1, start + 1, K, Z) + _coef_mur * (E.getNodeOld(start, start, K, Z) - E.getNode(start + 1, start + 1, K, Z));
        E.setNode(start, start, K, Z, val);

        val = E.getNodeOld(start + 1, len_y - 4, K, Z) + _coef_mur * (E.getNodeOld(start, len_y - 3, K, Z) - E.getNode(start + 1, len_y - 4, K, Z));
        E.setNode(start, len_y - 3, K, Z, val);

        val = E.getNodeOld(len_x - 4, start + 1, K, Z) + _coef_mur * (E.getNodeOld(len_x - 3, start, K, Z) - E.getNode(len_x - 4, start + 1, K, Z));
        E.setNode(len_x - 3, start, K, Z, val);

        val = E.getNodeOld(len_x - 4, len_y - 4, K, Z) + _coef_mur * (E.getNodeOld(len_x - 3, len_y - 3, K, Z) - E.getNode(len_x - 4, len_y - 4, K, Z));
        E.setNode(len_x - 3, len_y - 3, K, Z, val);
    }

    // РЕБРА, параллельные оси Y
    // std::cout << "8 РЕБРА, параллельные оси Y\n";

    for (int J = 2; J < len_y - 2; J++)
    {
        val = E.getNodeOld(start + 1, J, start + 1, Y) + _coef_mur * (E.getNodeOld(start, J, start, Y) - E.getNode(start + 1, J, start + 1, Y));
        E.setNode(start, J, start, Y, val);

        val = E.getNodeOld(start + 1, J, len_y - 4, Y) + _coef_mur * (E.getNodeOld(start, J, len_y - 3, Y) - E.getNode(start + 1, J, len_y - 4, Y));
        E.setNode(start, J, len_y - 3, Y, val);

        val = E.getNodeOld(len_x - 4, J, start + 1, Y) + _coef_mur * (E.getNodeOld(len_x - 3, J, start, Y) - E.getNode(len_x - 4, J, start + 1, Y));
        E.setNode(len_x - 3, J, start, Y, val);

        val = E.getNodeOld(len_x - 4, J, len_y - 4, Y) + _coef_mur * (E.getNodeOld(len_x - 3, J, len_y - 3, Y) - E.getNode(len_x - 4, J, len_y - 4, Y));
        E.setNode(len_x - 3, J, len_y - 3, Y, val);
    }

    // РЕБРА, параллельные оси X
    // std::cout << "9 РЕБРА, параллельные оси X\n";

    for (int I = 2; I < len_x - 2; I++)
    {
        // std::cout << I << "\n";

        val = E.getNodeOld(I, start + 1, start + 1, Y) + _coef_mur * (E.getNodeOld(I, start, start, Y) - E.getNode(I, start + 1, start + 1, Y));
        E.setNode(I, start, start, Y, val);

        val = E.getNodeOld(I, start + 1, len_y - 4, Y) + _coef_mur * (E.getNodeOld(I, start, len_y - 3, Y) - E.getNode(I, start + 1, len_y - 4, Y));
        E.setNode(I, start, len_y - 3, Y, val);

        val = E.getNodeOld(I, len_x - 4, start + 1, Y) + _coef_mur * (E.getNodeOld(I, len_x - 3, start, Y) - E.getNode(I, len_x - 4, start + 1, Y));
        E.setNode(I, len_x - 3, start, Y, val);

        val = E.getNodeOld(I, len_x - 4, len_y - 4, Y) + _coef_mur * (E.getNodeOld(I, len_x - 3, len_y - 3, Y) - E.getNode(I, len_x - 4, len_y - 4, Y));
        E.setNode(I, len_x - 3, len_y - 3, Y, val);
    }
}

void FDTD::GetBorderValues()
{
    E.copyNodes();
}
