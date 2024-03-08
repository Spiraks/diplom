
#include <vector>
#include "FDTD.h"

void MursFirstABC(Field &E)
{

    int NX = -1;
    int NX = -1;
    int NX = -1;
    vector EZN1;
    vector EZNN;

    std::vector<std::vector<std::vector<double>>>
        EZ(NX + 1, std::vector<std::vector<double>>(NY + 1, std::vector<double>(NZ + 1)));
    std::vector<std::vector<std::vector<double>>> EY(NX + 1, std::vector<std::vector<double>>(NY + 1, std::vector<double>(NZ + 1)));
    std::vector<std::vector<std::vector<double>>> EX(NX + 1, std::vector<std::vector<double>>(NY + 1, std::vector<double>(NZ + 1)));

    // УСЛОВИЯ ПОГЛОЩЕНИЯ. Массивы вида EZX1, EZXN и т.п. - массивы сохраненных на предыдущем шаге по времени значений поля Е в приграничной области. Процедура сохранения приводится ниже.
    // коэффициенты Mur_Tx, Mur_Ty и Mur_Tz равны c*dT/dX, c*dT/dY и c*dT/dZ соответственно. c-скорость света, dT- шаг по времени; dX, dY, dZ - шаги по пространству.

    // ЛЕВАЯ И ПРАВАЯ ПЛОСКОСТИ
    //  Ez для двух плоскостей YxZ (левой и правой)
    float val = 0;
    for (int J = 2; J <= NY - 1; J++)
    {
        for (int K = 1; K <= NZ - 1; K++)
        {
            val = EZX1[J][K] + (Mur_Tx - 1) / (Mur_Tx + 1) * (E.getNode(1, J, K, Z) - E.getNode(0, J, K, Z));
            E.setNode(0, J, K, Z, val);
            val = EZXN[J][K] + (Mur_Tx - 1) / (Mur_Tx + 1) * (E.getNode(NX - 2, J, K, Z) - E.getNode(NX - 1, J, K, Z));
            E.setNode(NX - 1, J, K, Z, val);
        }
    }

    // Ey для двух плоскостей YxZ [левой и правой]
    for (int J = 1; J <= NY - 1; J++)
    {
        for (int K = 2; K <= NZ - 1; K++)
        {
            EY[1][J][K] = EYX1[J][K] + (Mur_Tx - 1) / (Mur_Tx + 1) * (EY[2][J][K] - EY[1][J][K]);
            EY[NX][J][K] = EYXN[J][K] + (Mur_Tx - 1) / (Mur_Tx + 1) * (EY[NX - 1][J][K] - EY[NX][J][K]);
        }
    }

    // ВЕРХНЯЯ И НИЖНЯЯ ПЛОСКОСТИ
    //  Ex для двух плоскостей XxZ [верхней и нижней]
    for (K = 2; K <= NZ - 1; K++)
    {
        for (I = 1; I <= NX - 1; I++)
        {
            EX[I][1][K] = EXY1[I][K] + (Mur_Ty - 1) / (Mur_Ty + 1) * (EX[I][2][K] - EX[I][1][K]);
            EX[I][NY][K] = EXYN[I][K] + (Mur_Ty - 1) / (Mur_Ty + 1) * (EX[I][NY - 1][K] - EX[I][NY][K]);
        }
    }

    // Ez для двух плоскостей XxZ [верхней и нижней]
    for (K = 1; K <= NZ - 1; K++)
    {
        for (I = 2; I <= NX - 1; I++)
        {
            EZ[I][1][K] = EZY1[I][K] + (Mur_Ty - 1) / (Mur_Ty + 1) * (EZ[I][2][K] - EZ[I][1][K]);
            EZ[I][NY][K] = EZYN[I][K] + (Mur_Ty - 1) / (Mur_Ty + 1) * (EZ[I][NY - 1][K] - EZ[I][NY][K]);
        }
    }

    // БЛИЖНЯЯ И ДАЛЬНЯЯ ПЛОСКОСТИ
    // Ey для двух плоскостей XxY [ближней и дальней]
    for (I = 2; I <= NX - 1; I++)
    {
        for (J = 1; J <= NY - 1; J++)
        {
            EY[I][J][NZ] = EYZN[I][J] + (Mur_Tz - 1) / (Mur_Tz + 1) * (EY[I][J][NZ - 1] - EY[I][J][NZ]);
            EY[I][J][1] = EYZ1[I][J] + (Mur_Tz - 1) / (Mur_Tz + 1) * (EY[I][J][2] - EY[I][J][1]);
        }
    }

    // Ex для двух плоскостей XxY [ближней и дальней]
    for (I = 1; I <= NX - 1; I++)
    {
        for (J = 2; J <= NY - 1; J++)
        {
            EX[I][J][NZ] = EXZN[I][J] + (Mur_Tz - 1) / (Mur_Tz + 1) * (EX[I][J][NZ - 1] - EX[I][J][NZ]);
            EX[I][J][1] = EXZ1[I][J] + (Mur_Tz - 1) / (Mur_Tz + 1) * (EX[I][J][2] - EX[I][J][1]);
        }
    }

    // РЕБРА, параллельные оси Z
    for (K = 1; K <= NZ - 1; K++)
    {
        EZ[1][1][K] = EZX1[1][K] + (Mur_Tx - 1) / (Mur_Tx + 1) * (EZ[2][1][K] - EZ[1][1][K]);
        EZ[1][NY][K] = EZX1[NY][K] + (Mur_Tx - 1) / (Mur_Tx + 1) * (EZ[2][NY][K] - EZ[1][NY][K]);
        EZ[NX][1][K] = EZXN[1][K] + (Mur_Tx - 1) / (Mur_Tx + 1) * (EZ[NX - 1][1][K] - EZ[NX][1][K]);
        EZ[NX][NY][K] = EZXN[NY][K] + (Mur_Tx - 1) / (Mur_Tx + 1) * (EZ[NX - 1][NY][K] - EZ[NX][NY][K]);
    }

    // РЕБРА, параллельные оси Y
    for (J = 1; J <= NY - 1; J++)
    {
        EY[1][J][1] = EYX1[J][1] + (Mur_Tx - 1) / (Mur_Tx + 1) * (EY[2][J][1] - EY[1][J][1]);
        EY[1][J][NZ] = EYX1[J][NZ] + (Mur_Tx - 1) / (Mur_Tx + 1) * (EY[2][J][NZ] - EY[1][J][NZ]);
        EY[NX][J][1] = EYXN[J][1] + (Mur_Tx - 1) / (Mur_Tx + 1) * (EY[NX - 1][J][1] - EY[NX][J][1]);
        EY[NX][J][NZ] = EYXN[J][NZ] + (Mur_Tx - 1) / (Mur_Tx + 1) * (EY[NX - 1][J][NZ] - EY[NX][J][NZ]);
    }

    // РЕБРА, параллельные оси X
    for (I = 1; I <= NX - 1; I++)
    {
        EX[I][1][1] = EXY1[I][1] + (Mur_Ty - 1) / (Mur_Ty + 1) * (EX[I][2][1] - EX[I][1][1]);
        EX[I][1][NZ] = EXYN[I][NZ] + (Mur_Ty - 1) / (Mur_Ty + 1) * (EX[I][2][NZ] - EX[I][1][NZ]);
        EX[I][NY][1] = EXY1[I][1] + (Mur_Ty - 1) / (Mur_Ty + 1) * (EX[I][NY - 1][1] - EX[I][NY][1]);
        EX[I][NY][NZ] = EXYN[I][NZ] + (Mur_Ty - 1) / (Mur_Ty + 1) * (EX[I][NY - 1][NZ] - EX[I][NY][NZ]);
    }
}

void GetBorderValues()
{

    int NX, NY, NZ;
    std::vector<std::vector<int>> EZ(NX, std::vector<int>(NY, std::vector<int>(NZ)));
    std::vector<std::vector<int>> EY(NX, std::vector<int>(NY, std::vector<int>(NZ)));
    std::vector<std::vector<int>> EX(NX, std::vector<int>(NY, std::vector<int>(NZ)));
    std::vector<std::vector<int>> EZX1(NY, std::vector<int>(NZ));
    std::vector<std::vector<int>> EZXN(NY, std::vector<int>(NZ));
    std::vector<std::vector<int>> EYX1(NY, std::vector<int>(NZ));
    std::vector<std::vector<int>> EYXN(NY, std::vector<int>(NZ));
    std::vector<std::vector<int>> EZY1(NX, std::vector<int>(NZ));
    std::vector<std::vector<int>> EZYN(NX, std::vector<int>(NZ));
    std::vector<std::vector<int>> EXY1(NX, std::vector<int>(NZ));
    std::vector<std::vector<int>> EXYN(NX, std::vector<int>(NZ));
    std::vector<std::vector<int>> EYZ1(NX, std::vector<int>(NY - 1));
    std::vector<std::vector<int>> EYZN(NX, std::vector<int>(NY - 1));
    std::vector<std::vector<int>> EXZ1(NX - 1, std::vector<int>(NY));
    std::vector<std::vector<int>> EXZN(NX - 1, std::vector<int>(NY)) int I, J, K;

    // слева и справа
    for (J = 1; J <= NY; J++)
    {
        for (K = 1; K <= NZ - 1; K++)
        {
            EZX1[J][K] = EZ[2][J][K];
            EZXN[J][K] = EZ[NX - 1][J][K];
        }
    }
    for (J = 1; J <= NY - 1; J++)
    {
        for (K = 1; K <= NZ; K++)
        {
            EYX1[J][K] = EY[2][J][K];
            EYXN[J][K] = EY[NX - 1][J][K];
        }
    }

    // верх и низ
    for (I = 1; I <= NX; I++)
    {
        for (K = 1; K <= NZ - 1; K++)
        {
            EZY1[I][K] = EZ[I][2][K];
            EZYN[I][K] = EZ[I][NY - 1][K];
        }
    }
    for (I = 1; I <= NX - 1; I++)
    {
        for (K = 1; K <= NZ; K++)
        {
            EXY1[I][K] = EX[I][2][K];
            EXYN[I][K] = EX[I][NY - 1][K];
        }
    }

    // ближняя и дальняя
    for (I = 1; I <= NX; I++)
    {
        for (J = 1; J <= NY - 1; J++)
        {
            EYZ1[I][J] = EY[I][J][2];
            EYZN[I][J] = EY[I][J][NZ - 1];
        }
    }
    for (I = 1; I <= NX - 1; I++)
    {
        for (J = 1; J <= NY; J++)
        {
            EXZ1[I][J] = EX[I][J][2];
            EXZN[I][J] = EX[I][J][NZ - 1];
        }
    }
}
