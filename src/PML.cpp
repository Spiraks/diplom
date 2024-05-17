#include <vector>
#include <cmath>
#include "FDTD.h"
#ifdef PML
float FDTD::GetSigma(float D, float dStep, float pmlG, float pmlR,int N) {
    return -(eps0 * c) / (2 * dStep) * log(pmlG) / (pow(pmlG, N) - 1) * log(pmlR) * pow(pmlG, D / dStep);
}

void FDTD::Pml(){
    PML_E();
    PML_H();
}

void FDTD::InitPML()
{
    std::cout << "INIT PML\n";
    for (size_t K = 0; K < plmLayerN; ++K) {
        float sigma = this->GetSigma(K*dx,dx,2.15,0.001,K);
        pmlSigmaStarE[K] = 1;
        float sigma_s =(sigma *mu0)/eps0;
        std::cout << "sigma " <<sigma<< " | sigma_s " <<sigma_s<<" | "<<K<<"\n";
        pmlSigmaStarH[K] = 1;
        pmlExpSigmaStarE[K] = 1;
        pmlExpSigmaStarH[K] = 1;
    }
}

#define setPML(x, y, z, x1, y1, z1, mass) ({ \
    mass[x][y][z].Hzy = H._Z[x1][y1][z1]/2;\
    mass[x][y][z].Hxz = H._X[x1][y1][z1]/2;\
    mass[x][y][z].Hyz = H._Y[x1][y1][z1]/2;\
    mass[x][y][z].Ezy = E._Z[x1][y1][z1]/2;\
    mass[x][y][z].Exz = E._X[x1][y1][z1]/2;\
    mass[x][y][z].Eyz = E._Y[x1][y1][z1]/2;\
    mass[x][y][z].Hzx = mass[x][y][z].Hzy;\
    mass[x][y][z].Hxy = mass[x][y][z].Hxz;\
    mass[x][y][z].Hyx = mass[x][y][z].Hyz;\
    mass[x][y][z].Ezx = mass[x][y][z].Ezy;\
    mass[x][y][z].Exy = mass[x][y][z].Exz;\
    mass[x][y][z].Eyx = mass[x][y][z].Eyz; })

#define getPML(x, y, z, x1, y1, z1, mass) ({ \
    H._Z[x1][y1][z1] = mass[x][y][z].Hzx + mass[x][y][z].Hzy;\
    H._X[x1][y1][z1] = mass[x][y][z].Hxy + mass[x][y][z].Hxz;\
    H._Y[x1][y1][z1] = mass[x][y][z].Hyx + mass[x][y][z].Hyz;\
    E._Z[x1][y1][z1] = mass[x][y][z].Ezx + mass[x][y][z].Ezy;\
    E._X[x1][y1][z1] = mass[x][y][z].Exy + mass[x][y][z].Exz;\
    E._Y[x1][y1][z1] = mass[x][y][z].Eyx + mass[x][y][z].Eyz; })

void FDTD::SetBorderValuesPML(){

    size_t s = 1;
    size_t end = 1;
    size_t x = len_x - end;
    size_t y = len_y - end;
    size_t z = len_z - end;
    // ЛЕВАЯ И ПРАВАЯ ПЛОСКОСТИ
    //  Ez для двух плоскостей YxZ (левой и правой)
    // Ey для двух плоскостей YxZ [левой и правой]


    for (size_t J = s; J <= y; J++)
    {
        for (size_t K = s; K <= z; K++)
        {
            setPML(plmLayerNMinus1, J, K, s, J, K,pmlX);
            setPML(plmLayerNMinus1, J, K, x, J, K,pmlXN);
        }
    }
    // ВЕРХНЯЯ И НИЖНЯЯ ПЛОСКОСТИ
    //  Ex для двух плоскостей XxZ [верхней и нижней]
    // Ez для двух плоскостей XxZ [верхней и нижней]


    for (size_t K = s; K <= z; K++)
    {
        for (size_t I = s; I <= x; I++)
        {
            setPML(I, plmLayerNMinus1, K, I, s, K,pmlY);
            setPML(I, plmLayerNMinus1, K, I, y, K,pmlYN);
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
            setPML(I, J, plmLayerNMinus1, I, J, s,pmlZ);
            setPML(I, J, plmLayerNMinus1, I, J, x,pmlZN);
        }
    }
    

    // РЕБРА, параллельные оси Z
    // std::cout << "7 РЕБРА, параллельные оси Z\n";

    for (size_t K = s; K <= z; K++)
    {
        setPML(plmLayerNMinus1, plmLayerNMinus1, K, s, s, K,pmlXY);
        setPML(plmLayerNMinus1, plmLayerNMinus1, K, s, y, K,pmlXYN);
        setPML(plmLayerNMinus1, plmLayerNMinus1, K, x, s, K,pmlXNY);
        setPML(plmLayerNMinus1, plmLayerNMinus1, K, x, y, K,pmlXNYN);
    }
    
    // РЕБРА, параллельные оси Y
    // std::cout << "8 РЕБРА, параллельные оси Y\n";

    for (size_t J = s; J <= y; J++)
    {
        setPML(plmLayerNMinus1, J, plmLayerNMinus1, s, J, s,pmlZX);
        setPML(plmLayerNMinus1, J, plmLayerNMinus1, s, J, z,pmlZNX);
        setPML(plmLayerNMinus1, J, plmLayerNMinus1, x, J, s, pmlZXN);
        setPML(plmLayerNMinus1, J, plmLayerNMinus1, x, J, z, pmlZNXN);
    }


    // РЕБРА, параллельные оси X
    // std::cout << "9 РЕБРА, параллельные оси X\n";
    for (size_t I = s; I <= x; I++)
    {

        setPML(I, plmLayerNMinus1, plmLayerNMinus1, I, s, s,pmlZY);
        setPML(I, plmLayerNMinus1, plmLayerNMinus1, I, s, z,pmlZNY);
        setPML(I, plmLayerNMinus1, plmLayerNMinus1, I, y, s,pmlZYN);
        setPML(I, plmLayerNMinus1, plmLayerNMinus1, I, y, z,pmlZNYN);

    }
}

void FDTD::GetBorderValuesPML(){

    // std::cout << "Mur\n";
    size_t s = 1;
    size_t end = 1;
    size_t x = len_x - end;
    size_t y = len_y - end;
    size_t z = len_z - end;

    // ЛЕВАЯ И ПРАВАЯ ПЛОСКОСТИ
    // Ez для двух плоскостей YxZ (левой и правой)
    // Ey для двух плоскостей YxZ [левой и правой]

    for (size_t J = s; J <= y; J++)
    {
        for (size_t K = s; K <= z; K++)
        {
            getPML(plmLayerNMinus1, J, K, s, J, K,pmlX);
            getPML(plmLayerNMinus1, J, K, x, J, K,pmlXN);
        }
    }
    // ВЕРХНЯЯ И НИЖНЯЯ ПЛОСКОСТИ
    //  Ex для двух плоскостей XxZ [верхней и нижней]
    // Ez для двух плоскостей XxZ [верхней и нижней]

    for (size_t K = s; K <= z; K++)
    {
        for (size_t I = s; I <= x; I++)
        {
            getPML(I, plmLayerNMinus1, K, I, s, K,pmlY);
            getPML(I, plmLayerNMinus1, K, I, y, K,pmlYN);
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
            getPML( I, J, plmLayerNMinus1, I, J, s,pmlZ);
            getPML( I, J, plmLayerNMinus1, I, J, x,pmlZN);
        }
    }

    // РЕБРА, параллельные оси Z
    // std::cout << "7 РЕБРА, параллельные оси Z\n";

    for (size_t K = s; K <= z; K++)
    {
        getPML(plmLayerNMinus1, plmLayerNMinus1, K, s, s, K,pmlXY);
        getPML(plmLayerNMinus1, plmLayerNMinus1, K, s, y, K,pmlXYN);
        getPML(plmLayerNMinus1, plmLayerNMinus1, K, x, s, K,pmlXNY);
        getPML(plmLayerNMinus1, plmLayerNMinus1, K, x, y, K,pmlXNYN);
    }

    // РЕБРА, параллельные оси Y
    // std::cout << "8 РЕБРА, параллельные оси Y\n";

    for (size_t J = s; J <= y; J++)
    {
        getPML(plmLayerNMinus1, J, plmLayerNMinus1, s, J, s,pmlZX);
        getPML(plmLayerNMinus1, J, plmLayerNMinus1, s, J, z,pmlZNX);
        getPML(plmLayerNMinus1, J, plmLayerNMinus1, x, J, s, pmlZXN);
        getPML(plmLayerNMinus1, J, plmLayerNMinus1, x, J, z, pmlZNXN);
    }

    // РЕБРА, параллельные оси X
    // std::cout << "9 РЕБРА, параллельные оси X\n";
    for (size_t I = s; I <= x; I++)
    {
        getPML(I, plmLayerNMinus1, plmLayerNMinus1, I, s, s,pmlZY);
        getPML(I, plmLayerNMinus1, plmLayerNMinus1, I, s, z,pmlZNY);
        getPML(I, plmLayerNMinus1, plmLayerNMinus1, I, y, s,pmlZYN);
        getPML(I, plmLayerNMinus1, plmLayerNMinus1, I, y, z,pmlZNYN);
    }
}
#endif