#include <vector>
#include <cmath>
#include "FDTD.h"
#ifdef PML

#define START 1
#define END 2

float FDTD::GetSigma(float D, float dStep, float pmlG, float pmlR,int N) {
    return -(eps0 * c) / (2 * dStep) * log(pmlG) / (pow(pmlG, N) - 1) * log(pmlR) * pow(pmlG, D / dStep);
}

void FDTD::Pml(){
    PML_H();
    PML_E();
}

void FDTD::InitPML()
{
    std::cout << "INIT PML\n";
    for (size_t K = 1; K <= plmLayerN; ++K) {
        float _sigma = this->GetSigma(K*dx,dx,1.385,0.001,K);
        pmlSigmaStarE[K-1] = (1-0,5 * _sigma *dt)/ (1+0,5 *_sigma *dt);
        float _mu =(_sigma *mu0)/eps0;
        std::cout << "_sigma = " <<_sigma<< "| mu = " <<_mu<<"|"<<K<<"\n";
        pmlSigmaStarH[K-1] = (1-0,5 * _mu *dt)/ (1+0,5 *_mu *dt);
        pmlExpSigmaStarE[K-1] = dt/(1 + 0,5 * _sigma * dt)/dx;
        pmlExpSigmaStarH[K-1] = dt/(1 + 0,5 * _mu * dt)/dx;
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

    size_t s = START;
    size_t end = END;
    size_t x = len_x - end;
    size_t y = len_y - end;
    size_t z = len_z - end;

    // setPML(plmLayerN1, plmLayerN1, plmLayerN1, s, s, s,pmlXYZ);
    // setPML(plmLayerN1, plmLayerN1, plmLayerN1, x, s, s,pmlXNYZ);            
    // setPML(plmLayerN1, plmLayerN1, plmLayerN1, s, s, z,pmlXYZN);
    // setPML(plmLayerN1, plmLayerN1, plmLayerN1, s, y, s,pmlXYNZ);            
    // setPML(plmLayerN1, plmLayerN1, plmLayerN1, x, s, z,pmlXNYZN);
    // setPML(plmLayerN1, plmLayerN1, plmLayerN1, s, y, z,pmlXYNZN);            
    // setPML(plmLayerN1, plmLayerN1, plmLayerN1, x, y, s,pmlXNYNZ);
    // setPML(plmLayerN1, plmLayerN1, plmLayerN1, x, y, z,pmlXNYNZN);
    
    // ЛЕВАЯ И ПРАВАЯ ПЛОСКОСТИ
    //  Ez для двух плоскостей YxZ (левой и правой)
    // Ey для двух плоскостей YxZ [левой и правой]

    for (size_t J = s; J <= y; J++)
    {
        for (size_t K = s; K <= z; K++)
        {
            setPML(plmLayerN1, J, K, s, J, K,pmlX);
            setPML(plmLayerN1, J, K, x, J, K,pmlXN);
        }
    }
    // ВЕРХНЯЯ И НИЖНЯЯ ПЛОСКОСТИ
    //  Ex для двух плоскостей XxZ [верхней и нижней]
    // Ez для двух плоскостей XxZ [верхней и нижней]


    for (size_t K = s; K <= z; K++)
    {
        for (size_t I = s; I <= x; I++)
        {
            setPML(I, plmLayerN1, K, I, s, K,pmlY);
            setPML(I, plmLayerN1, K, I, y, K,pmlYN);
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
            setPML(I, J, plmLayerN1, I, J, s,pmlZ);
            setPML(I, J, plmLayerN1, I, J, z,pmlZN);
        }
    }
    

    // РЕБРА, параллельные оси Z
    // std::cout << "7 РЕБРА, параллельные оси Z\n";

    for (size_t K = s; K <= z; K++)
    {
        setPML(plmLayerN1, plmLayerN1, K, s, s, K,pmlXY);
        setPML(plmLayerN1, plmLayerN1, K, s, y, K,pmlXYN);
        setPML(plmLayerN1, plmLayerN1, K, x, s, K,pmlXNY);
        setPML(plmLayerN1, plmLayerN1, K, x, y, K,pmlXNYN);
    }
    
    // РЕБРА, параллельные оси Y
    // std::cout << "8 РЕБРА, параллельные оси Y\n";

    for (size_t J = s; J <= y; J++)
    {
        setPML(plmLayerN1, J, plmLayerN1, s, J, s,pmlZX);
        setPML(plmLayerN1, J, plmLayerN1, s, J, z,pmlZNX);
        setPML(plmLayerN1, J, plmLayerN1, x, J, s, pmlZXN);
        setPML(plmLayerN1, J, plmLayerN1, x, J, z, pmlZNXN);
    }


    // РЕБРА, параллельные оси X
    // std::cout << "9 РЕБРА, параллельные оси X\n";
    for (size_t I = s; I <= x; I++)
    {

        setPML(I, plmLayerN1, plmLayerN1, I, s, s,pmlZY);
        setPML(I, plmLayerN1, plmLayerN1, I, s, z,pmlZNY);
        setPML(I, plmLayerN1, plmLayerN1, I, y, s,pmlZYN);
        setPML(I, plmLayerN1, plmLayerN1, I, y, z,pmlZNYN);

    }
}

void FDTD::GetBorderValuesPML(){

    // std::cout << "Mur\n";
    int s = 0;
    int s1 = 1;
    int end = 1;
    int x = len_x - end;
    int y = len_y - end;
    int z = len_z - end;
    int x1 = x - 1;
    int y1 = y - 1;
    int z1 = z - 1;

    // getPML(plmLayerN1, plmLayerN1, plmLayerN1, s, s, s,pmlXYZ);
    // getPML(plmLayerN1, plmLayerN1, plmLayerN1, x, s, s,pmlXNYZ);            
    // getPML(plmLayerN1, plmLayerN1, plmLayerN1, s, s, z,pmlXYZN);
    // getPML(plmLayerN1, plmLayerN1, plmLayerN1, s, y, s,pmlXYNZ);            
    // getPML(plmLayerN1, plmLayerN1, plmLayerN1, x, s, z,pmlXNYZN);
    // getPML(plmLayerN1, plmLayerN1, plmLayerN1, s, y, z,pmlXYNZN);            
    // getPML(plmLayerN1, plmLayerN1, plmLayerN1, x, y, s,pmlXNYNZ);
    // getPML(plmLayerN1, plmLayerN1, plmLayerN1, x, y, z,pmlXNYNZN);
    // ЛЕВАЯ И ПРАВАЯ ПЛОСКОСТИ
    // Ez для двух плоскостей YxZ (левой и правой)
    // Ey для двух плоскостей YxZ [левой и правой]

    for (size_t J = s1; J <= y1; J++)
    {
        for (size_t K = s1; K <= z1; K++)
        {
            getPML(plmLayerN1, J, K, s, J, K,pmlX);
            getPML(plmLayerN1, J, K, x, J, K,pmlXN);
        }
    }
    // ВЕРХНЯЯ И НИЖНЯЯ ПЛОСКОСТИ
    //  Ex для двух плоскостей XxZ [верхней и нижней]
    // Ez для двух плоскостей XxZ [верхней и нижней]

    for (size_t K = s1; K <= z1; K++)
    {
        for (size_t I = s1; I <= x1; I++)
        {
            getPML(I, plmLayerN1, K, I, s, K,pmlY);
            getPML(I, plmLayerN1, K, I, y, K,pmlYN);
        }
    }

    // БЛИЖНЯЯ И ДАЛЬНЯЯ ПЛОСКОСТИ
    // Ey для двух плоскостей XxY [ближней и дальней]
    // Ex для двух плоскостей XxY [ближней и дальней]
    // std::cout << "5 Ey для двух плоскостей XxY [ближней и дальней]\n";
    for (size_t I = s1; I <= x1; I++)
    {
        for (size_t J = s1; J <= y1; J++)
        {
            getPML( I, J, plmLayerN1, I, J, s,pmlZ);
            getPML( I, J, plmLayerN1, I, J, x,pmlZN);
        }
    }

    // РЕБРА, параллельные оси Z
    // std::cout << "7 РЕБРА, параллельные оси Z\n";

    for (size_t K = s1; K <= z1; K++)
    {
        getPML(plmLayerN1, plmLayerN1, K, s, s, K,pmlXY);
        getPML(plmLayerN1, plmLayerN1, K, s, y, K,pmlXYN);
        getPML(plmLayerN1, plmLayerN1, K, x, s, K,pmlXNY);
        getPML(plmLayerN1, plmLayerN1, K, x, y, K,pmlXNYN);
    }

    // РЕБРА, параллельные оси Y
    // std::cout << "8 РЕБРА, параллельные оси Y\n";

    for (size_t J = s1; J <= y1; J++)
    {
        getPML(plmLayerN1, J, plmLayerN1, s, J, s,pmlZX);
        getPML(plmLayerN1, J, plmLayerN1, s, J, z,pmlZNX);
        getPML(plmLayerN1, J, plmLayerN1, x, J, s, pmlZXN);
        getPML(plmLayerN1, J, plmLayerN1, x, J, z, pmlZNXN);
    }

    // РЕБРА, параллельные оси X
    // std::cout << "9 РЕБРА, параллельные оси X\n";
    for (size_t I = s1; I <= x1; I++)
    {
        getPML(I, plmLayerN1, plmLayerN1, I, s, s,pmlZY);
        getPML(I, plmLayerN1, plmLayerN1, I, s, z,pmlZNY);
        getPML(I, plmLayerN1, plmLayerN1, I, y, s,pmlZYN);
        getPML(I, plmLayerN1, plmLayerN1, I, y, z,pmlZNYN);
    }
}
#endif