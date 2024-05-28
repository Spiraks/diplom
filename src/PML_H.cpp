#include "FDTD.h"
#ifdef PML
void FDTD::PML_H(){

    //подсчет граней
    PML_UpdateHXN();
    PML_UpdateHX();
    PML_UpdateHYN();
    PML_UpdateHY();
    PML_UpdateHZN();
    PML_UpdateHZ();

    // //подсчет ребер
    PML_UpdateHXNYN();
    PML_UpdateHXYN();
    PML_UpdateHXNY();
    PML_UpdateHXY();

    PML_UpdateHZNYN();
    PML_UpdateHZYN();
    PML_UpdateHZNY();
    PML_UpdateHZY();

    PML_UpdateHZNXN();
    PML_UpdateHZXN();
    PML_UpdateHZNX();
    PML_UpdateHZX();

    // подсчет углов
    PML_UpdateHXYZ();
    PML_UpdateHXNYZ();
    PML_UpdateHXYZN();
    PML_UpdateHXYNZ();

    PML_UpdateHXNYZN();
    PML_UpdateHXYNZN();
    PML_UpdateHXNYNZ();
    PML_UpdateHXNYNZN();
}

#define _HYX(mass,x, y, z, _s1, _s2) ({ \
    mass[x][y][z].Hyx =_s1 * mass[x][y][z].Hyx + _s2 * \
    (mass[x + 1][y][z].Ezx + mass[x + 1][y][z].Ezy \
    - mass[x][y][z].Ezx - mass[x][y][z].Ezy);})

#define _HXY(mass,x, y, z, _s1, _s2) ({ \
    mass[x][y][z].Hxy = _s1 * mass[x][y][z].Hxy + _s2 * \
    (mass[x][y+1][z].Ezx + mass[x][y+1][z].Ezy \
    - mass[x][y][z].Ezx - mass[x][y][z].Ezy);})

#define _HXZ(mass,x, y, z, _s1, _s2) ({ \
    mass[x][y][z].Hxz =_s1 * mass[x][y][z].Hxz + _s2 * \
    (mass[x][y][z+1].Eyz + mass[x][y][z+1].Eyx \
    - mass[x][y][z].Eyz - mass[x][y][z].Eyx);})

#define _HZX(mass,x, y, z, _s1, _s2) ({ \
    mass[x][y][z].Hzx = _s1 * mass[x][y][z].Hzx + _s2 * \
    (mass[x + 1][y][z].Eyz + mass[x + 1][y][z].Eyx \
    - mass[x][y][z].Eyz - mass[x][y][z].Eyx);})

#define _HZY(mass,x, y, z, _s1, _s2) ({ \
    mass[x][y][z].Hzy = _s1 * mass[x][y][z].Hzy + _s2 * \
    (mass[x][y+1][z].Exy + mass[x][y+1][z].Exz \
    - mass[x][y][z].Exy - mass[x][y][z].Exz);})

#define _HYZ(mass,x, y, z, _s1, _s2) ({ \
    mass[x][y][z].Hyz =_s1 * mass[x][y][z].Hyz + _s2 * \
    (mass[x][y][z+1].Exy + mass[x][y][z+1].Exz \
    - mass[x][y][z].Exy - mass[x][y][z].Exz);})

//____________________________
#define _0HYX(mass,x, y, z, _s1, _s2) ({ \
    mass[x][y][z].Hyx =_s1 * mass[x][y][z].Hyx + _s2 * \
    (0 - mass[x][y][z].Ezx - mass[x][y][z].Ezy);})

#define _0HXY(mass,x, y, z, _s1, _s2) ({ \
    mass[x][y][z].Hxy = _s1 * mass[x][y][z].Hxy + _s2 * \
    (0 - mass[x][y][z].Ezx - mass[x][y][z].Ezy);})

#define _0HXZ(mass,x, y, z, _s1, _s2) ({ \
    mass[x][y][z].Hxz =_s1 * mass[x][y][z].Hxz + _s2 * \
    (0 - mass[x][y][z].Eyz - mass[x][y][z].Eyx);})

#define _0HZX(mass,x, y, z, _s1, _s2) ({ \
    mass[x][y][z].Hzx = _s1 * mass[x][y][z].Hzx + _s2 * \
    (0 - mass[x][y][z].Eyz - mass[x][y][z].Eyx);})

#define _0HZY(mass,x, y, z, _s1, _s2) ({ \
    mass[x][y][z].Hzy = _s1 * mass[x][y][z].Hzy + _s2 * \
    (0 - mass[x][y][z].Exy - mass[x][y][z].Exz);})

#define _0HYZ(mass,x, y, z, _s1, _s2) ({ \
    mass[x][y][z].Hyz =_s1 * mass[x][y][z].Hyz + _s2 * \
    (0 - mass[x][y][z].Exy - mass[x][y][z].Exz);})


void FDTD::PML_UpdateHXYZ() {
    size_t I, J, K;

    I = plmLayerN1;
    J = plmLayerN1;
    K = plmLayerN1;
    _0HZX(pmlXYZ,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);
    _0HYX(pmlXYZ,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);
    _0HXY(pmlXYZ,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);
    _0HZY(pmlXYZ,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]); 
    _0HYZ(pmlXYZ,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);
    _0HXZ(pmlXYZ,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);


    for (I = 0; I < plmLayerN2; ++I) {
        for (J = 0; J < plmLayerN2; ++J) {
            for (K = 0; K < plmLayerN2; ++K) {
    // Hyz
    // std::cout << "Hxz\n";  
                _HYZ(pmlXYZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXZ(pmlXYZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    // Hyx
    // std::cout << "Hzx\n";  
                _HZX(pmlXYZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HYX(pmlXYZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    // Hzy
    // Hxy
                _HZY(pmlXYZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXY(pmlXYZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);

            }
        }
    }
}
void FDTD::PML_UpdateHXNYZ() {
    size_t I, J, K;

    I = plmLayerN1;
    J = plmLayerN1;
    K = plmLayerN1;
    _0HZX(pmlXNYZ,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);
    _0HYX(pmlXNYZ,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);

    _0HZY(pmlXNYZ,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);
    _0HXY(pmlXNYZ,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);

    _0HYZ(pmlXNYZ,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);
    _0HXZ(pmlXNYZ,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);


    for (I = 0; I < plmLayerN2; ++I) {
        for (J = 0; J < plmLayerN2; ++J) {
            for (K = 0; K < plmLayerN2; ++K) {
    // Hyz
    //std::cout << "Hxz\n";  
                _HYZ(pmlXNYZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXZ(pmlXNYZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    // Hyx
    //std::cout << "Hzx\n";  
                _HZX(pmlXNYZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HYX(pmlXNYZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    // Hzy
    // Hxy
                _HZY(pmlXNYZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXY(pmlXNYZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);

            }
        }
    }
}
void FDTD::PML_UpdateHXYNZ() {
    size_t I, J, K;

    I = plmLayerN1;
    J = plmLayerN1;
    K = plmLayerN1;
    _0HZX(pmlXYNZ,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);
    _0HYX(pmlXYNZ,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);

    _0HZY(pmlXYNZ,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);
    _0HXY(pmlXYNZ,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);

    _0HYZ(pmlXYNZ,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);
    _0HXZ(pmlXYNZ,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);


    for (I = 0; I < plmLayerN2; ++I) {
        for (J = 0; J < plmLayerN2; ++J) {
            for (K = 0; K < plmLayerN2; ++K) {
    // Hyz
    //std::cout << "Hxz\n";  
                _HYZ(pmlXYNZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXZ(pmlXYNZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    // Hyx
    //std::cout << "Hzx\n";  
                _HZX(pmlXYNZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HYX(pmlXYNZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    // Hzy
    // Hxy
                _HZY(pmlXYNZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXY(pmlXYNZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);

            }
        }
    }
}
void FDTD::PML_UpdateHXYZN() {
        size_t I, J, K;

    I = plmLayerN1;
    J = plmLayerN1;
    K = plmLayerN1;
    _0HZX(pmlXYZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    _0HYX(pmlXYZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);

    _0HZY(pmlXYZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    _0HXY(pmlXYZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);

    _0HYZ(pmlXYZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    _0HXZ(pmlXYZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);


    for (I = 0; I < plmLayerN2; ++I) {
        for (J = 0; J < plmLayerN2; ++J) {
            for (K = 0; K < plmLayerN2; ++K) {
    // Hyz
    //std::cout << "Hxz\n";  
                _HYZ(pmlXYZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXZ(pmlXYZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    // Hyx
    //std::cout << "Hzx\n";  
                _HZX(pmlXYZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HYX(pmlXYZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    // Hzy
    // Hxy
                _HZY(pmlXYZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXY(pmlXYZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);

            }
        }
    }
}
void FDTD::PML_UpdateHXNYNZ() {
        size_t I, J, K;

    I = plmLayerN1;
    J = plmLayerN1;
    K = plmLayerN1;
    _0HZX(pmlXNYNZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    _0HYX(pmlXNYNZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);

    _0HZY(pmlXNYNZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    _0HXY(pmlXNYNZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);

    _0HYZ(pmlXNYNZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    _0HXZ(pmlXNYNZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);


    for (I = 0; I < plmLayerN2; ++I) {
        for (J = 0; J < plmLayerN2; ++J) {
            for (K = 0; K < plmLayerN2; ++K) {
    // Hyz
    //std::cout << "Hxz\n";  
                _HYZ(pmlXNYNZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXZ(pmlXNYNZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    // Hyx
    //std::cout << "Hzx\n";  
                _HZX(pmlXNYNZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HYX(pmlXNYNZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    // Hzy
    // Hxy
                _HZY(pmlXNYNZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXY(pmlXNYNZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);

            }
        }
    }
}
void FDTD::PML_UpdateHXNYZN() {
        size_t I, J, K;

    I = plmLayerN1;
    J = plmLayerN1;
    K = plmLayerN1;
    _0HZX(pmlXNYZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    _0HYX(pmlXNYZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);

    _0HZY(pmlXNYZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    _0HXY(pmlXNYZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);

    _0HYZ(pmlXNYZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    _0HXZ(pmlXNYZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);


    for (I = 0; I < plmLayerN2; ++I) {
        for (J = 0; J < plmLayerN2; ++J) {
            for (K = 0; K < plmLayerN2; ++K) {
    // Hyz
    //std::cout << "Hxz\n";  
                _HYZ(pmlXNYZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXZ(pmlXNYZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    // Hyx
    //std::cout << "Hzx\n";  
                _HZX(pmlXNYZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HYX(pmlXNYZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    // Hzy
    // Hxy
                _HZY(pmlXNYZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXY(pmlXNYZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);

            }
        }
    }
}
void FDTD::PML_UpdateHXYNZN() {
        size_t I, J, K;

    I = plmLayerN1;
    J = plmLayerN1;
    K = plmLayerN1;
    _0HZX(pmlXYNZN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);
    _0HYX(pmlXYNZN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);

    _0HZY(pmlXYNZN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);
    _0HXY(pmlXYNZN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);

    _0HYZ(pmlXYNZN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);
    _0HXZ(pmlXYNZN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);


    for (I = 0; I < plmLayerN2; ++I) {
        for (J = 0; J < plmLayerN2; ++J) {
            for (K = 0; K < plmLayerN2; ++K) {
    // Hyz
    //std::cout << "Hxz\n";  
                _HYZ(pmlXYNZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXZ(pmlXYNZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    // Hyx
    //std::cout << "Hzx\n";  
                _HZX(pmlXYNZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HYX(pmlXYNZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    // Hzy
    // Hxy
                _HZY(pmlXYNZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXY(pmlXYNZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);

            }
        }
    }
}
void FDTD::PML_UpdateHXNYNZN() {
        size_t I, J, K;

    I = plmLayerN1;
    J = plmLayerN1;
    K = plmLayerN1;
    _0HZX(pmlXNYNZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    _0HYX(pmlXNYNZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);

    _0HZY(pmlXNYNZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    _0HXY(pmlXNYNZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);

    _0HYZ(pmlXNYNZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    _0HXZ(pmlXNYNZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);


    for (I = 0; I < plmLayerN2; ++I) {
        for (J = 0; J < plmLayerN2; ++J) {
            for (K = 0; K < plmLayerN2; ++K) {
    // Hyz
    //std::cout << "Hxz\n";  
                _HYZ(pmlXNYNZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXZ(pmlXNYNZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    // Hyx
    //std::cout << "Hzx\n";  
                _HZX(pmlXNYNZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HYX(pmlXNYNZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    // Hzy
    // Hxy
                _HZY(pmlXNYNZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXY(pmlXNYNZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);

            }
        }
    }
}
//подсчет ребер
void FDTD::PML_UpdateHXNYN() {
    size_t I, J, K;

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < plmLayerN1; ++I) {
        for (J = 0; J < plmLayerN1; ++J) {
            for (K = 0; K < len_z -1; ++K) {
                _HYZ(pmlXNYN,I,J,K,1,pmlBHZ);
                _HXZ(pmlXNYN,I,J,K,1,pmlBHZ);
            }
        }
    }

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I < plmLayerN2; ++I) {
        for (J = 0; J < plmLayerN2; ++J) {
            for (K = 0; K < len_y; ++K) {
                _HZX(pmlXNYN,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
                _HYX(pmlXNYN,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);

            }
        }
    }

    I = plmLayerN1;
    J = plmLayerN1;
    for (K = 0; K < len_z; ++K) {
        _0HZX(pmlXNYN,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
        _0HYX(pmlXNYN,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
    }

    // Hzy
    // Hxy
    for (I = 0; I < plmLayerN2; ++I) {
        for (J = 0; J < plmLayerN2; ++J) {
            for (K = 0; K < len_z; ++K) {
                _HZY(pmlXNYN,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
                _HXY(pmlXNYN,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);

            }
        }
    }

    I = plmLayerN1;
    J = plmLayerN1;
    for (K = 0; K < len_z; ++K) {
        _0HZY(pmlXNYN,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
        _0HXY(pmlXNYN,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
    }
}
void FDTD::PML_UpdateHXYN() {
       size_t I, J, K;

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < plmLayerN1; ++I) {
        for (J = 0; J < plmLayerN1; ++J) {
            for (K = 0; K < len_z -1; ++K) {
                _HYZ(pmlXYN,I,J,K,1,pmlBHZ);
                _HXZ(pmlXYN,I,J,K,1,pmlBHZ);
            }
        }
    }

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I < plmLayerN2; ++I) {
        for (J = 0; J < plmLayerN2; ++J) {
            for (K = 0; K < len_y; ++K) {
                _HZX(pmlXYN,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
                _HYX(pmlXYN,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);

            }
        }
    }

    I = plmLayerN1;
    J = plmLayerN1;
    for (K = 0; K < len_z; ++K) {
        _0HZX(pmlXYN,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
        _0HYX(pmlXYN,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
    }

    // Hzy
    // Hxy
    for (I = 0; I < plmLayerN2; ++I) {
        for (J = 0; J < plmLayerN2; ++J) {
            for (K = 0; K < len_z; ++K) {
                _HZY(pmlXYN,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
                _HXY(pmlXYN,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);

            }
        }
    }

    I = plmLayerN1;
    J = plmLayerN1;
    for (K = 0; K < len_z; ++K) {
        _0HZY(pmlXYN,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
        _0HXY(pmlXYN,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
    }
}
void FDTD::PML_UpdateHXNY() {
        size_t I, J, K;

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < plmLayerN1; ++I) {
        for (J = 0; J < plmLayerN1; ++J) {
            for (K = 0; K < len_z -1; ++K) {
                _HYZ(pmlXNY,I,J,K,1,pmlBHZ);
                _HXZ(pmlXNY,I,J,K,1,pmlBHZ);
            }
        }
    }

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I < plmLayerN2; ++I) {
        for (J = 0; J < plmLayerN2; ++J) {
            for (K = 0; K < len_y; ++K) {
                _HZX(pmlXNY,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
                _HYX(pmlXNY,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);

            }
        }
    }

    I = plmLayerN1;
    J = plmLayerN1;
    for (K = 0; K < len_z; ++K) {
        _0HZX(pmlXNY,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
        _0HYX(pmlXNY,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
    }

    // Hzy
    // Hxy
    for (I = 0; I < plmLayerN2; ++I) {
        for (J = 0; J < plmLayerN2; ++J) {
            for (K = 0; K < len_z; ++K) {
                _HZY(pmlXNY,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
                _HXY(pmlXNY,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);

            }
        }
    }

    I = plmLayerN1;
    J = plmLayerN1;
    for (K = 0; K < len_z; ++K) {
        _0HZY(pmlXNY,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
        _0HXY(pmlXNY,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
    }
}
void FDTD::PML_UpdateHXY() {
    size_t I, J, K;

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < plmLayerN1; ++I) {
        for (J = 0; J < plmLayerN1; ++J) {
            for (K = 0; K < len_z -1; ++K) {
                _HYZ(pmlXY,I,J,K,1,pmlBHZ);
                _HXZ(pmlXY,I,J,K,1,pmlBHZ);
            }
        }
    }

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I < plmLayerN2; ++I) {
        for (J = 0; J < plmLayerN2; ++J) {
            for (K = 0; K < len_y; ++K) {
                _HZX(pmlXY,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
                _HYX(pmlXY,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);

            }
        }
    }

    I = plmLayerN1;
    J = plmLayerN1;
    for (K = 0; K < len_z; ++K) {
        _0HZX(pmlXY,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
        _0HYX(pmlXY,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
    }

    // Hzy
    // Hxy
    for (I = 0; I < plmLayerN2; ++I) {
        for (J = 0; J < plmLayerN2; ++J) {
            for (K = 0; K < len_z; ++K) {
                _HZY(pmlXY,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
                _HXY(pmlXY,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);

            }
        }
    }

    I = plmLayerN1;
    J = plmLayerN1;
    for (K = 0; K < len_z; ++K) {
        _0HZY(pmlXY,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
        _0HXY(pmlXY,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
    }
}

void FDTD::PML_UpdateHZY() {
    size_t I, J, K;

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I <len_x-1; ++I) {
        for (J = 0; J < plmLayerN1; ++J) {
            for (K = 0; K < plmLayerN1; ++K) {
                _HYX(pmlZY,I,J,K,1,pmlBHX);
                _HZX(pmlZY,I,J,K,1,pmlBHX);
            }
        }
    }

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerN2; ++J) {
            for (K = 0; K < plmLayerN2; ++K) {
                _HYZ(pmlZY,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXZ(pmlZY,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
            }
        }
    }

    J = plmLayerN1;
    K = plmLayerN1;
    for (I = 0; I < len_x; ++I) {
        _0HYZ(pmlZY,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
        _0HXZ(pmlZY,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    }

    // Hzy
    // Hxy
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerN2; ++J) {
            for (K = 0; K < plmLayerN2; ++K) {
                _HZY(pmlZY,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXY(pmlZY,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
            }
        }
    }
    J = plmLayerN1;
    K = plmLayerN1;
    for (I = 0; I < len_x; ++I) {
        _0HZY(pmlZY,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
        _0HXY(pmlZY,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    }
}
void FDTD::PML_UpdateHZYN() {
    size_t I, J, K;

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I <len_x-1; ++I) {
        for (J = 0; J < plmLayerN1; ++J) {
            for (K = 0; K < plmLayerN1; ++K) {
                _HYX(pmlZYN,I,J,K,1,pmlBHX);
                _HZX(pmlZYN,I,J,K,1,pmlBHX);
            }
        }
    }

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerN2; ++J) {
            for (K = 0; K < plmLayerN2; ++K) {
                _HYZ(pmlZYN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXZ(pmlZYN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
            }
        }
    }

    J = plmLayerN1;
    K = plmLayerN1;
    for (I = 0; I < len_x; ++I) {
        _0HYZ(pmlZYN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
        _0HXZ(pmlZYN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    }

    // Hzy
    // Hxy
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerN2; ++J) {
            for (K = 0; K < plmLayerN2; ++K) {
                _HZY(pmlZYN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXY(pmlZYN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
            }
        }
    }
    J = plmLayerN1;
    K = plmLayerN1;
    for (I = 0; I < len_x; ++I) {
        _0HZY(pmlZYN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
        _0HXY(pmlZYN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    }
}
void FDTD::PML_UpdateHZNY() {
    size_t I, J, K;

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I <len_x-1; ++I) {
        for (J = 0; J < plmLayerN1; ++J) {
            for (K = 0; K < plmLayerN1; ++K) {
                _HYX(pmlZNY,I,J,K,1,pmlBHX);
                _HZX(pmlZNY,I,J,K,1,pmlBHX);
            }
        }
    }

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerN2; ++J) {
            for (K = 0; K < plmLayerN2; ++K) {
                _HYZ(pmlZNY,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXZ(pmlZNY,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
            }
        }
    }

    J = plmLayerN1;
    K = plmLayerN1;
    for (I = 0; I < len_x; ++I) {
        _0HYZ(pmlZNY,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
        _0HXZ(pmlZNY,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    }

    // Hzy
    // Hxy
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerN2; ++J) {
            for (K = 0; K < plmLayerN2; ++K) {
                _HZY(pmlZNY,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXY(pmlZNY,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
            }
        }
    }
    J = plmLayerN1;
    K = plmLayerN1;
    for (I = 0; I < len_x; ++I) {
        _0HZY(pmlZNY,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
        _0HXY(pmlZNY,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    }
}
void FDTD::PML_UpdateHZNYN() {
size_t I, J, K;

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I <len_x-1; ++I) {
        for (J = 0; J < plmLayerN1; ++J) {
            for (K = 0; K < plmLayerN1; ++K) {
                _HYX(pmlZNYN,I,J,K,1,pmlBHX);
                _HZX(pmlZNYN,I,J,K,1,pmlBHX);
            }
        }
    }

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerN2; ++J) {
            for (K = 0; K < plmLayerN2; ++K) {
                _HYZ(pmlZNYN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXZ(pmlZNYN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
            }
        }
    }

    J = plmLayerN1;
    K = plmLayerN1;
    for (I = 0; I < len_x; ++I) {
        _0HYZ(pmlZNYN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
        _0HXZ(pmlZNYN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    }

    // Hzy
    // Hxy
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerN2; ++J) {
            for (K = 0; K < plmLayerN2; ++K) {
                _HZY(pmlZNYN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXY(pmlZNYN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
            }
        }
    }
    J = plmLayerN1;
    K = plmLayerN1;
    for (I = 0; I < len_x; ++I) {
        _0HZY(pmlZNYN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
        _0HXY(pmlZNYN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    }
}

void FDTD::PML_UpdateHZNXN() {
    size_t I, J, K;

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I <plmLayerN2; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerN2; ++K) {
                _HYX(pmlZNXN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HZX(pmlZNXN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
            }
        }
    }

    I = plmLayerN1;
    K = plmLayerN1;
    for (J = 0; J < len_y; ++J) {
        _0HYX(pmlZNXN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
        _0HZX(pmlZNXN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    }

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < plmLayerN2; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerN2; ++K) {
                _HYZ(pmlZNXN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXZ(pmlZNXN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
            }
        }
    }

    I = plmLayerN1;
    K = plmLayerN1;
    for (J = 0; J < len_y; ++J) {
        _0HYZ(pmlZNXN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
        _0HXZ(pmlZNXN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    }
    

    // Hzy
    // Hxy
    for (I = 0; I < plmLayerN1; ++I) {
        for (J = 0; J < len_y-1; ++J) {
            for (K = 0; K < plmLayerN1; ++K) {
                _HZY(pmlZNXN,I,J,K,1,pmlBHY);
                _HXY(pmlZNXN,I,J,K,1,pmlBHY);
            }
        }
    }
}
void FDTD::PML_UpdateHZNX() {
    size_t I, J, K;

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I <plmLayerN2; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerN2; ++K) {
                _HYX(pmlZNX,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HZX(pmlZNX,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
            }
        }
    }

    I = plmLayerN1;
    K = plmLayerN1;
    for (J = 0; J < len_y; ++J) {
        _0HYX(pmlZNX,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
        _0HZX(pmlZNX,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    }

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < plmLayerN2; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerN2; ++K) {
                _HYZ(pmlZNX,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXZ(pmlZNX,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
            }
        }
    }

    I = plmLayerN1;
    K = plmLayerN1;
    for (J = 0; J < len_y; ++J) {
        _0HYZ(pmlZNX,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
        _0HXZ(pmlZNX,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    }
    

    // Hzy
    // Hxy
    for (I = 0; I < plmLayerN1; ++I) {
        for (J = 0; J < len_y-1; ++J) {
            for (K = 0; K < plmLayerN1; ++K) {
                _HZY(pmlZNX,I,J,K,1,pmlBHY);
                _HXY(pmlZNX,I,J,K,1,pmlBHY);
            }
        }
    }
}
void FDTD::PML_UpdateHZXN() {
    size_t I, J, K;

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I <plmLayerN2; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerN2; ++K) {
                _HYX(pmlZXN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HZX(pmlZXN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
            }
        }
    }

    I = plmLayerN1;
    K = plmLayerN1;
    for (J = 0; J < len_y; ++J) {
        _0HYX(pmlZXN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
        _0HZX(pmlZXN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    }

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < plmLayerN2; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerN2; ++K) {
                _HYZ(pmlZXN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXZ(pmlZXN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
            }
        }
    }

    I = plmLayerN1;
    K = plmLayerN1;
    for (J = 0; J < len_y; ++J) {
        _0HYZ(pmlZXN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
        _0HXZ(pmlZXN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    }
    

    // Hzy
    // Hxy
    for (I = 0; I < plmLayerN1; ++I) {
        for (J = 0; J < len_y-1; ++J) {
            for (K = 0; K < plmLayerN1; ++K) {
                _HZY(pmlZXN,I,J,K,1,pmlBHY);
                _HXY(pmlZXN,I,J,K,1,pmlBHY);
            }
        }
    }
}
void FDTD::PML_UpdateHZX() {
    size_t I, J, K;

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I <plmLayerN2; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerN2; ++K) {
                _HYX(pmlZX,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HZX(pmlZX,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
            }
        }
    }

    I = plmLayerN1;
    K = plmLayerN1;
    for (J = 0; J < len_y; ++J) {
        _0HYX(pmlZX,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
        _0HZX(pmlZX,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    }

    // Hyz
    //std::cout << "Hxz\n";
    for (I = 0; I < plmLayerN2; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerN2; ++K) {
                _HYZ(pmlZX,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXZ(pmlZX,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
            }
        }
    }

    I = plmLayerN1;
    K = plmLayerN1;
    for (J = 0; J < len_y; ++J) {
        _0HYZ(pmlZX,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
        _0HXZ(pmlZX,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
    }
    
    // Hzy
    // Hxy
    for (I = 0; I < plmLayerN1; ++I) {
        for (J = 0; J < len_y-1; ++J) {
            for (K = 0; K < plmLayerN1; ++K) {
                _HZY(pmlZX,I,J,K,1,pmlBHY);
                _HXY(pmlZX,I,J,K,1,pmlBHY);
            }
        }
    }
}

//подсчет граней
void FDTD::PML_UpdateHZ() {
    size_t I, J, K;

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I <len_x-1; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerN1; ++K) {
                _HYX(pmlZ,I,J,K,1,pmlBHX);
                _HZX(pmlZ,I,J,K,1,pmlBHX);
            }
        }
    }

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerN2; ++K) {
                _HYZ(pmlZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXZ(pmlZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
            }
        }
    }

    K = plmLayerN1;
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < len_y; ++J) {
            _0HYZ(pmlZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
            _0HXZ(pmlZ,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
        }
    }

    // Hzy
    // Hxy
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < len_y-1; ++J) {
            for (K = 0; K < plmLayerN1; ++K) {
                _HZY(pmlZ,I,J,K,1,pmlBHY);
                _HXY(pmlZ,I,J,K,1,pmlBHY);
            }
        }
    }
}
void FDTD::PML_UpdateHZN() {
size_t I, J, K;

    // Hyx
    //std::cout << "Hzx\n";  
    
    for (I = 0; I < len_x-1; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerN1; ++K) {
                _HYX(pmlZN,I,J,K,1,pmlBHX);
                _HZX(pmlZN,I,J,K,1,pmlBHX);
            }
        }
    }

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerN2; ++K) {
                _HYZ(pmlZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
                _HXZ(pmlZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
            }
        }
    }

    K = plmLayerN1;
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < len_y; ++J) {
            _0HYZ(pmlZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
            _0HXZ(pmlZN,I,J,K,pmlSigmaStarH[K],pmlExpSigmaStarH[K]);
        }
    }

    // Hzy
    // Hxy
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < len_y-1; ++J) {
            for (K = 0; K < plmLayerN1; ++K) {
                _HZY(pmlZN,I,J,K,1,pmlBHY);
                _HXY(pmlZN,I,J,K,1,pmlBHY);
            }
        }
    }
}
void FDTD::PML_UpdateHY() {
    size_t I, J, K;
    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I <len_x -1; ++I) {
        for (J = 0; J < plmLayerN1; ++J) {
            for (K = 0; K < len_z; ++K) {
                _HYX(pmlY,I,J,K,1,pmlBHX);
                _HZX(pmlY,I,J,K,1,pmlBHX);
            }
        }
    }

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerN1; ++J) {
            for (K = 0; K < len_z -1; ++K) {
                _HYZ(pmlY,I,J,K,1,pmlBHZ);
                _HXZ(pmlY,I,J,K,1,pmlBHZ);

            }
        }
    }

    // Hzy
    // Hxy
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerN2; ++J) {
            for (K = 0; K < len_z; ++K) {
                _HZY(pmlY,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
                _HXY(pmlY,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
            }
        }
    }

    J = plmLayerN1;
    for (I = 0; I < len_x; ++I) {
        for (K = 0; K < len_z; ++K) {
            _0HZY(pmlY,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
            _0HXY(pmlY,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
        }
    }
}
void FDTD::PML_UpdateHYN() {
    size_t I, J, K;

    // Hyx
    //std::cout << "Hzx\n";  

    for (I = 0; I <len_x - 1 ;++I) {
        for (J = 0; J < plmLayerN1; ++J) {
            for (K = 0; K < len_y; ++K) {
                _HYX(pmlYN,I,J,K,1,pmlBHX);
                _HZX(pmlYN,I,J,K,1,pmlBHX);
            }
        }
    }

    // Hyz
    //std::cout << "Hxz\n";  

    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerN1; ++J) {
            for (K = 0; K < len_z -1; ++K) {
                _HYZ(pmlYN,I,J,K,1,pmlBHZ);
                _HXZ(pmlYN,I,J,K,1,pmlBHZ);

            }
        }
    }

    // Hzy
    // Hxy

    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerN2; ++J) {
            for (K = 0; K < len_y; ++K) {
                _HZY(pmlYN,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
                _HXY(pmlYN,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
            }
        }
    }

    J = plmLayerN1;
    for (I = 0; I < len_x; ++I) {
        for (K = 0; K < len_z; ++K) {
            _0HZY(pmlYN,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
            _0HXY(pmlYN,I,J,K,pmlSigmaStarH[J],pmlExpSigmaStarH[J]);
        }
    }
}
void FDTD::PML_UpdateHX() {
    size_t I, J, K;

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < plmLayerN1; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < len_z -1; ++K) {
                _HYZ(pmlX,I,J,K,1,pmlBHZ);
                _HXZ(pmlX,I,J,K,1,pmlBHZ);

            }
        }
    }

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I < plmLayerN2; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < len_y; ++K) {
                _HZX(pmlX,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);
                _HYX(pmlX,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);

            }
        }
    }

    I = plmLayerN1;
    for (J = 0; J < len_y; ++J) {
        for (K = 0; K < len_y; ++K) {
            _0HZX(pmlX,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);
            _0HYX(pmlX,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);
        }
    }

    // Hzy
    // Hxy
    for (I = 0; I < plmLayerN1; ++I) {
        for (J = 0; J < len_y -1; ++J) {
            for (K = 0; K < len_y; ++K) {
                _HZY(pmlX,I,J,K,1,pmlBHY);
                _HXY(pmlX,I,J,K,1,pmlBHY);

            }
        }
    }
}
void FDTD::PML_UpdateHXN() {
    size_t I, J, K;

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < plmLayerN1; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < len_z -1; ++K) {
                _HYZ(pmlXN,I,J,K,1,pmlBHZ);
                _HXZ(pmlXN,I,J,K,1,pmlBHZ);
            }
        }
    }

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I < plmLayerN2; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < len_y; ++K) {
                _HZX(pmlXN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);
                _HYX(pmlXN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);

            }
        }
    }

    I = plmLayerN1;
    for (J = 0; J < len_y; ++J) {
        for (K = 0; K < len_y; ++K) {
            _0HZX(pmlXN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);
            _0HYX(pmlXN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarH[I]);
        }
    }

    // Hzy
    // Hxy
    for (I = 0; I < plmLayerN1; ++I) {
        for (J = 0; J < len_y -1; ++J) {
            for (K = 0; K < len_y; ++K) {
                _HZY(pmlXN,I,J,K,1,pmlBHY);
                _HXY(pmlXN,I,J,K,1,pmlBHY);

            }
        }
    }
}
#endif