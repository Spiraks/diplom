#include "FDTD.h"
#ifdef PML
#define _EYX(mass,x, y, z, _s1, _s2) ({ \
    mass[x][y][z].Eyx =_s1 * mass[x][y][z].Eyx + _s2 * \
    (mass[x + 1][y][z].Hzx + mass[x + 1][y][z].Hzy \
    - mass[x][y][z].Hzx - mass[x][y][z].Hzy);})

#define _EXY(mass,x, y, z, _s1, _s2) ({ \
    mass[x][y][z].Exy = _s1 * mass[x][y][z].Exy + _s2 * \
    (mass[x][y+1][z].Hzx + mass[x][y+1][z].Hzy \
    - mass[x][y][z].Hzx - mass[x][y][z].Hzy);})

#define _EXZ(mass,x, y, z, _s1, _s2) ({ \
    mass[x][y][z].Exz =_s1 * mass[x][y][z].Exz + _s2 * \
    (mass[x][y][z+1].Hyz + mass[x][y][z+1].Hyx \
    - mass[x][y][z].Hyz - mass[x][y][z].Hyx);})

#define _EZX(mass,x, y, z, _s1, _s2) ({ \
    mass[x][y][z].Ezx = _s1 * mass[x][y][z].Ezx + _s2 * \
    (mass[x + 1][y][z].Hyz + mass[x + 1][y][z].Hyx \
    - mass[x][y][z].Hyz - mass[x][y][z].Hyx);})

#define _EZY(mass,x, y, z, _s1, _s2) ({ \
    mass[x][y][z].Ezy = _s1 * mass[x][y][z].Ezy + _s2 * \
    (mass[x][y+1][z].Hxy + mass[x][y+1][z].Hxz \
    - mass[x][y][z].Hxy - mass[x][y][z].Hxz);})

#define _EYZ(mass,x, y, z, _s1, _s2) ({ \
    mass[x][y][z].Eyz =_s1 * mass[x][y][z].Eyz + _s2 * \
    (mass[x][y][z+1].Hxy + mass[x][y][z+1].Hxz \
    - mass[x][y][z].Hxy - mass[x][y][z].Hxz);})

//____________________________
#define _0EYX(mass,x, y, z, _s1, _s2) ({ \
    mass[x][y][z].Eyx =_s1 * mass[x][y][z].Eyx + _s2 * \
    (0 - mass[x][y][z].Hzx - mass[x][y][z].Hzy);})

#define _0EXY(mass,x, y, z, _s1, _s2) ({ \
    mass[x][y][z].Exy = _s1 * mass[x][y][z].Exy + _s2 * \
    (0 - mass[x][y][z].Hzx - mass[x][y][z].Hzy);})

#define _0EXZ(mass,x, y, z, _s1, _s2) ({ \
    mass[x][y][z].Exz =_s1 * mass[x][y][z].Exz + _s2 * \
    (0 - mass[x][y][z].Hyz - mass[x][y][z].Hyx);})

#define _0EZX(mass,x, y, z, _s1, _s2) ({ \
    mass[x][y][z].Ezx = _s1 * mass[x][y][z].Ezx + _s2 * \
    (0 - mass[x][y][z].Hyz - mass[x][y][z].Hyx);})

#define _0EZY(mass,x, y, z, _s1, _s2) ({ \
    mass[x][y][z].Ezy = _s1 * mass[x][y][z].Ezy + _s2 * \
    (0 - mass[x][y][z].Hxy - mass[x][y][z].Hxz);})

#define _0EYZ(mass,x, y, z, _s1, _s2) ({ \
    mass[x][y][z].Eyz =_s1 * mass[x][y][z].Eyz + _s2 * \
    (0 - mass[x][y][z].Hxy - mass[x][y][z].Hxz);})

void FDTD::PML_E(){
    //подсчет граней
        // std::cout << "подсчет граней\n";
    // std::cout << "EXN\n";
    PML_UpdateEXN();

    // std::cout << "EX\n";
    PML_UpdateEX();

    // std::cout << "EYN\n";
    PML_UpdateEYN();

    // std::cout << "EY\n";
    PML_UpdateEY();

    // std::cout << "EZN\n";
    PML_UpdateEZN();
    // std::cout << "EZ\n";

    PML_UpdateEZ();
    //подсчет ребер
        // std::cout << "подсчет ребер\n";
    // std::cout << "EXNYN\n";
    PML_UpdateEXNYN();
    // std::cout << "EXYN\n";
    PML_UpdateEXYN();
    // std::cout << "EXNY\n";
    PML_UpdateEXNY();
    // std::cout << "EXY\n";
    PML_UpdateEXY();

    // std::cout << "PML_UpdateEZY\n";
    PML_UpdateEZY();
    // std::cout << "PML_UpdateEZYN\n";
    PML_UpdateEZYN();
    // std::cout << "PML_UpdateEZNY\n";
    PML_UpdateEZNY();
    // std::cout << "PML_UpdateEZNYN\n";
    PML_UpdateEZNYN();

    // std::cout << "PML_UpdateEZNXN\n";
    PML_UpdateEZNXN();
    // std::cout << "PML_UpdateEZNX\n";
    PML_UpdateEZNX();
    // std::cout << "PML_UpdateEZXN\n";
    PML_UpdateEZXN();
    // std::cout << "PML_UpdateEZX\n";
    PML_UpdateEZX();


    // подсчет углов
    // std::cout << "подсчет углов\n";

    // std::cout << "PML_UpdateEXYZ\n";
    PML_UpdateEXYZ();
        // std::cout << "PML_UpdateEXNYZ\n";
    PML_UpdateEXNYZ();

    // std::cout << "PML_UpdateEXYNZ\n";
    PML_UpdateEXYNZ();

    // std::cout << "PML_UpdateEXYZN\n";
    PML_UpdateEXYZN();

    // std::cout << "PML_UpdateEXNYNZ\n";
    PML_UpdateEXNYNZ();

    // std::cout << "PML_UpdateEXNYZN\n";
    PML_UpdateEXNYZN();

    // std::cout << "PML_UpdateEXYNZN\n";
    PML_UpdateEXYNZN();

    // std::cout << "PML_UpdateEXNYNZN\n";
    PML_UpdateEXNYNZN();
}
void FDTD::PML_UpdateEXYZ() {
    size_t I, J, K;

    I = plmLayerNMinus1;
    J = plmLayerNMinus1;
    K = plmLayerNMinus1;
    _0EZX(pmlXYZ,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);
    _0EYX(pmlXYZ,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);
    _0EXY(pmlXYZ,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);
    _0EZY(pmlXYZ,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]); 
    _0EYZ(pmlXYZ,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);
    _0EXZ(pmlXYZ,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);


    for (I = 0; I < plmLayerNMinus2; ++I) {
        for (J = 0; J < plmLayerNMinus2; ++J) {
            for (K = 0; K < plmLayerNMinus2; ++K) {
    // Hyz
    // std::cout << "Hxz\n";  
                _EYZ(pmlXYZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXZ(pmlXYZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    // Hyx
    // std::cout << "Hzx\n";  
                _EZX(pmlXYZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EYX(pmlXYZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    // Hzy
    // Hxy
                _EZY(pmlXYZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXY(pmlXYZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);

            }
        }
    }
}
void FDTD::PML_UpdateEXNYZ() {
    size_t I, J, K;

    I = plmLayerNMinus1;
    J = plmLayerNMinus1;
    K = plmLayerNMinus1;
    _0EZX(pmlXNYZ,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);
    _0EYX(pmlXNYZ,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);

    _0EZY(pmlXNYZ,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);
    _0EXY(pmlXNYZ,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);

    _0EYZ(pmlXNYZ,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);
    _0EXZ(pmlXNYZ,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);


    for (I = 0; I < plmLayerNMinus2; ++I) {
        for (J = 0; J < plmLayerNMinus2; ++J) {
            for (K = 0; K < plmLayerNMinus2; ++K) {
    // Hyz
    //std::cout << "Hxz\n";  
                _EYZ(pmlXNYZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXZ(pmlXNYZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    // Hyx
    //std::cout << "Hzx\n";  
                _EZX(pmlXNYZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EYX(pmlXNYZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    // Hzy
    // Hxy
                _EZY(pmlXNYZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXY(pmlXNYZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);

            }
        }
    }
}
void FDTD::PML_UpdateEXYNZ() {
    size_t I, J, K;

    I = plmLayerNMinus1;
    J = plmLayerNMinus1;
    K = plmLayerNMinus1;
    _0EZX(pmlXYNZ,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);
    _0EYX(pmlXYNZ,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);

    _0EZY(pmlXYNZ,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);
    _0EXY(pmlXYNZ,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);

    _0EYZ(pmlXYNZ,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);
    _0EXZ(pmlXYNZ,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);


    for (I = 0; I < plmLayerNMinus2; ++I) {
        for (J = 0; J < plmLayerNMinus2; ++J) {
            for (K = 0; K < plmLayerNMinus2; ++K) {
    // Hyz
    //std::cout << "Hxz\n";  
                _EYZ(pmlXYNZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXZ(pmlXYNZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    // Hyx
    //std::cout << "Hzx\n";  
                _EZX(pmlXYNZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EYX(pmlXYNZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    // Hzy
    // Hxy
                _EZY(pmlXYNZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXY(pmlXYNZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);

            }
        }
    }
}
void FDTD::PML_UpdateEXYZN() {
        size_t I, J, K;

    I = plmLayerNMinus1;
    J = plmLayerNMinus1;
    K = plmLayerNMinus1;
    _0EZX(pmlXYZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    _0EYX(pmlXYZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);

    _0EZY(pmlXYZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    _0EXY(pmlXYZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);

    _0EYZ(pmlXYZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    _0EXZ(pmlXYZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);


    for (I = 0; I < plmLayerNMinus2; ++I) {
        for (J = 0; J < plmLayerNMinus2; ++J) {
            for (K = 0; K < plmLayerNMinus2; ++K) {
    // Hyz
    //std::cout << "Hxz\n";  
                _EYZ(pmlXYZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXZ(pmlXYZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    // Hyx
    //std::cout << "Hzx\n";  
                _EZX(pmlXYZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EYX(pmlXYZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    // Hzy
    // Hxy
                _EZY(pmlXYZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXY(pmlXYZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);

            }
        }
    }
}
void FDTD::PML_UpdateEXNYNZ() {
        size_t I, J, K;

    I = plmLayerNMinus1;
    J = plmLayerNMinus1;
    K = plmLayerNMinus1;
    _0EZX(pmlXNYNZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    _0EYX(pmlXNYNZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);

    _0EZY(pmlXNYNZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    _0EXY(pmlXNYNZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);

    _0EYZ(pmlXNYNZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    _0EXZ(pmlXNYNZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);


    for (I = 0; I < plmLayerNMinus2; ++I) {
        for (J = 0; J < plmLayerNMinus2; ++J) {
            for (K = 0; K < plmLayerNMinus2; ++K) {
    // Hyz
    //std::cout << "Hxz\n";  
                _EYZ(pmlXNYNZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXZ(pmlXNYNZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    // Hyx
    //std::cout << "Hzx\n";  
                _EZX(pmlXNYNZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EYX(pmlXNYNZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    // Hzy
    // Hxy
                _EZY(pmlXNYNZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXY(pmlXNYNZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);

            }
        }
    }
}
void FDTD::PML_UpdateEXNYZN() {
        size_t I, J, K;

    I = plmLayerNMinus1;
    J = plmLayerNMinus1;
    K = plmLayerNMinus1;
    _0EZX(pmlXNYZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    _0EYX(pmlXNYZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);

    _0EZY(pmlXNYZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    _0EXY(pmlXNYZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);

    _0EYZ(pmlXNYZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    _0EXZ(pmlXNYZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);


    for (I = 0; I < plmLayerNMinus2; ++I) {
        for (J = 0; J < plmLayerNMinus2; ++J) {
            for (K = 0; K < plmLayerNMinus2; ++K) {
    // Hyz
    //std::cout << "Hxz\n";  
                _EYZ(pmlXNYZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXZ(pmlXNYZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    // Hyx
    //std::cout << "Hzx\n";  
                _EZX(pmlXNYZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EYX(pmlXNYZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    // Hzy
    // Hxy
                _EZY(pmlXNYZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXY(pmlXNYZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);

            }
        }
    }
}
void FDTD::PML_UpdateEXYNZN() {
        size_t I, J, K;

    I = plmLayerNMinus1;
    J = plmLayerNMinus1;
    K = plmLayerNMinus1;
    _0EZX(pmlXYNZN,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);
    _0EYX(pmlXYNZN,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);

    _0EZY(pmlXYNZN,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);
    _0EXY(pmlXYNZN,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);

    _0EYZ(pmlXYNZN,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);
    _0EXZ(pmlXYNZN,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);


    for (I = 0; I < plmLayerNMinus2; ++I) {
        for (J = 0; J < plmLayerNMinus2; ++J) {
            for (K = 0; K < plmLayerNMinus2; ++K) {
    // Hyz
    //std::cout << "Hxz\n";  
                _EYZ(pmlXYNZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXZ(pmlXYNZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    // Hyx
    //std::cout << "Hzx\n";  
                _EZX(pmlXYNZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EYX(pmlXYNZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    // Hzy
    // Hxy
                _EZY(pmlXYNZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXY(pmlXYNZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);

            }
        }
    }
}
void FDTD::PML_UpdateEXNYNZN() {
        size_t I, J, K;

    I = plmLayerNMinus1;
    J = plmLayerNMinus1;
    K = plmLayerNMinus1;
    _0EZX(pmlXNYNZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    _0EYX(pmlXNYNZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);

    _0EZY(pmlXNYNZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    _0EXY(pmlXNYNZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);

    _0EYZ(pmlXNYNZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    _0EXZ(pmlXNYNZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);


    for (I = 0; I < plmLayerNMinus2; ++I) {
        for (J = 0; J < plmLayerNMinus2; ++J) {
            for (K = 0; K < plmLayerNMinus2; ++K) {
    // Hyz
    //std::cout << "Hxz\n";  
                _EYZ(pmlXNYNZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXZ(pmlXNYNZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    // Hyx
    //std::cout << "Hzx\n";  
                _EZX(pmlXNYNZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EYX(pmlXNYNZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    // Hzy
    // Hxy
                _EZY(pmlXNYNZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXY(pmlXNYNZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);

            }
        }
    }
}
//подсчет ребер
void FDTD::PML_UpdateEXNYN() {
    size_t I, J, K;

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < plmLayerNMinus1; ++I) {
        for (J = 0; J < plmLayerNMinus1; ++J) {
            for (K = 0; K < len_z -1; ++K) {
                _EYZ(pmlXNYN,I,J,K,1,pmlBHZ);
                _EXZ(pmlXNYN,I,J,K,1,pmlBHZ);
            }
        }
    }

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I < plmLayerNMinus2; ++I) {
        for (J = 0; J < plmLayerNMinus2; ++J) {
            for (K = 0; K < len_y; ++K) {
                _EZX(pmlXNYN,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
                _EYX(pmlXNYN,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);

            }
        }
    }

    I = plmLayerNMinus1;
    J = plmLayerNMinus1;
    for (K = 0; K < len_z; ++K) {
        _0EZX(pmlXNYN,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
        _0EYX(pmlXNYN,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
    }

    // Hzy
    // Hxy
    for (I = 0; I < plmLayerNMinus2; ++I) {
        for (J = 0; J < plmLayerNMinus2; ++J) {
            for (K = 0; K < len_z; ++K) {
                _EZY(pmlXNYN,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
                _EXY(pmlXNYN,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);

            }
        }
    }

    I = plmLayerNMinus1;
    J = plmLayerNMinus1;
    for (K = 0; K < len_z; ++K) {
        _0EZY(pmlXNYN,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
        _0EXY(pmlXNYN,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
    }
}
void FDTD::PML_UpdateEXYN() {
       size_t I, J, K;

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < plmLayerNMinus1; ++I) {
        for (J = 0; J < plmLayerNMinus1; ++J) {
            for (K = 0; K < len_z -1; ++K) {
                _EYZ(pmlXYN,I,J,K,1,pmlBHZ);
                _EXZ(pmlXYN,I,J,K,1,pmlBHZ);
            }
        }
    }

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I < plmLayerNMinus2; ++I) {
        for (J = 0; J < plmLayerNMinus2; ++J) {
            for (K = 0; K < len_y; ++K) {
                _EZX(pmlXYN,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
                _EYX(pmlXYN,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);

            }
        }
    }

    I = plmLayerNMinus1;
    J = plmLayerNMinus1;
    for (K = 0; K < len_z; ++K) {
        _0EZX(pmlXYN,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
        _0EYX(pmlXYN,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
    }

    // Hzy
    // Hxy
    for (I = 0; I < plmLayerNMinus2; ++I) {
        for (J = 0; J < plmLayerNMinus2; ++J) {
            for (K = 0; K < len_z; ++K) {
                _EZY(pmlXYN,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
                _EXY(pmlXYN,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);

            }
        }
    }

    I = plmLayerNMinus1;
    J = plmLayerNMinus1;
    for (K = 0; K < len_z; ++K) {
        _0EZY(pmlXYN,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
        _0EXY(pmlXYN,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
    }
}
void FDTD::PML_UpdateEXNY() {
        size_t I, J, K;

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < plmLayerNMinus1; ++I) {
        for (J = 0; J < plmLayerNMinus1; ++J) {
            for (K = 0; K < len_z -1; ++K) {
                _EYZ(pmlXNY,I,J,K,1,pmlBHZ);
                _EXZ(pmlXNY,I,J,K,1,pmlBHZ);
            }
        }
    }

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I < plmLayerNMinus2; ++I) {
        for (J = 0; J < plmLayerNMinus2; ++J) {
            for (K = 0; K < len_y; ++K) {
                _EZX(pmlXNY,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
                _EYX(pmlXNY,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);

            }
        }
    }

    I = plmLayerNMinus1;
    J = plmLayerNMinus1;
    for (K = 0; K < len_z; ++K) {
        _0EZX(pmlXNY,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
        _0EYX(pmlXNY,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
    }

    // Hzy
    // Hxy
    for (I = 0; I < plmLayerNMinus2; ++I) {
        for (J = 0; J < plmLayerNMinus2; ++J) {
            for (K = 0; K < len_z; ++K) {
                _EZY(pmlXNY,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
                _EXY(pmlXNY,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);

            }
        }
    }

    I = plmLayerNMinus1;
    J = plmLayerNMinus1;
    for (K = 0; K < len_z; ++K) {
        _0EZY(pmlXNY,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
        _0EXY(pmlXNY,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
    }
}
void FDTD::PML_UpdateEXY() {
    size_t I, J, K;

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < plmLayerNMinus1; ++I) {
        for (J = 0; J < plmLayerNMinus1; ++J) {
            for (K = 0; K < len_z -1; ++K) {
                _EYZ(pmlXY,I,J,K,1,pmlBHZ);
                _EXZ(pmlXY,I,J,K,1,pmlBHZ);
            }
        }
    }

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I < plmLayerNMinus2; ++I) {
        for (J = 0; J < plmLayerNMinus2; ++J) {
            for (K = 0; K < len_y; ++K) {
                _EZX(pmlXY,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
                _EYX(pmlXY,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);

            }
        }
    }

    I = plmLayerNMinus1;
    J = plmLayerNMinus1;
    for (K = 0; K < len_z; ++K) {
        _0EZX(pmlXY,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
        _0EYX(pmlXY,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
    }

    // Hzy
    // Hxy
    for (I = 0; I < plmLayerNMinus2; ++I) {
        for (J = 0; J < plmLayerNMinus2; ++J) {
            for (K = 0; K < len_z; ++K) {
                _EZY(pmlXY,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
                _EXY(pmlXY,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);

            }
        }
    }

    I = plmLayerNMinus1;
    J = plmLayerNMinus1;
    for (K = 0; K < len_z; ++K) {
        _0EZY(pmlXY,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
        _0EXY(pmlXY,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
    }
}

void FDTD::PML_UpdateEZY() {
    size_t I, J, K;

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I <len_x-1; ++I) {
        for (J = 0; J < plmLayerNMinus1; ++J) {
            for (K = 0; K < plmLayerNMinus1; ++K) {
                _EYX(pmlZY,I,J,K,1,pmlBHX);
                _EZX(pmlZY,I,J,K,1,pmlBHX);
            }
        }
    }

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerNMinus2; ++J) {
            for (K = 0; K < plmLayerNMinus2; ++K) {
                _EYZ(pmlZY,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXZ(pmlZY,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
            }
        }
    }

    J = plmLayerNMinus1;
    K = plmLayerNMinus1;
    for (I = 0; I < len_x; ++I) {
        _0EYZ(pmlZY,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
        _0EXZ(pmlZY,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    }

    // Hzy
    // Hxy
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerNMinus2; ++J) {
            for (K = 0; K < plmLayerNMinus2; ++K) {
                _EZY(pmlZY,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXY(pmlZY,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
            }
        }
    }
    J = plmLayerNMinus1;
    K = plmLayerNMinus1;
    for (I = 0; I < len_x; ++I) {
        _0EZY(pmlZY,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
        _0EXY(pmlZY,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    }
}
void FDTD::PML_UpdateEZYN() {
    size_t I, J, K;

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I <len_x-1; ++I) {
        for (J = 0; J < plmLayerNMinus1; ++J) {
            for (K = 0; K < plmLayerNMinus1; ++K) {
                _EYX(pmlZYN,I,J,K,1,pmlBHX);
                _EZX(pmlZYN,I,J,K,1,pmlBHX);
            }
        }
    }

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerNMinus2; ++J) {
            for (K = 0; K < plmLayerNMinus2; ++K) {
                _EYZ(pmlZYN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXZ(pmlZYN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
            }
        }
    }

    J = plmLayerNMinus1;
    K = plmLayerNMinus1;
    for (I = 0; I < len_x; ++I) {
        _0EYZ(pmlZYN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
        _0EXZ(pmlZYN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    }

    // Hzy
    // Hxy
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerNMinus2; ++J) {
            for (K = 0; K < plmLayerNMinus2; ++K) {
                _EZY(pmlZYN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXY(pmlZYN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
            }
        }
    }
    J = plmLayerNMinus1;
    K = plmLayerNMinus1;
    for (I = 0; I < len_x; ++I) {
        _0EZY(pmlZYN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
        _0EXY(pmlZYN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    }
}
void FDTD::PML_UpdateEZNY() {
    size_t I, J, K;

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I <len_x-1; ++I) {
        for (J = 0; J < plmLayerNMinus1; ++J) {
            for (K = 0; K < plmLayerNMinus1; ++K) {
                _EYX(pmlZNY,I,J,K,1,pmlBHX);
                _EZX(pmlZNY,I,J,K,1,pmlBHX);
            }
        }
    }

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerNMinus2; ++J) {
            for (K = 0; K < plmLayerNMinus2; ++K) {
                _EYZ(pmlZNY,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXZ(pmlZNY,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
            }
        }
    }

    J = plmLayerNMinus1;
    K = plmLayerNMinus1;
    for (I = 0; I < len_x; ++I) {
        _0EYZ(pmlZNY,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
        _0EXZ(pmlZNY,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    }

    // Hzy
    // Hxy
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerNMinus2; ++J) {
            for (K = 0; K < plmLayerNMinus2; ++K) {
                _EZY(pmlZNY,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXY(pmlZNY,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
            }
        }
    }
    J = plmLayerNMinus1;
    K = plmLayerNMinus1;
    for (I = 0; I < len_x; ++I) {
        _0EZY(pmlZNY,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
        _0EXY(pmlZNY,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    }
}
void FDTD::PML_UpdateEZNYN() {
size_t I, J, K;

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I <len_x-1; ++I) {
        for (J = 0; J < plmLayerNMinus1; ++J) {
            for (K = 0; K < plmLayerNMinus1; ++K) {
                _EYX(pmlZNYN,I,J,K,1,pmlBHX);
                _EZX(pmlZNYN,I,J,K,1,pmlBHX);
            }
        }
    }

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerNMinus2; ++J) {
            for (K = 0; K < plmLayerNMinus2; ++K) {
                _EYZ(pmlZNYN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXZ(pmlZNYN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
            }
        }
    }

    J = plmLayerNMinus1;
    K = plmLayerNMinus1;
    for (I = 0; I < len_x; ++I) {
        _0EYZ(pmlZNYN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
        _0EXZ(pmlZNYN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    }

    // Hzy
    // Hxy
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerNMinus2; ++J) {
            for (K = 0; K < plmLayerNMinus2; ++K) {
                _EZY(pmlZNYN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXY(pmlZNYN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
            }
        }
    }
    J = plmLayerNMinus1;
    K = plmLayerNMinus1;
    for (I = 0; I < len_x; ++I) {
        _0EZY(pmlZNYN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
        _0EXY(pmlZNYN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    }
}

void FDTD::PML_UpdateEZNXN() {
    size_t I, J, K;

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I <plmLayerNMinus2; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerNMinus2; ++K) {
                _EYX(pmlZNXN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EZX(pmlZNXN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
            }
        }
    }

    I = plmLayerNMinus1;
    K = plmLayerNMinus1;
    for (J = 0; J < len_y; ++J) {
        _0EYX(pmlZNXN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
        _0EZX(pmlZNXN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    }

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < plmLayerNMinus2; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerNMinus2; ++K) {
                _EYZ(pmlZNXN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXZ(pmlZNXN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
            }
        }
    }

    I = plmLayerNMinus1;
    K = plmLayerNMinus1;
    for (J = 0; J < len_y; ++J) {
        _0EYZ(pmlZNXN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
        _0EXZ(pmlZNXN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    }
    

    // Hzy
    // Hxy
    for (I = 0; I < plmLayerNMinus1; ++I) {
        for (J = 0; J < len_y-1; ++J) {
            for (K = 0; K < plmLayerNMinus1; ++K) {
                _EZY(pmlZNXN,I,J,K,1,pmlBHY);
                _EXY(pmlZNXN,I,J,K,1,pmlBHY);
            }
        }
    }
}
void FDTD::PML_UpdateEZNX() {
    size_t I, J, K;

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I <plmLayerNMinus2; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerNMinus2; ++K) {
                _EYX(pmlZNX,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EZX(pmlZNX,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
            }
        }
    }

    I = plmLayerNMinus1;
    K = plmLayerNMinus1;
    for (J = 0; J < len_y; ++J) {
        _0EYX(pmlZNX,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
        _0EZX(pmlZNX,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    }

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < plmLayerNMinus2; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerNMinus2; ++K) {
                _EYZ(pmlZNX,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXZ(pmlZNX,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
            }
        }
    }

    I = plmLayerNMinus1;
    K = plmLayerNMinus1;
    for (J = 0; J < len_y; ++J) {
        _0EYZ(pmlZNX,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
        _0EXZ(pmlZNX,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    }
    

    // Hzy
    // Hxy
    for (I = 0; I < plmLayerNMinus1; ++I) {
        for (J = 0; J < len_y-1; ++J) {
            for (K = 0; K < plmLayerNMinus1; ++K) {
                _EZY(pmlZNX,I,J,K,1,pmlBHY);
                _EXY(pmlZNX,I,J,K,1,pmlBHY);
            }
        }
    }
}
void FDTD::PML_UpdateEZXN() {
    size_t I, J, K;

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I <plmLayerNMinus2; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerNMinus2; ++K) {
                _EYX(pmlZXN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EZX(pmlZXN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
            }
        }
    }

    I = plmLayerNMinus1;
    K = plmLayerNMinus1;
    for (J = 0; J < len_y; ++J) {
        _0EYX(pmlZXN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
        _0EZX(pmlZXN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    }

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < plmLayerNMinus2; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerNMinus2; ++K) {
                _EYZ(pmlZXN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXZ(pmlZXN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
            }
        }
    }

    I = plmLayerNMinus1;
    K = plmLayerNMinus1;
    for (J = 0; J < len_y; ++J) {
        _0EYZ(pmlZXN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
        _0EXZ(pmlZXN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    }
    

    // Hzy
    // Hxy
    for (I = 0; I < plmLayerNMinus1; ++I) {
        for (J = 0; J < len_y-1; ++J) {
            for (K = 0; K < plmLayerNMinus1; ++K) {
                _EZY(pmlZXN,I,J,K,1,pmlBHY);
                _EXY(pmlZXN,I,J,K,1,pmlBHY);
            }
        }
    }
}
void FDTD::PML_UpdateEZX() {
    size_t I, J, K;

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I <plmLayerNMinus2; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerNMinus2; ++K) {
                _EYX(pmlZX,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EZX(pmlZX,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
            }
        }
    }

    I = plmLayerNMinus1;
    K = plmLayerNMinus1;
    for (J = 0; J < len_y; ++J) {
        _0EYX(pmlZX,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
        _0EZX(pmlZX,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    }

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < plmLayerNMinus2; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerNMinus2; ++K) {
                _EYZ(pmlZX,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXZ(pmlZX,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
            }
        }
    }

    I = plmLayerNMinus1;
    K = plmLayerNMinus1;
    for (J = 0; J < len_y; ++J) {
        _0EYZ(pmlZX,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
        _0EXZ(pmlZX,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
    }
    

    // Hzy
    // Hxy
    for (I = 0; I < plmLayerNMinus1; ++I) {
        for (J = 0; J < len_y-1; ++J) {
            for (K = 0; K < plmLayerNMinus1; ++K) {
                _EZY(pmlZX,I,J,K,1,pmlBHY);
                _EXY(pmlZX,I,J,K,1,pmlBHY);
            }
        }
    }
}
//подсчет граней
void FDTD::PML_UpdateEZ() {
    size_t I, J, K;

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I <len_x-1; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerNMinus1; ++K) {
                _EYX(pmlZ,I,J,K,1,pmlBHX);
                _EZX(pmlZ,I,J,K,1,pmlBHX);
            }
        }
    }

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerNMinus2; ++K) {
                _EYZ(pmlZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXZ(pmlZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
            }
        }
    }

    K = plmLayerNMinus1;
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < len_y; ++J) {
            _0EYZ(pmlZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
            _0EXZ(pmlZ,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
        }
    }

    // Hzy
    // Hxy
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < len_y-1; ++J) {
            for (K = 0; K < plmLayerNMinus1; ++K) {
                _EZY(pmlZ,I,J,K,1,pmlBHY);
                _EXY(pmlZ,I,J,K,1,pmlBHY);
            }
        }
    }
}
void FDTD::PML_UpdateEZN() {
size_t I, J, K;

    // Hyx
    //std::cout << "Hzx\n";  
    
    for (I = 0; I < len_x-1; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerNMinus1; ++K) {
                _EYX(pmlZN,I,J,K,1,pmlBHX);
                _EZX(pmlZN,I,J,K,1,pmlBHX);
            }
        }
    }

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerNMinus2; ++K) {
                _EYZ(pmlZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
                _EXZ(pmlZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
            }
        }
    }

    K = plmLayerNMinus1;
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < len_y; ++J) {
            _0EYZ(pmlZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
            _0EXZ(pmlZN,I,J,K,pmlSigmaStarE[K],pmlExpSigmaStarE[K]);
        }
    }

    // Hzy
    // Hxy
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < len_y-1; ++J) {
            for (K = 0; K < plmLayerNMinus1; ++K) {
                _EZY(pmlZN,I,J,K,1,pmlBHY);
                _EXY(pmlZN,I,J,K,1,pmlBHY);
            }
        }
    }
}
void FDTD::PML_UpdateEY() {
    size_t I, J, K;
    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I <len_x -1; ++I) {
        for (J = 0; J < plmLayerNMinus1; ++J) {
            for (K = 0; K < len_z; ++K) {
                _EYX(pmlY,I,J,K,1,pmlBHX);
                _EZX(pmlY,I,J,K,1,pmlBHX);
            }
        }
    }

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerNMinus1; ++J) {
            for (K = 0; K < len_z -1; ++K) {
                _EYZ(pmlY,I,J,K,1,pmlBHZ);
                _EXZ(pmlY,I,J,K,1,pmlBHZ);

            }
        }
    }

    // Hzy
    // Hxy
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerNMinus2; ++J) {
            for (K = 0; K < len_z; ++K) {
                _EZY(pmlY,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
                _EXY(pmlY,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
            }
        }
    }

    J = plmLayerNMinus1;
    for (I = 0; I < len_x; ++I) {
        for (K = 0; K < len_z; ++K) {
            _0EZY(pmlY,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
            _0EXY(pmlY,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
        }
    }
}
void FDTD::PML_UpdateEYN() {
    size_t I, J, K;

    // Hyx
    //std::cout << "Hzx\n";  

    for (I = 0; I <len_x - 1 ;++I) {
        for (J = 0; J < plmLayerNMinus1; ++J) {
            for (K = 0; K < len_y; ++K) {
                _EYX(pmlYN,I,J,K,1,pmlBHX);
                _EZX(pmlYN,I,J,K,1,pmlBHX);
            }
        }
    }

    // Hyz
    //std::cout << "Hxz\n";  

    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerNMinus1; ++J) {
            for (K = 0; K < len_z -1; ++K) {
                _EYZ(pmlYN,I,J,K,1,pmlBHZ);
                _EXZ(pmlYN,I,J,K,1,pmlBHZ);

            }
        }
    }

    // Hzy
    // Hxy

    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerNMinus2; ++J) {
            for (K = 0; K < len_y; ++K) {
                _EZY(pmlYN,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
                _EXY(pmlYN,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
            }
        }
    }

    J = plmLayerNMinus1;
    for (I = 0; I < len_x; ++I) {
        for (K = 0; K < len_z; ++K) {
            _0EZY(pmlYN,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
            _0EXY(pmlYN,I,J,K,pmlSigmaStarE[J],pmlExpSigmaStarE[J]);
        }
    }
}
void FDTD::PML_UpdateEX() {
    size_t I, J, K;

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < plmLayerNMinus1; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < len_z -1; ++K) {
                _EYZ(pmlX,I,J,K,1,pmlBHZ);
                _EXZ(pmlX,I,J,K,1,pmlBHZ);

            }
        }
    }

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I < plmLayerNMinus2; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < len_y; ++K) {
                _EZX(pmlX,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);
                _EYX(pmlX,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);

            }
        }
    }

    I = plmLayerNMinus1;
    for (J = 0; J < len_y; ++J) {
        for (K = 0; K < len_y; ++K) {
            _0EZX(pmlX,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);
            _0EYX(pmlX,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);
        }
    }

    // Hzy
    // Hxy
    for (I = 0; I < plmLayerNMinus1; ++I) {
        for (J = 0; J < len_y -1; ++J) {
            for (K = 0; K < len_y; ++K) {
                _EZY(pmlX,I,J,K,1,pmlBHY);
                _EXY(pmlX,I,J,K,1,pmlBHY);

            }
        }
    }
}
void FDTD::PML_UpdateEXN() {
    size_t I, J, K;

    // Hyz
    //std::cout << "Hxz\n";  
    for (I = 0; I < plmLayerNMinus1; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < len_z -1; ++K) {
                _EYZ(pmlXN,I,J,K,1,pmlBHZ);
                _EXZ(pmlXN,I,J,K,1,pmlBHZ);

            }
        }
    }

    // Hyx
    //std::cout << "Hzx\n";  
    for (I = 0; I < plmLayerNMinus2; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < len_y; ++K) {
                _EZX(pmlXN,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);
                _EYX(pmlXN,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);

            }
        }
    }

    I = plmLayerNMinus1;
    for (J = 0; J < len_y; ++J) {
        for (K = 0; K < len_y; ++K) {
            _0EZX(pmlXN,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);
            _0EYX(pmlXN,I,J,K,pmlSigmaStarE[I],pmlExpSigmaStarE[I]);
        }
    }

    // Hzy
    // Hxy
    for (I = 0; I < plmLayerNMinus1; ++I) {
        for (J = 0; J < len_y -1; ++J) {
            for (K = 0; K < len_y; ++K) {
                _EZY(pmlXN,I,J,K,1,pmlBHY);
                _EXY(pmlXN,I,J,K,1,pmlBHY);

            }
        }
    }
}
#endif