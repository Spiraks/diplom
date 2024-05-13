#include "FDTD.h"

void FDTD::PML_E(){
    //подсчет граней
    PML_UpdateEXN();
    PML_UpdateEX();
    PML_UpdateEYN();
    PML_UpdateEY();
    PML_UpdateEZN();
    PML_UpdateEZ();
    // //подсчет ребер
    // PML_UpdateEXNYN();
    // PML_UpdateEXYN();
    // PML_UpdateEXNY();
    // PML_UpdateEXY();

    // PML_UpdateEZNYN();
    // PML_UpdateEZYN();
    // PML_UpdateEZNY();
    // PML_UpdateEZY();

    // PML_UpdateEZNXN();
    // PML_UpdateEZXN();
    // PML_UpdateEZNX();
    // PML_UpdateEZX();

    // // подсчет углов
    // PML_UpdateEXYZ();
    // PML_UpdateEXNYZ();
    // PML_UpdateEXYZN();
    // PML_UpdateEXYNZ();

    // PML_UpdateEXNYZN();
    // PML_UpdateEXYNZN();
    // PML_UpdateEXNYNZ();

    // PML_UpdateEXNYNZN();
}
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


void FDTD::PML_UpdateEZ() {
    size_t I, J, K;

    // Hyx
    // Hzx
    for (I = 0; I <len_x; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerNumberMinus1; ++K) {
                _EYX(pmlZ,I,J,K,1,pmlBHX);
                _EZX(pmlZ,I,J,K,1,pmlBHX);
            }
        }
    }

    // Hyz
    // Hxz
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerNumberMinus2; ++K) {
                _EYZ(pmlZ,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
                _EXZ(pmlZ,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
            }
        }
    }

    K = plmLayerNumberMinus1;
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < len_y; ++J) {
            _0EYZ(pmlZ,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
            _0EXZ(pmlZ,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
        }
    }

    // Hzy
    // Hxy
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < len_y-1; ++J) {
            for (K = 0; K < plmLayerNumberMinus1; ++K) {
                _EZY(pmlZ,I,J,K,1,pmlBHY);
                _EXY(pmlZ,I,J,K,1,pmlBHY);
            }
        }
    }
}
void FDTD::PML_UpdateEZN() {
size_t I, J, K;

    // Hyx
    // Hzx
    for (I = 0; I <len_x; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerNumberMinus1; ++K) {
                _EYX(pmlZN,I,J,K,1,pmlBHX);
                _EZX(pmlZN,I,J,K,1,pmlBHX);
            }
        }
    }

    // Hyz
    // Hxz
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerNumberMinus2; ++K) {
                _EYZ(pmlZN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
                _EXZ(pmlZN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
            }
        }
    }

    K = plmLayerNumberMinus1;
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < len_y; ++J) {
            _0EYZ(pmlZN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
            _0EXZ(pmlZN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
        }
    }

    // Hzy
    // Hxy
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < len_y-1; ++J) {
            for (K = 0; K < plmLayerNumberMinus1; ++K) {
                _EZY(pmlZN,I,J,K,1,pmlBHY);
                _EXY(pmlZN,I,J,K,1,pmlBHY);
            }
        }
    }
}
void FDTD::PML_UpdateEY() {
    size_t I, J, K;
    // Hyx
    // Hzx
    for (I = 0; I <len_x -1; ++I) {
        for (J = 0; J < plmLayerNumberMinus1; ++J) {
            for (K = 0; K < len_z; ++K) {
                _EYX(pmlY,I,J,K,1,pmlBHX);
                _EZX(pmlY,I,J,K,1,pmlBHX);
            }
        }
    }

    // Hyz
    // Hxz
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerNumberMinus1; ++J) {
            for (K = 0; K < len_z -1; ++K) {
                _EYZ(pmlY,I,J,K,1,pmlBHZ);
                _EXZ(pmlY,I,J,K,1,pmlBHZ);

            }
        }
    }

    // Hzy
    // Hxy
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerNumberMinus2; ++J) {
            for (K = 0; K < len_z; ++K) {
                _EZY(pmlY,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
                _EXY(pmlY,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
            }
        }
    }

    J = plmLayerNumberMinus1;
    for (J = 0; J < len_y; ++J) {
        for (K = 0; K < len_z; ++K) {
            _0EZY(pmlY,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
            _0EXY(pmlY,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
        }
    }
}
void FDTD::PML_UpdateEYN() {
    size_t I, J, K;

    // Hyx
    // Hzx
    for (I = 0; I <len_x - 1 ;++I) {
        for (J = 0; J < plmLayerNumberMinus1; ++J) {
            for (K = 0; K < len_y; ++K) {
                _EYX(pmlYN,I,J,K,1,pmlBHX);
                _EZX(pmlYN,I,J,K,1,pmlBHX);
            }
        }
    }

    // Hyz
    // Hxz
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerNumberMinus1; ++J) {
            for (K = 0; K < len_z -1; ++K) {
                _EYZ(pmlYN,I,J,K,1,pmlBHZ);
                _EXZ(pmlYN,I,J,K,1,pmlBHZ);

            }
        }
    }

    // Hzy
    // Hxy
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerNumberMinus2; ++J) {
            for (K = 0; K < len_y; ++K) {
                _EZY(pmlYN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
                _EXY(pmlYN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
            }
        }
    }

    J = plmLayerNumberMinus1;
    for (J = 0; J < len_y; ++J) {
        for (K = 0; K < len_y; ++K) {
            _0EZY(pmlYN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
            _0EXY(pmlYN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
        }
    }
}
void FDTD::PML_UpdateEX() {
    size_t I, J, K;

    // Hyz
    // Hxz
    for (I = 0; I <= plmLayerNumberMinus1; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < len_z -1; ++K) {
                _EYZ(pmlX,I,J,K,1,pmlBHZ);
                _EXZ(pmlX,I,J,K,1,pmlBHZ);

            }
        }
    }

    // Hyx
    // Hzx
    for (I = 0; I <= plmLayerNumberMinus2; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < len_y; ++K) {
                _EZX(pmlX,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
                _EYX(pmlX,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);

            }
        }
    }

    I = plmLayerNumberMinus1;
    for (J = 0; J < len_y; ++J) {
        for (K = 0; K < len_y; ++K) {
            _0EZX(pmlX,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
            _0EYX(pmlX,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
        }
    }

    // Hzy
    // Hxy
    for (I = 0; I <= plmLayerNumberMinus1; ++I) {
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
    // Hxz
    for (I = 0; I <= plmLayerNumberMinus1; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < len_z -1; ++K) {
                _EYZ(pmlXN,I,J,K,1,pmlBHZ);
                _EXZ(pmlXN,I,J,K,1,pmlBHZ);

            }
        }
    }

    // Hyx
    // Hzx
    for (I = 0; I <= plmLayerNumberMinus2; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < len_y; ++K) {
                _EZX(pmlXN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
                _EYX(pmlXN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);

            }
        }
    }

    I = plmLayerNumberMinus1;
    for (J = 0; J < len_y; ++J) {
        for (K = 0; K < len_y; ++K) {
            _0EZX(pmlXN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
            _0EYX(pmlXN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
        }
    }

    // Hzy
    // Hxy
    for (I = 0; I <= plmLayerNumberMinus1; ++I) {
        for (J = 0; J < len_y -1; ++J) {
            for (K = 0; K < len_y; ++K) {
                _EZY(pmlXN,I,J,K,1,pmlBHY);
                _EXY(pmlXN,I,J,K,1,pmlBHY);

            }
        }
    }
}
