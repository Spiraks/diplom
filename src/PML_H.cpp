#include "FDTD.h"

void FDTD::PML_H(){
    //подсчет граней
    PML_UpdateHXN();
    PML_UpdateHX();
    PML_UpdateHYN();
    PML_UpdateHY();
    PML_UpdateHZN();
    PML_UpdateHZ();
    // //подсчет ребер
    // PML_UpdateHXNYN();
    // PML_UpdateHXYN();
    // PML_UpdateHXNY();
    // PML_UpdateHXY();

    // PML_UpdateHZNYN();
    // PML_UpdateHZYN();
    // PML_UpdateHZNY();
    // PML_UpdateHZY();

    // PML_UpdateHZNXN();
    // PML_UpdateHZXN();
    // PML_UpdateHZNX();
    // PML_UpdateHZX();

    // подсчет углов
    // PML_UpdateHXYZ();
    // PML_UpdateHXNYZ();
    // PML_UpdateHXYZN();
    // PML_UpdateHXYNZ();

    // PML_UpdateHXNYZN();
    // PML_UpdateHXYNZN();
    // PML_UpdateHXNYNZ();
    // PML_UpdateHXNYNZN();
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

void FDTD::PML_UpdateHZ() {
    size_t I, J, K;

    // Hyx
    // Hzx
    for (I = 0; I <len_x; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerNumberMinus1; ++K) {
                _HYX(pmlZ,I,J,K,1,pmlBHX);
                _HZX(pmlZ,I,J,K,1,pmlBHX);
            }
        }
    }

    // Hyz
    // Hxz
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerNumberMinus2; ++K) {
                _HYZ(pmlZ,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
                _HXZ(pmlZ,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
            }
        }
    }

    K = plmLayerNumberMinus1;
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < len_y; ++J) {
            _0HYZ(pmlZ,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
            _0HXZ(pmlZ,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
        }
    }

    // Hzy
    // Hxy
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < len_y-1; ++J) {
            for (K = 0; K < plmLayerNumberMinus1; ++K) {
                _HZY(pmlZ,I,J,K,1,pmlBHY);
                _HXY(pmlZ,I,J,K,1,pmlBHY);
            }
        }
    }
}
void FDTD::PML_UpdateHZN() {
size_t I, J, K;

    // Hyx
    // Hzx
    for (I = 0; I <len_x; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerNumberMinus1; ++K) {
                _HYX(pmlZN,I,J,K,1,pmlBHX);
                _HZX(pmlZN,I,J,K,1,pmlBHX);
            }
        }
    }

    // Hyz
    // Hxz
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < plmLayerNumberMinus2; ++K) {
                _HYZ(pmlZN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
                _HXZ(pmlZN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
            }
        }
    }

    K = plmLayerNumberMinus1;
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < len_y; ++J) {
            _0HYZ(pmlZN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
            _0HXZ(pmlZN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
        }
    }

    // Hzy
    // Hxy
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < len_y-1; ++J) {
            for (K = 0; K < plmLayerNumberMinus1; ++K) {
                _HZY(pmlZN,I,J,K,1,pmlBHY);
                _HXY(pmlZN,I,J,K,1,pmlBHY);
            }
        }
    }
}
void FDTD::PML_UpdateHY() {
    size_t I, J, K;
    // Hyx
    // Hzx
    for (I = 0; I <len_x -1; ++I) {
        for (J = 0; J < plmLayerNumberMinus1; ++J) {
            for (K = 0; K < len_z; ++K) {
                _HYX(pmlY,I,J,K,1,pmlBHX);
                _HZX(pmlY,I,J,K,1,pmlBHX);
            }
        }
    }

    // Hyz
    // Hxz
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerNumberMinus1; ++J) {
            for (K = 0; K < len_z -1; ++K) {
                _HYZ(pmlY,I,J,K,1,pmlBHZ);
                _HXZ(pmlY,I,J,K,1,pmlBHZ);

            }
        }
    }

    // Hzy
    // Hxy
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerNumberMinus2; ++J) {
            for (K = 0; K < len_z; ++K) {
                _HZY(pmlY,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
                _HXY(pmlY,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
            }
        }
    }

    J = plmLayerNumberMinus1;
    for (J = 0; J < len_y; ++J) {
        for (K = 0; K < len_z; ++K) {
            _0HZY(pmlY,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
            _0HXY(pmlY,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
        }
    }
}
void FDTD::PML_UpdateHYN() {
    size_t I, J, K;

    // Hyx
    // Hzx
    for (I = 0; I <len_x - 1 ;++I) {
        for (J = 0; J < plmLayerNumberMinus1; ++J) {
            for (K = 0; K < len_y; ++K) {
                _HYX(pmlYN,I,J,K,1,pmlBHX);
                _HZX(pmlYN,I,J,K,1,pmlBHX);
            }
        }
    }

    // Hyz
    // Hxz
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerNumberMinus1; ++J) {
            for (K = 0; K < len_z -1; ++K) {
                _HYZ(pmlYN,I,J,K,1,pmlBHZ);
                _HXZ(pmlYN,I,J,K,1,pmlBHZ);

            }
        }
    }

    // Hzy
    // Hxy
    for (I = 0; I < len_x; ++I) {
        for (J = 0; J < plmLayerNumberMinus2; ++J) {
            for (K = 0; K < len_y; ++K) {
                _HZY(pmlYN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
                _HXY(pmlYN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
            }
        }
    }

    J = plmLayerNumberMinus1;
    for (J = 0; J < len_y; ++J) {
        for (K = 0; K < len_y; ++K) {
            _0HZY(pmlYN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
            _0HXY(pmlYN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
        }
    }
}
void FDTD::PML_UpdateHX() {
    size_t I, J, K;

    // Hzx
    // Hyx

    for (I = 0; I <= plmLayerNumberMinus2; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < len_y; ++K) {
                _HZX(pmlX,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);\
                _HYX(pmlX,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
            }
        }
    }

    I = plmLayerNumberMinus1;
    for (J = 0; J < len_y; ++J) {
        for (K = 0; K < len_y; ++K) {
            _0HYX(pmlX,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
            _0HZX(pmlX,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
        }   
    }

    // Hzy
    // Hxy

    for (I = 0; I <= plmLayerNumberMinus1; ++I) {
        for (J = 0; J < len_y -1; ++J) {
            for (K = 0; K < len_y; ++K) {
                _HZY(pmlX,I,J,K,1,pmlBHY);
                _HXY(pmlX,I,J,K,1,pmlBHY);
            }
        }
    }

    // Hyz
    // Hxz
    for (I = 0; I <= plmLayerNumberMinus1; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < len_z -1; ++K) {
                _HYZ(pmlX,I,J,K,1,pmlBHZ);
                _HXZ(pmlX,I,J,K,1,pmlBHZ);
            }
        }
    }
}
void FDTD::PML_UpdateHXN() {
    size_t I, J, K;
    // Здесь должен быть код для заполнения pmlSigmaStarH соответствующими значениями
    for (I = 0; I <= plmLayerNumberMinus2; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < len_y; ++K) {
                _HZX(pmlXN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);\
                _HYX(pmlXN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
            }
        }
    }

    I = plmLayerNumberMinus1;
    for (J = 0; J < len_y; ++J) {
        for (K = 0; K < len_y; ++K) {
            _0HYX(pmlXN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
            _0HZX(pmlXN,I,J,K,pmlSigmaStarH[I],pmlExpSigmaStarX[I]);
        }   
    }

    // Hzy
    // Hxy

    for (I = 0; I <= plmLayerNumberMinus1; ++I) {
        for (J = 0; J < len_y -1; ++J) {
            for (K = 0; K < len_y; ++K) {
                _HZY(pmlXN,I,J,K,1,pmlBHY);
                _HXY(pmlXN,I,J,K,1,pmlBHY);
            }
        }
    }

    // Hyz
    // Hxz
    for (I = 0; I <= plmLayerNumberMinus1; ++I) {
        for (J = 0; J < len_y; ++J) {
            for (K = 0; K < len_z -1; ++K) {
                _HYZ(pmlXN,I,J,K,1,pmlBHZ);
                _HXZ(pmlXN,I,J,K,1,pmlBHZ);
            }
        }
    }
}
