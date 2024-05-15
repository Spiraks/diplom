#include <iostream>
#include <fstream>
#include <unistd.h>
#include <math.h>
#include <vector>

#define X 0
#define Y 1
#define Z 2
// #define MUR
// #define PML

#define VECTOR std::vector<float>
#define VECTOR_2D std::vector<VECTOR>
#define VECTOR_3D std::vector<VECTOR_2D>
#define VECTOR_4D std::vector<VECTOR_3D>

#define V_PML std::vector<PMLData>
#define V_PML_2D std::vector<V_PML>
#define V_PML_3D std::vector<V_PML_2D>

// Константы и параметры для моделирования
const float c = 2.997925010e8;   //Скорость света (м/с)
const float eps0 = 8.854e-12;    //Вакуумная диэлектрическая проницаемость (Ф/м)
const float mu0 = 4 * M_PI * 1.0e-7; // Вакуумная магнитная проницаемость (Гн/м)
const float dt = 1.0e-12;            // Временной шаг (с)
const float dstep = 0.3e-2;
const float dx = dstep; // Пространственный шаг по x (м)
const float dy = dstep; // Пространственный шаг по y (м)
const float dz = dstep; // Пространственный шаг по z (м)
const float sigma = 0;  // Коэффициент проводимости для потерь

// Для граничных условий мура
const float Mur_Tx = c * dt / dx;
const float Mur_Ty = c * dt / dy;
const float Mur_Tz = c * dt / dz;

struct PMLData {
    float Hyx, Hyz, Hzx, Hzy, Hxy, Hxz;
    float Exy, Exz, Eyx, Eyz, Ezx, Ezy;
};

struct Mesh3D
{
    int len_x, len_y, len_z;
};

class Field
{
public:
    Field() {}
    Field(Mesh3D &mesh, float value = 0)
    {
        len_x = mesh.len_x;
        len_y = mesh.len_y;
        len_z = mesh.len_z;
        this->_value = value;
        nodes = VECTOR_4D(len_x, VECTOR_3D(len_y, VECTOR_2D(len_z, VECTOR(3, _value))));
        nodes_1 = VECTOR_4D(len_x, VECTOR_3D(len_y, VECTOR_2D(len_z, VECTOR(3, _value))));
    }

    float gt(size_t x, size_t y, size_t z, size_t comp);
    void setNode(size_t x, size_t y, size_t z, size_t comp, float value);
    float gtO(size_t x, size_t y, size_t z, size_t comp);

    void updateE(Field &H, Field &J, Field &c1, Field &c2);
    void updateH(Field &E, Field &c0);

    void saveExToFileZ(const std::string &filename, int comp);
    void saveExToFileX(const std::string &filename, int comp);
    void saveExToFileY(const std::string &filename, int comp);
    void saveToFile(const std::string &filename);

    void saveGraphZ(const std::string &filename, int comp);
    void saveRangeY();

    void copyNodes();

    size_t getSizeX()
    {
        return len_x;
    }
    size_t getSizeY()
    {
        return len_y;
    }
    size_t getSizeZ()
    {
        return len_z;
    }
    float getMinX() const
    {
        return min_x;
    }
    float getMinY() const
    {
        return min_y;
    }
    float getMinZ() const
    {
        return min_z;
    }
    float getMaxX() const
    {
        return max_x;
    }
    float getMaxY() const
    {
        return max_y;
    }
    float getMaxZ() const
    {
        return max_z;
    }

private:
    int len_x = 0;
    int len_y = 0;
    int len_z = 0;
    float _value = 0;
    float min_x = 0;
    float min_y = 0;
    float min_z = 0;
    float max_x = 0;
    float max_y = 0;
    float max_z = 0;
    VECTOR_4D nodes;
    //for mur
    VECTOR_4D nodes_1;
};

class FDTD
{
private:
    Field c0;
    Field c1;
    Field c2;
    float _time = 0;
    int len_x = 0;
    int len_y = 0;
    int len_z = 0;
        #ifdef PML
            #ifdef MUR
                #error
            #else
        //for pml
    float pmlBHX = 1;
    float pmlBHZ = 1;
    float pmlBHY = 1;
    int plmLayerN = 15;
    size_t plmLayerNMinus1 = plmLayerN - 1;
    size_t plmLayerNMinus2 = plmLayerN - 2;
    std::vector<float> pmlSigmaStarH;
    std::vector<float> pmlSigmaStarE;
    std::vector<float> pmlExpSigmaStarH; 
    std::vector<float> pmlExpSigmaStarE;

    // по граням
    V_PML_3D pmlXN;
    V_PML_3D pmlX ;
    V_PML_3D pmlYN;
    V_PML_3D pmlY ;
    V_PML_3D pmlZN;
    V_PML_3D pmlZ ;

    // по ребрам
    V_PML_3D pmlXNYN;
    V_PML_3D pmlXYN;
    V_PML_3D pmlXNY;
    V_PML_3D pmlXY;

    V_PML_3D pmlZNYN;
    V_PML_3D pmlZYN;
    V_PML_3D pmlZNY;
    V_PML_3D pmlZY;

    V_PML_3D pmlZNXN;
    V_PML_3D pmlZXN;
    V_PML_3D pmlZNX;
    V_PML_3D pmlZX;

    // по углам
    V_PML_3D pmlXYZ;
    V_PML_3D pmlXNYZ;
    V_PML_3D pmlXYZN;
    V_PML_3D pmlXYNZ;
    V_PML_3D pmlXNYZN;
    V_PML_3D pmlXYNZN;
    V_PML_3D pmlXNYNZ;
    V_PML_3D pmlXNYNZN;
                #endif    
        #endif

public:



    Field E;
    Field H;
    Field Sigm; // 2 -14
    Field J;
    Field Eps; // 1.00057
    Field Mu;  // 0.999991

    FDTD(Mesh3D &mesh)
    {
        E = Field(mesh);
        H = Field(mesh);
        J = Field(mesh);
        Sigm = Field(mesh, 2.0e-5); // 2 -14
        Eps = Field(mesh, 1.00057);  // 1.00057
        Mu = Field(mesh, 0.999991);  // 0.999991
        c0 = Field(mesh);
        c1 = Field(mesh);
        c2 = Field(mesh);
        len_x = mesh.len_x;
        len_y = mesh.len_y;
        len_z = mesh.len_z;
        std::cout << "SIZE FIELD: \n" << "len_x = " << len_x << "\n"<<"len_y = " << len_y << "\n" <<"len_z = " << len_z << "\n";

        
        initCoeffi();
        #ifdef PML
            #ifdef MUR
                #error
            #else
                pmlSigmaStarE = VECTOR(plmLayerN);
                pmlSigmaStarH = VECTOR(plmLayerN);
                pmlExpSigmaStarE = VECTOR(plmLayerN,1);
                pmlExpSigmaStarH = VECTOR(plmLayerN,1);

                pmlXN = V_PML_3D(plmLayerN,V_PML_2D(len_y, V_PML(len_z)));
                pmlX  = V_PML_3D(plmLayerN,V_PML_2D(len_y, V_PML(len_z)));
                pmlYN = V_PML_3D(len_x,V_PML_2D(plmLayerN, V_PML(len_z)));
                pmlY  = V_PML_3D(len_x,V_PML_2D(plmLayerN, V_PML(len_z)));
                pmlZN = V_PML_3D(len_x,V_PML_2D(len_y, V_PML(plmLayerN)));
                pmlZ  = V_PML_3D(len_x,V_PML_2D(len_y, V_PML(plmLayerN)));

                pmlXNYN =  V_PML_3D(plmLayerN,V_PML_2D(plmLayerN, V_PML(len_z)));
                pmlXYN =   V_PML_3D(plmLayerN,V_PML_2D(plmLayerN, V_PML(len_z)));
                pmlXNY =   V_PML_3D(plmLayerN,V_PML_2D(plmLayerN, V_PML(len_z)));
                pmlXY =    V_PML_3D(plmLayerN,V_PML_2D(plmLayerN, V_PML(len_z)));

                pmlZNYN =  V_PML_3D(len_x,V_PML_2D(plmLayerN, V_PML(plmLayerN)));
                pmlZYN =   V_PML_3D(len_x,V_PML_2D(plmLayerN, V_PML(plmLayerN)));
                pmlZNY =   V_PML_3D(len_x,V_PML_2D(plmLayerN, V_PML(plmLayerN)));
                pmlZY =    V_PML_3D(len_x,V_PML_2D(plmLayerN, V_PML(plmLayerN)));

                pmlZNXN =  V_PML_3D(plmLayerN,V_PML_2D(len_y, V_PML(plmLayerN)));
                pmlZXN =   V_PML_3D(plmLayerN,V_PML_2D(len_y, V_PML(plmLayerN)));
                pmlZNX =   V_PML_3D(plmLayerN,V_PML_2D(len_y, V_PML(plmLayerN)));
                pmlZX =    V_PML_3D(plmLayerN,V_PML_2D(len_y, V_PML(plmLayerN)));

                pmlXYZ  =  V_PML_3D(plmLayerN,V_PML_2D(plmLayerN, V_PML(plmLayerN)));
                pmlXNYZ  =  V_PML_3D(plmLayerN,V_PML_2D(plmLayerN, V_PML(plmLayerN)));
                pmlXYZN  =  V_PML_3D(plmLayerN,V_PML_2D(plmLayerN, V_PML(plmLayerN)));
                pmlXYNZ  =  V_PML_3D(plmLayerN,V_PML_2D(plmLayerN, V_PML(plmLayerN)));
                pmlXNYZN  =  V_PML_3D(plmLayerN,V_PML_2D(plmLayerN, V_PML(plmLayerN)));
                pmlXYNZN  =  V_PML_3D(plmLayerN,V_PML_2D(plmLayerN, V_PML(plmLayerN)));
                pmlXNYNZ  =  V_PML_3D(plmLayerN,V_PML_2D(plmLayerN, V_PML(plmLayerN)));
                pmlXNYNZN  =  V_PML_3D(plmLayerN,V_PML_2D(plmLayerN, V_PML(plmLayerN)));
                InitPML();
            #endif    
        #endif
    }
    void update(float time);
    void initCoeffi();

    void Mur();
    void GetBorderValuesMur();

    void Pml();
    void PML_E();
    void PML_H();
    void InitPML();
    float GetSigma(float D, float dStep, float pmlG, float pmlR,int N);
    void GetBorderValuesPML();
    void SetBorderValuesPML();



    //подсчет граней
    void PML_UpdateEXN();
    void PML_UpdateEX();
    void PML_UpdateEYN();
    void PML_UpdateEY();
    void PML_UpdateEZN();
    void PML_UpdateEZ();
    //подсчет ребер
    void PML_UpdateEXNYN();
    void PML_UpdateEXYN();
    void PML_UpdateEXNY();
    void PML_UpdateEXY();

    void PML_UpdateEZNYN();
    void PML_UpdateEZYN();
    void PML_UpdateEZNY();
    void PML_UpdateEZY();

    void PML_UpdateEZNXN();
    void PML_UpdateEZXN();
    void PML_UpdateEZNX();
    void PML_UpdateEZX();

    // подсчет углов
    void PML_UpdateEXYZ();
    void PML_UpdateEXNYZ();
    void PML_UpdateEXYZN();
    void PML_UpdateEXYNZ();

    void PML_UpdateEXNYZN();
    void PML_UpdateEXYNZN();
    void PML_UpdateEXNYNZ();
    void PML_UpdateEXNYNZN();

    //подсчет граней
    void PML_UpdateHXN();
    void PML_UpdateHX();
    void PML_UpdateHYN();
    void PML_UpdateHY();
    void PML_UpdateHZN();
    void PML_UpdateHZ();
    //подсчет ребер
    void PML_UpdateHXNYN();
    void PML_UpdateHXYN();
    void PML_UpdateHXNY();
    void PML_UpdateHXY();

    void PML_UpdateHZNYN();
    void PML_UpdateHZYN();
    void PML_UpdateHZNY();
    void PML_UpdateHZY();

    void PML_UpdateHZNXN();
    void PML_UpdateHZXN();
    void PML_UpdateHZNX();
    void PML_UpdateHZX();

    // подсчет углов
    void PML_UpdateHXYZ();
    void PML_UpdateHXNYZ();
    void PML_UpdateHXYZN();
    void PML_UpdateHXYNZ();

    void PML_UpdateHXNYZN();
    void PML_UpdateHXYNZN();
    void PML_UpdateHXNYNZ();
    void PML_UpdateHXNYNZN();
};